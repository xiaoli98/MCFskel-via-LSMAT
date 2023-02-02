//
// Created by boy on 04/12/22.
//

#ifndef MCFSKET_COLLAPSER_H
#define MCFSKET_COLLAPSER_H

#define COLLAPESER_DEBUG 0
#define MESH_DEBUG 1

#include <map>
#include <vcg/complex/algorithms/edge_collapse.h>
#include "utils.hpp"

class Collapser{
    typedef MyMesh::PerVertexAttributeHandle<double> omega;
    typedef vector<tuple<MyMesh::VertexType*, MyMesh::VertexType*, double>> laplacian_triple;
    typedef tri::BasicVertexPair<MyMesh::VertexType> VertexPair;

private:
    MyMesh *m;
    double omega_L_0;
    omega omega_L;
    double omega_H_0;
    omega omega_H;
    double omega_M_0;
    omega omega_M;
    double zero_TH;

    Eigen::SparseMatrix<double> LHS;
    Eigen::MatrixXd RHS;
    Eigen::MatrixXd X;

    int nrows, ncols;

    MyMesh::PerVertexAttributeHandle<MyMesh::CoordType> vert_mat;
    MyMesh::PerVertexAttributeHandle<bool> isFixed;
    MyMesh::PerVertexAttributeHandle<bool> isSplitted;
    MyMesh::PerMeshAttributeHandle<laplacian_triple> laplacian;
    MyMesh::PerMeshAttributeHandle<double> edge_threshold;
    tri::EdgeCollapser<MyMesh, VertexPair> collapser;

public:
    Collapser(MyMesh* m){
        this->m = m;
        omega_L = tri::Allocator<MyMesh>::GetPerVertexAttribute<double>(*m, string("omega_L"));
        omega_L_0 = 1;
        omega_H = tri::Allocator<MyMesh>::GetPerVertexAttribute<double>(*m, string("omega_H"));
        omega_H_0 = 20;
        omega_M = tri::Allocator<MyMesh>::GetPerVertexAttribute<double>(*m, string("omega_M"));
        omega_M_0 = 40;
        for (auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            omega_L[vi] = omega_L_0;
            omega_H[vi] = omega_H_0;
            omega_M[vi] = omega_M_0;
        }
        zero_TH = 1e-7;

        vert_mat = tri::Allocator<MyMesh>::GetPerVertexAttribute<MyMesh::CoordType>(*m, string("vert_mat"));
        isFixed = tri::Allocator<MyMesh>::GetPerVertexAttribute<bool>(*m, string("isFixed"));
        for(auto vi = m->vert.begin(); vi!=m->vert.end(); vi++)
            isFixed[vi] = false;
        isSplitted = tri::Allocator<MyMesh>::GetPerVertexAttribute<bool>(*m, string("isSplitted"));
        edge_threshold = tri::Allocator<MyMesh>::GetPerMeshAttribute<double>(*m, string("edge_threshold"));
        edge_threshold() = 0.002*m->bbox.Diag();
    }

    void compute(){
        createLHS();
        createRHS();
        solveLS();
        MyMesh::PerVertexAttributeHandle<MyMesh::CoordType> sol = tri::Allocator<MyMesh>::GetPerVertexAttribute<MyMesh::CoordType>(*m, string("solution"));

        for (auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
            vi->P() = sol[vi];
        }
//        update_omega();
        update_topology();
        detect_degeneracies();
        tri::io::ExporterOFF<MyMesh>::Save(*m, "prova.off", tri::io::Mask::IOM_FACECOLOR);
    }
    void createLHS(){
        laplacian = tri::Allocator<MyMesh>::GetPerMeshAttribute<laplacian_triple>(*m, string("laplacian"));
        nrows = 3 * m->VN();
        ncols = m->VN();
        LHS.setZero();
        LHS.resize(nrows, ncols);
        RHS = Eigen::MatrixXd::Zero(nrows, 3);
        X = Eigen::MatrixXd::Zero(ncols, 3);

        vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(ncols*9);

        MyMesh::VertexPointer v0,v1;
        double w;
        for (auto l = laplacian().begin(); l != laplacian().end(); l++){
            w = get<2>(*l);
            v0 = get<0>(*l);
            v1 = get<1>(*l);
            if(v0 != v1)
                triplets.emplace_back(Eigen::Triplet<double>(tri::Index(*m, v0), tri::Index(*m, v1), w * omega_L[v0]));
            else
                triplets.emplace_back(Eigen::Triplet<double>(tri::Index(*m, v0), tri::Index(*m, v0), w));
        }

        for(auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
            triplets.emplace_back(Eigen::Triplet<double>(tri::Index(*m, *vi)+ncols, tri::Index(*m, *vi), omega_H[*vi]));
        }
        for(auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
            triplets.emplace_back(Eigen::Triplet<double>(tri::Index(*m, *vi)+2*ncols, tri::Index(*m, *vi), omega_M[*vi]));
        }
#if COLLAPESER_DEBUG
        cout << triplets.size()<<endl;
        cout << "LHS rows:" << LHS.rows()<<endl;
        cout << "LHS cols:" << LHS.cols()<<endl;
        for (auto it = triplets.begin(); it != triplets.end(); it++){
            if (it->row()<0 || it->row()>=LHS.rows() || it->col()<0 || it->col()>=LHS.cols())
                printf("col: %d\trow:%d\t value: %f\n", it->col(), it->row(), it->value());
        }

#endif
        LHS.setFromTriplets(triplets.begin(), triplets.end());
    }

    void createRHS(){
        for(auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
            MyMesh::CoordType c = (float)omega_H[vi] * vi->P();
            RHS.row(ncols + tri::Index(*m, *vi)) = Eigen::Vector3d(c.X(), c.Y(), c.Z());
        }

        for(auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
            MyMesh::CoordType c = (float) omega_M[vi] * vert_mat[vi];
            RHS.row(2*ncols + tri::Index(*m, *vi)) = Eigen::Vector3d(c.X(), c.Y(), c.Z());
        }
    }

    void solveLS(){
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(LHS.transpose() * LHS);

        X.col(0) = solver.solve(LHS.transpose() * RHS.col(0));
        X.col(1) = solver.solve(LHS.transpose() * RHS.col(1));
        X.col(2) = solver.solve(LHS.transpose() * RHS.col(2));

//        printf("LHS:\n");
//        cout << Eigen::MatrixXd(LHS)<<endl;

        if (!std::isfinite(X.norm())){
            throw EIGEN_EXCEPTIONS("Problem with linear least square solution.");
        }

        MyMesh::PerVertexAttributeHandle<MyMesh::CoordType> sol = tri::Allocator<MyMesh>::GetPerVertexAttribute<MyMesh::CoordType>(*m, string("solution"));
        for(auto vi=m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
            Eigen::Vector3d p = X.row(tri::Index(*m, *vi));
            sol[vi] = MyMesh::CoordType(p[0], p[1], p[2]);
        }
//        for(auto vi=m->vert.begin(); vi != m->vert.end(); vi++){
//            printf("sol[vi]: ");
//            PRINTP(sol[vi])
//        }
    }

    void update_omega(){
        printf("update_omega\n");
        for (auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
            if(isFixed[vi]){
                printf("setting omega_L to zero: %p\n", vi.base());
                omega_L[vi] = 0;
                omega_H[vi] = 1.0/zero_TH;
                omega_M[vi] = 0;
                continue;
            }

            omega_L[vi] = omega_L_0;
            omega_H[vi] = omega_H_0;
            omega_M[vi] = omega_M_0;

            if(isSplitted[vi]){
                printf("setting splitted vert: %p\n", vi.base());
                omega_L[vi] = omega_L_0;
                omega_H[vi] = omega_H_0;
                omega_M[vi] = 0;
            }
        }
    }

    typedef tuple<MyMesh::VertexType*,MyMesh::VertexType*,MyMesh::VertexType*> triVert;
    void update_topology(){
//        edge_split();
//        tri::io::ExporterOFF<MyMesh>::Save(*m, "after split_1.off", tri::io::Mask::IOM_FACECOLOR);
        int total_collapse = 0;
        int n_collapse = 0;
        do{
            n_collapse = edge_collapse();
            total_collapse += n_collapse;
        }while(n_collapse > 0);
//        tri::io::ExporterOFF<MyMesh>::Save(*m, "after collapse.off", tri::io::Mask::IOM_FACECOLOR);
        int total_split = 0;
        int n_split = 0;
        do{
            n_split = edge_split();
            total_split += n_split;
        }while(n_split > 0);
//        tri::io::ExporterOFF<MyMesh>::Save(*m, "after split_2.off", tri::io::Mask::IOM_FACECOLOR);
//        edge_split();
//        tri::io::ExporterOFF<MyMesh>::Save(*m, "after split_3.off", tri::io::Mask::IOM_FACECOLOR);

//        cout << "topology updated, face and vertex cleaned"<<endl;
    }

    int edge_collapse(){
        printf("edge_collapse: ");
        int count=0;
//        printf("edge_threshold: %f\n", edge_threshold());
        MyMesh::VertexType *V0, *V1;
        double edge_length;
        tri::UpdateFlags<MyMesh>::VertexClearV(*m);

        vector<MyMesh::VertexPointer> adjVertices;
        for(auto vit = m->vert.begin(); vit != m->vert.end(); vit++) {
            if(!vit->IsD()) {
                face::VVStarVF<MyMesh::FaceType>(vit.base(), adjVertices);
                for (auto adjv = adjVertices.begin(); adjv != adjVertices.end(); adjv++) {
                    if ((*adjv)->IsD()) continue;
                    if(isFixed[vit.base()] && isFixed[*adjv.base()]) continue;
                    MyMesh::FacePointer fp = vit->VFp();
                    V0 = vit.base();
                    V1 = (*adjv);
                    VertexPair vpair(V0, V1);
                    if(collapser.LinkConditions(vpair)){
                        edge_length = Distance(V0->P(), V1->P());
                        if (edge_length < edge_threshold()) {
//                            cout << "\tV0: "<<V0->P().X() <<"\t"<< V0->P().Y() <<"\t"<< V0->P().Z();
//                            cout << "\tV1: "<<V1->P().X() <<"\t"<< V1->P().Y() <<"\t"<< V1->P().Z();
//                            printf("\t edge_length: %f\t threshold: %f\n", edge_length, edge_threshold());
                            V1->P() = (V1->P() + V0->P()) / 2;

                            auto mat0 = vert_mat[V0];
                            auto mat1 = vert_mat[V1];
                            //find the closest medial axis to V1
                            double dist1 = Distance(mat0, V1->P());
                            double dist2 = Distance(mat1, V1->P());

                            if (dist1 < dist2)
                                vert_mat[V1] = vert_mat[V0];
                            else
                                vert_mat[V1] = vert_mat[V1];

                            collapser.Do(*m, vpair, V1->P());
                            count++;
                            break;
                        }
                    }
                }
            }
        }


//        tri::UpdateTopology<MyMesh>::FaceFace(*m);
        printf("collapsed %d times\n", count);
        tri::Allocator<MyMesh>::CompactFaceVector(*m);
        tri::Allocator<MyMesh>::CompactVertexVector(*m);
        tri::UpdateTopology<MyMesh>::VertexFace(*m);
        printf("\t FN: %d\n", m->FN());
        printf("\t VN: %d\n", m->VN());
        return count;
    }

    class EdgeInfo{
    public:
//        vector<double> angles;
        vector<MyMesh::VertexPointer> oppositeVert;
        vector<MyMesh::FacePointer> faces;
        vector<bool> mask;
        vector<int> faceV0Idx;
//        MyMesh::VertexPointer edge_verts[2];
        int count;
        EdgeInfo(){
            count = 1;
            faces.reserve(2);
            oppositeVert.reserve(2);
            mask.reserve(2);
        }
        void addFace(MyMesh::FacePointer f){faces.emplace_back(f);}
        void addOppositeVert(MyMesh::VertexPointer v){oppositeVert.emplace_back(v);}
        void addMask(bool ma){ mask.emplace_back(ma);}
        void addV0Idx(int i){faceV0Idx.emplace_back(i);}
        MyMesh::FacePointer getFace(int i){
            assert(i<faces.size() && i>=0);
            return faces[i];
        }
        MyMesh::VertexPointer getOppositeVert(int i){
            assert(i < oppositeVert.size() && i >= 0);
            return oppositeVert[i];
        }
        bool getMask(int i){
            assert(i < mask.size() && i>=0);
            return mask[i];
        }
        int getV0Idx(int i){
            assert(i < faceV0Idx.size() && i >= 0);
            return faceV0Idx[i];
        }
    };

    struct EdgeComp{
        bool operator()(array<MyMesh::VertexPointer,2> e1, array<MyMesh::VertexPointer,2> e2){
            if(e1[0] != e2[0]) return e1[0] < e2[0];
            else return e1[1] < e2[1];
        }
    };

    int edge_split(){
        printf("edge split: ");
        vector<triVert> splittingFaceVerts;
        vector<tuple<MyMesh::CoordType, MyMesh::CoordType>> newCoords;
        vector<bool> mask; //a mask for identify if the face uses the same splitting with the previous face
        bool flag = false;
        double angle_threshold = 110*(3.14/180); //110° in radiant
        map<array<MyMesh::VertexPointer,2>, EdgeInfo, EdgeComp> allEdges;
        for(auto fit = m->face.begin(); fit != m->face.end(); fit++){
            for(int i = 0; i < fit->VN(); i++){
                array<MyMesh::VertexPointer, 2> vertPair = {fit->V0(i), fit->V1(i)}; //aka edge
                if(vertPair[0] > vertPair[1]) swap(vertPair[0], vertPair[1]);
                if(allEdges.find(vertPair) != allEdges.end()) {
                    EdgeInfo *edgeInfo = &allEdges[vertPair];
                    if(edgeInfo->count > 1){
                        printf("arco strambo!\n");
                    }
                    edgeInfo->addFace(fit.base());
                    edgeInfo->addV0Idx(i);
                    edgeInfo->addOppositeVert(fit->V2(i));

                    if(edgeInfo->getOppositeVert(0) == edgeInfo->getOppositeVert(1)){
                        printf("faccia sovrapposta!\n");
                        edgeInfo->addMask(false);
                    }
                    else
                        edgeInfo->addMask(true);
                    edgeInfo->count++;
                }
                else {
                    EdgeInfo edgeInfo;
                    edgeInfo.addFace(fit.base());
                    edgeInfo.addV0Idx(i);
                    edgeInfo.addOppositeVert(fit->V2(i));
                    edgeInfo.addMask(true);
                    allEdges.insert(pair<array<MyMesh::VertexPointer, 2>, EdgeInfo>(vertPair, edgeInfo));
                }
            }
        }

        for(auto pair : allEdges){
            array<MyMesh::VertexPointer,2> vertPair = pair.first;
            EdgeInfo *info = &pair.second;
            if(info->count < 2) continue;
            double length = Distance(vertPair[0]->cP(), vertPair[1]->cP());
//            printf("length: %f\n", length);
            if(info->count > 1 && length >= zero_TH && !info->getFace(0)->IsD() && !info->getFace(1)->IsD()){
//            if(count > 1 && info.length >= zero_TH){
                double angle0, angle1;
                int idx[2], i = 0; //per determinare le 2 faccie non sovrapposte
                for(int j = 0; j < info->mask.size(); j++){
                    if(info->getMask(j)) idx[i++] = j;
                }
                angle0 = Angle(vertPair[0]->cP() - info->getOppositeVert(idx[0])->cP(), vertPair[1]->cP() - info->getOppositeVert(idx[0])->cP());
                angle1 = Angle(vertPair[0]->cP() - info->getOppositeVert(idx[1])->cP(), vertPair[1]->cP() - info->getOppositeVert(idx[1])->cP());
//                printf("angle0: %f\tangle1:%f\t threshold:%f\n", angle0, angle1, angle_threshold);
                if(angle0 >= angle_threshold && angle1 >= angle_threshold) { // ill-formed triangle
                    flag = !flag;
                    MyMesh::VertexPointer side_vert = angle0 > angle1 ? info->getOppositeVert(idx[0]) : info->getOppositeVert(idx[1]);
                    MyMesh::CoordType projector = (vertPair[0]->cP() - vertPair[1]->cP()).normalized();
                    MyMesh::CoordType projectee = side_vert->cP() - vertPair[0]->cP();
                    double t = projector.dot(projectee);    //a coefficient for the new vertex coordinates
                    MyMesh::CoordType mat0 = vert_mat[vertPair[0]];
                    MyMesh::CoordType mat1 = vert_mat[vertPair[1]];
                    MyMesh::CoordType mat_projector = (mat1 - mat0).normalized();
                    //first is the coordinate of the vertex, the second is the MAT coordinate associated with
                    newCoords.emplace_back(tuple<MyMesh::CoordType, MyMesh::CoordType>(vertPair[0]->cP() + projector * t, mat0 + mat_projector * t));
                    for(auto f = info->faces.begin(); f != info->faces.end(); f++){
                        face::VFDetach(**f);
                        tri::Allocator<MyMesh>::DeleteFace(*m, **f);
                    }
#if COLLAPESER_DEBUG
                    printf("-------\n");
                    PRINTP(info.edge_verts[0]->cP())
                    PRINTP(info.edge_verts[1]->cP())
                    PRINTP(info.oppositeVert[0]->cP())
                    PRINTP(info.oppositeVert[1]->cP())
                    assert(info.edge_verts[0] != info.edge_verts[1] && info.edge_verts[0] != info.oppositeVert[0] && info.edge_verts[0] != info.oppositeVert[1]);
                    assert(info.edge_verts[1] != info.oppositeVert[0] && info.edge_verts[1] != info.oppositeVert[1]);
                    assert(info.oppositeVert[0] != info.oppositeVert[1]);
#endif
                    for(int j = 0; j < info->faces.size(); j++){
                        int V0idx = info->getV0Idx(j);
                        splittingFaceVerts.emplace_back(triVert(info->getFace(j)->V0(V0idx),info->getFace(j)->V1(V0idx),info->getFace(j)->V2(V0idx)));
                        mask.emplace_back(flag);
                    }
                }
            }
        }

        if(newCoords.size() > 0) {
            tri::UpdateTopology<MyMesh>::VertexFace(*m);
            //allocating new vertices and new faces
            vector<triVert> newFaces;
            tri::Allocator<MyMesh>::PointerUpdater<MyMesh::VertexPointer> vpu;
            MyMesh::VertexIterator vi = tri::Allocator<MyMesh>::AddVertices(*m, newCoords.size(), vpu);
            assert(mask.size() == splittingFaceVerts.size());
            int sptfaceCount = 1;
            int coordCount = 1;
            auto current_flag = mask.begin();
            auto splittingFaceIt = splittingFaceVerts.begin();
            auto coordsIt = newCoords.begin();

            MyMesh::CoordType coord = get<0>(*newCoords.begin());
            MyMesh::CoordType new_mat = get<1>(*newCoords.begin());
            vi->P() = coord;

            int n = 0;
            while (vi != m->vert.end()) {
                MyMesh::VertexPointer v0 = get<0>(*splittingFaceIt);
                MyMesh::VertexPointer v1 = get<1>(*splittingFaceIt);
                MyMesh::VertexPointer v2 = get<2>(*splittingFaceIt);

                if (vpu.NeedUpdate()) {
                    vpu.Update(v0);
                    vpu.Update(v1);
                    vpu.Update(v2);
                }
//            assert(tri::Index(*m, v0) < m->vert.size() && tri::Index(*m, v0) >= 0);
//            assert(tri::Index(*m, v1) < m->vert.size() && tri::Index(*m, v1) >= 0);
//            assert(tri::Index(*m, v2) < m->vert.size() && tri::Index(*m, v2) >= 0);
//            assert(tri::Index(*m, vi.base()) < m->vert.size() && tri::Index(*m, vi.base()) >= 0);
//            printf("%p\t%p\t%p\n", v0, v1, vi.base());
                newFaces.emplace_back(triVert(vi.base(), v2, v0));
//            printf("%p\t%p\t%p\n", vi.base(), v1, v2);
                newFaces.emplace_back(triVert(vi.base(), v1, v2));
//
//            printf("current flag: ");
//            printf(*current_flag ? "1\t":"0\t");
//            printf("next flag: ");
//            printf(*next(current_flag) ? "1\n":"0\n");
                if (next(current_flag) == mask.end() || *next(current_flag) != *current_flag) {
                    vert_mat[vi] = new_mat;
                    isSplitted[vi] = true;
                    coord = get<0>(*coordsIt);
                    new_mat = get<1>(*coordsIt);
                    vi->P() = coord;
                    //setting properties of vi
                    omega_L[vi] = omega_L_0;
                    omega_H[vi] = omega_H_0;
                    omega_M[vi] = 0;

                    coordsIt++;
                    vi++;
                    n++;
                }
//            previous_flag = *current_flag;
                splittingFaceIt++;
                current_flag++;
                sptfaceCount++;
//            printf("split Face count: %d\n",sptfaceCount++);

            }

//            cout << "new faces size: " << newFaces.size() << endl;
//            cout << "coordinates size: " << newCoords.size() << endl;

            MyMesh::FaceIterator fi = tri::Allocator<MyMesh>::AddFaces(*m, newFaces.size());
            for (auto nfi = newFaces.begin(); nfi != newFaces.end() && fi != m->face.end(); nfi++, fi++) {
                fi->V(0) = get<0>(*nfi);
                fi->V(1) = get<1>(*nfi);
                fi->V(2) = get<2>(*nfi);
                fi->C() = Color4b::Red;
            }
        }
        printf("splitted %lu edges\n", newCoords.size());
        tri::UpdateTopology<MyMesh>::VertexFace(*m);
//        printf("vert %lu, %lu\n", m->vert.size(), m->VN());
//        printf("face %lu, %lu\n", m->face.size(), m->FN());
        tri::Allocator<MyMesh>::CompactFaceVector(*m);
        tri::Allocator<MyMesh>::CompactVertexVector(*m);
//        printf("--------end split--------\n");
        return newCoords.size();
    }

    void detect_degeneracies(){
        printf("detect degeneracies: ");
        double length = edge_threshold()/10;
        for (auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if (isFixed[vi] || vi->IsD()) continue;

            int counter = 0;
            vector<MyMesh::VertexPointer> adjacentVertices;
            vcg::face::VVStarVF<MyMesh::FaceType>(vi.base(), adjacentVertices);

            for(auto adjV = adjacentVertices.begin(); adjV != adjacentVertices.end(); adjV++){
                double edge_length = Distance(vi->P(), adjV.operator*()->P());
                VertexPair vpair(vi.base(), (*adjV.base()));
                if(edge_length <= length && collapser.LinkConditions(vpair)){
//                if(edge_length <= length){
                    counter++;
                }
            }

            if (counter >= 2) {
                isFixed[vi] = true;
                omega_L[vi] = 0;
                omega_H[vi] = 1.0/zero_TH;
                omega_M[vi] = 0;
            } else {
                isFixed[vi] = false;
            }
        }
        int n = 0;
        for(auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if(isFixed[vi.base()])
                n++;
        }
//        printf("-------------detect degeretaties-------------\n");
        printf("total fixed vert: %d\n", n);
    }
};


#endif //MCFSKET_COLLAPSER_H
