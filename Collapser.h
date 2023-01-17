//
// Created by boy on 04/12/22.
//

#ifndef MCFSKET_COLLAPSER_H
#define MCFSKET_COLLAPSER_H

#define COLLAPESER_DEBUG 1

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
        update_omega();
        update_topology();
        tri::io::ExporterOFF<MyMesh>::Save(*m, "prova.off", tri::io::Mask::IOM_FACECOLOR);
        detect_degeneracies();

//        int counter = 0;
//        for(auto fit = m->face.begin(); fit != m->face.end(); fit++){
//            printf("-------------FACE %lu----------------\n", tri::Index(*m, fit.base()));
//            if(fit->IsD()) continue;
//            for(int i=0; i<fit->VN(); i++){
//                PRINTP(fit->V(i)->P());
//                f == f->FFp(e)->FFp(f->FFi(e))
//                MyMesh::FacePointer temp = fit->FFp(i)->FFp(fit->FFi(i));
//
//                if(fit.base() != temp) counter++;
//            }
//        }
//        cout << "non manifold face: "<<counter<<endl;
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

        for (auto l = laplacian().begin(); l != laplacian().end(); l++){
            triplets.emplace_back(Eigen::Triplet<double>(tri::Index(*m, get<0>(*l)), tri::Index(*m, get<1>(*l)), get<2>(*l)));
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

        if (!std::isfinite(X.norm())){
            throw EIGEN_EXCEPTIONS("Problem with linear least square solution.");
        }

        MyMesh::PerVertexAttributeHandle<MyMesh::CoordType> sol = tri::Allocator<MyMesh>::GetPerVertexAttribute<MyMesh::CoordType>(*m, string("solution"));
        for(auto vi=m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
            Eigen::Vector3d p = X.row(tri::Index(*m, *vi));
            sol[vi] = MyMesh::CoordType(p[0], p[1], p[2]);
        }
    }

    void update_omega(){
        for (auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
            if(isFixed[vi]){
                omega_L[vi] = 0;
                omega_H[vi] = 1.0/zero_TH;
                omega_M[vi] = 0;
                continue;
            }

            omega_L[vi] = omega_L_0;
            omega_H[vi] = omega_H_0;
            omega_M[vi] = omega_M_0;

            if(isSplitted[vi]){
                omega_L[vi] = omega_L_0;
                omega_H[vi] = omega_H_0;
                omega_M[vi] = 0;
            }
        }
    }

    typedef tuple<MyMesh::VertexType*,MyMesh::VertexType*,MyMesh::VertexType*> triVert;
    typedef tuple<MyMesh::VertexType*,MyMesh::VertexType*,MyMesh::VertexType*, MyMesh::VertexType*> quadVert;
    void update_topology(){
        edge_collapse();
        edge_split1();
        cout << "topology updated, face and vertex cleaned"<<endl;
    }

    void edge_collapse(){
        tri::EdgeCollapser<MyMesh, VertexPair> collapser;
        int count=0;
        MyMesh::VertexType *V0, *V1;
        double edge_length;
        tri::UpdateFlags<MyMesh>::VertexClearV(*m);

        vector<MyMesh::VertexPointer> adjVertices;
        for(auto vit = m->vert.begin(); vit != m->vert.end(); vit++) {
            if(!vit->IsD()) {
                face::VVStarVF<MyMesh::FaceType>(vit.base(), adjVertices);
                for (auto adjv = adjVertices.begin(); adjv != adjVertices.end(); adjv++) {
                    if ((*adjv)->IsD()) continue;
                    MyMesh::FacePointer fp = vit->VFp();
                    V0 = vit.base();
                    V1 = (*adjv);

                    edge_length = Distance(V0->P(), V1->P());
                    if (edge_length < edge_threshold()) {
#if COLLAPESER_DEBUG
                        cout << "edge_length:" << edge_length << " threshold: " << edge_threshold() << endl;
#endif
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

                        VertexPair vpair(V0, V1);
                        if(collapser.LinkConditions(vpair)){
                            collapser.Do(*m, vpair, V1->P());
                            count++;
                            break;
                        }
                    }
                }
            }
        }

        tri::UpdateTopology<MyMesh>::VertexFace(*m);
//        tri::UpdateTopology<MyMesh>::FaceFace(*m);
        printf("collapsed %d times\n", count);
        tri::Allocator<MyMesh>::CompactFaceVector(*m);
        tri::Allocator<MyMesh>::CompactVertexVector(*m);
    }

    class EdgeInfo{
    public:
        double angles[2];
        double length;
        MyMesh::VertexPointer oppositeVert[2];
        MyMesh::FacePointer faces[2];
        MyMesh::VertexPointer edge_verts[2];
        EdgeInfo(MyMesh::VertexPointer v1, MyMesh::VertexPointer v2){
            edge_verts[0] = v1;
            edge_verts[1] = v2;
            sort();
        }
        void sort(){
            if(edge_verts[0] > edge_verts[1]) swap(edge_verts[0], edge_verts[1]);
        }
    };

    struct EdgeComp{
        bool operator()(const EdgeInfo e1, const EdgeInfo e2){
            if(e1.edge_verts[0] < e2.edge_verts[0])
                return true;
            else if (e1.edge_verts[0] == e2.edge_verts[0]){
                if(e1.edge_verts[1] < e2.edge_verts[1])
                    return true;
                else
                    return false;
            }
            else
                return false;
        }
    };

    void edge_split1(){
        vector<quadVert> splittingFaceVerts;
        vector<tuple<MyMesh::CoordType, double>> newCoords;
        double angle_threshold = 110*(3.14/180); //110° in radiant
        map<EdgeInfo, int, EdgeComp> allEdges;
        for(auto fit = m->face.begin(); fit != m->face.end(); fit++){
            for(int i = 0; i < fit->VN(); i++){
                EdgeInfo edgeInfo(fit->V0(i), fit->V1(i));
                edgeInfo.length = Distance(fit->V0(i)->cP(), fit->V1(i)->cP());
                auto exist = allEdges.find(edgeInfo);
                if(exist != allEdges.end()) {
                    edgeInfo.faces[1] = fit.base();
                    edgeInfo.oppositeVert[1] = fit->V2(i);
                    edgeInfo.angles[1] = Angle(fit->V0(i)->cP() - fit->V2(i)->cP(), fit->V1(i)->cP() - fit->V2(i)->cP());
                    allEdges[edgeInfo]++;
                }
                else {
                    edgeInfo.faces[0] = fit.base();
                    edgeInfo.oppositeVert[0] = fit->V2(i);
                    edgeInfo.angles[0] = Angle(fit->V0(i)->cP() - fit->V2(i)->cP(), fit->V1(i)->cP() - fit->V2(i)->cP());
                    allEdges.insert(std::make_pair(edgeInfo, 1));
                }
            }
        }
        printf("total edges: %lu\n", allEdges.size());
        for(auto edge : allEdges){
            EdgeInfo info = edge.first;
            int count =edge.second;
            if(count > 1 && info.length >= zero_TH && !info.faces[0]->IsD() && !info.faces[1]->IsD()){
                if(info.angles[0] >= angle_threshold && info.angles[1] >= angle_threshold) { // ill-formed triangle
                    MyMesh::VertexPointer side_vert = info.angles[0] > info.angles[1] ? info.oppositeVert[0] : info.oppositeVert[1];
                    MyMesh::CoordType projector = (info.edge_verts[0]->cP() - info.edge_verts[1]->cP()).Normalize();
                    MyMesh::CoordType projectee = side_vert->cP() - info.edge_verts[0]->cP();
                    double t = projector.dot(projectee);    //a coefficient for the new vertex coordinates
                    newCoords.emplace_back(tuple<MyMesh::CoordType, double>(info.edge_verts[0]->cP() + projector * t, t));
                    tri::Allocator<MyMesh>::DeleteFace(*m, *info.faces[0]);
                    tri::Allocator<MyMesh>::DeleteFace(*m, *info.faces[1]);
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
                    splittingFaceVerts.emplace_back(quadVert(info.edge_verts[0], info.edge_verts[1], info.oppositeVert[0], info.oppositeVert[1]));
                }
            }
        }

        //allocating new vertices and new faces
        vector<triVert> newFaces;
        tri::Allocator<MyMesh>::PointerUpdater<MyMesh::VertexPointer> vpu;
        MyMesh::VertexIterator vi = tri::Allocator<MyMesh>::AddVertices(*m, newCoords.size(), vpu);

        assert(splittingFaceVerts.size() == newCoords.size());

        auto splittingFaceIt = splittingFaceVerts.begin();
        auto coordsIt = newCoords.begin();
        while(vi != m->vert.end()){
            MyMesh::CoordType coord = get<0>(*coordsIt);
            double t = get<1>(*coordsIt);
            vi->P() = coord;
            MyMesh::VertexPointer v0 = get<0>(*splittingFaceIt);
            MyMesh::VertexPointer v1 = get<1>(*splittingFaceIt);
            MyMesh::VertexPointer v2 = get<2>(*splittingFaceIt);
            MyMesh::VertexPointer v3 = get<3>(*splittingFaceIt);
            if(vpu.NeedUpdate()){
                vpu.Update(v0);
                vpu.Update(v1);
                vpu.Update(v2);
                vpu.Update(v3);
            }
            newFaces.emplace_back(triVert(v0, vi.base(), v2));
            newFaces.emplace_back(triVert(vi.base(), v1, v2));
            newFaces.emplace_back(triVert(v0, v3, vi.base()));
            newFaces.emplace_back(triVert(vi.base(), v3, v1));

            MyMesh::CoordType mat0 = vert_mat[v0];
            MyMesh::CoordType mat1 = vert_mat[v1];
            MyMesh::CoordType mat_projector = (mat1 - mat0).normalized();
            vert_mat[vi] = mat0 + mat_projector * t;
            isSplitted[vi] = true;

            vi++;
            splittingFaceIt++;
            coordsIt++;
        }
        MyMesh::FaceIterator  fi = tri::Allocator<MyMesh>::AddFaces(*m, newFaces.size());
        for(auto nfi = newFaces.begin(); nfi != newFaces.end() && fi!=m->face.end(); nfi++, fi++){
            fi->V(0) = get<0>(*nfi);
            fi->V(1) = get<1>(*nfi);
            fi->V(2) = get<2>(*nfi);
            fi->C() = Color4b::Red;
        }

        printf("splitted %lu edges\n", newCoords.size());
        tri::UpdateTopology<MyMesh>::VertexFace(*m);
        tri::Allocator<MyMesh>::CompactFaceVector(*m);
        tri::Allocator<MyMesh>::CompactVertexVector(*m);
    }

    void edge_split(){
        tri::UpdateFlags<MyMesh>::VertexClearV(*m);
        double angle_threshold = 110*(3.14/180); //110° in radiant
        vector<triVert> newFaces;               //contains 3 vertices where faces should be created
        vector<quadVert> splittingFaceVerts;     //deleted faces vertices
        vector<tuple<MyMesh::CoordType, double>> newCoords;
        MyMesh::FacePointer faceA, faceB;

        vector<MyMesh::FacePointer> fromVec;
        vector<int> idxFrom;

        double oppositeAngleA, oppositeAngleB;
        int vn;

        for(auto vit = m->vert.begin(); vit != m->vert.end(); vit++){
            vit->SetV();
            face::VFIterator<MyMesh::FaceType> vfi(vit.base());
            for(;!vfi.End(); ++vfi) {
                faceA = vfi.F();
//                if(faceA->IsD()) continue;
                vn = vfi.F()->VN();
                MyMesh::VertexPointer from_vert = faceA->V(vfi.I());
//                assert(vit.base() == faceA->V(vfi.I()));
                MyMesh::VertexPointer to_vert = faceA->V((vfi.I() + 1) % vn);
                if (to_vert->IsV()) continue; //edge already processed

                MyMesh::VertexPointer edge_oppositeVertA = faceA->V((vfi.I() + 2) % vn);
                MyMesh::VertexPointer edge_oppositeVertB;

                face::VFStarVF(from_vert, fromVec, idxFrom);
                auto idx = idxFrom.begin();
                //find the adjacent face and the vertex opposite to the edge
                bool faceFound = false;
                for (auto f = fromVec.begin(); f != fromVec.end(); f++, idx++) {
                    vn = (*f.base())->VN();
                    if ((to_vert == (*(f.base()))->V((*idx + 1) % vn) ||
                         to_vert == (*(f.base()))->V((*idx + 2) % vn)) && *f.base() != faceA) { faceB = *(f.base()); faceFound = true;}
                }
                for(int i = 0; i < vn; i++){
                    printf("face found: %d\n", faceFound);
                    if(faceB->V(i) != from_vert && faceB->V(i) != to_vert){
                        edge_oppositeVertB = faceB->V((i));
                        break;
                    }
                }
                oppositeAngleB = Angle(from_vert->cP() - edge_oppositeVertB->cP(), to_vert->cP() - edge_oppositeVertB->cP());

                if(Distance(from_vert->P(), to_vert->P()) < zero_TH) continue;//too short edge, should be collapsed

                oppositeAngleA = Angle(from_vert->cP() - edge_oppositeVertA->cP(), to_vert->cP() - edge_oppositeVertA->cP());
#if COLLAPESER_DEBUG
//                printf("Angle A: %f \t Angle B: %f \t threshold: %f\n", oppositeAngleA, oppositeAngleB, angle_threshold);
#endif
                if(oppositeAngleA < angle_threshold || oppositeAngleB < angle_threshold) continue;

                MyMesh::VertexPointer side_vert = oppositeAngleA > oppositeAngleB ?  edge_oppositeVertA : edge_oppositeVertB;
                MyMesh::CoordType projector = (to_vert->P()-from_vert->P()).Normalize();
                MyMesh::CoordType projectee = side_vert->P() - from_vert->P();
                double t = projector.dot(projectee);    //a coefficient for the new vertex coordinates

                newCoords.emplace_back(tuple<MyMesh::CoordType, double>(from_vert->P() + projector * t, t));
                if(!faceA->IsD()) tri::Allocator<MyMesh>::DeleteFace(*m, *faceA);
                if(!faceB->IsD()) tri::Allocator<MyMesh>::DeleteFace(*m, *faceB);
                splittingFaceVerts.emplace_back(quadVert (from_vert, to_vert, edge_oppositeVertA, edge_oppositeVertB));
            }
        }

        //allocating new vertices and new faces
        tri::Allocator<MyMesh>::PointerUpdater<MyMesh::VertexPointer> vpu;
        MyMesh::VertexIterator vi = tri::Allocator<MyMesh>::AddVertices(*m, newCoords.size(), vpu);

        assert(splittingFaceVerts.size() == newCoords.size());

        auto splittingFaceIt = splittingFaceVerts.begin();
        auto coordsIt = newCoords.begin();
        while(vi != m->vert.end()){
            MyMesh::CoordType coord = get<0>(*coordsIt);
            double t = get<1>(*coordsIt);
            vi->P() = coord;

            MyMesh::VertexPointer v0 = get<0>(*splittingFaceIt);
            MyMesh::VertexPointer v1 = get<1>(*splittingFaceIt);
            MyMesh::VertexPointer v2 = get<2>(*splittingFaceIt);
            MyMesh::VertexPointer v3 = get<3>(*splittingFaceIt);
            if(vpu.NeedUpdate()){
                vpu.Update(v0);
                vpu.Update(v1);
                vpu.Update(v2);
                vpu.Update(v3);
            }
            newFaces.emplace_back(triVert(v0, vi.base(), v2));
            newFaces.emplace_back(triVert(vi.base(), v1, v2));
            newFaces.emplace_back(triVert(v0, v3, vi.base()));
            newFaces.emplace_back(triVert(vi.base(), v3, v1));

            MyMesh::CoordType mat0 = vert_mat[v0];
            MyMesh::CoordType mat1 = vert_mat[v1];
            MyMesh::CoordType mat_projector = (mat1 - mat0).normalized();
            vert_mat[vi] = mat0 + mat_projector * t;
            isSplitted[vi] = true;

            vi++;
            splittingFaceIt++;
            coordsIt++;
        }
        MyMesh::FaceIterator  fi = tri::Allocator<MyMesh>::AddFaces(*m, newFaces.size());
        for(auto nfi = newFaces.begin(); nfi != newFaces.end() && fi!=m->face.end(); nfi++, fi++){
            fi->V(0) = get<0>(*nfi);
            fi->V(1) = get<1>(*nfi);
            fi->V(2) = get<2>(*nfi);
            fi->C() = Color4b::Red;
        }

        printf("splitted %lu edges\n", newCoords.size());
        tri::UpdateTopology<MyMesh>::VertexFace(*m);
        tri::Allocator<MyMesh>::CompactFaceVector(*m);
        tri::Allocator<MyMesh>::CompactVertexVector(*m);
    }

    void detect_degeneracies(){
        double length = edge_threshold()/10;
        for (auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if (isFixed[vi] || vi->IsD()) continue;

            int counter = 0;
            vector<MyMesh::VertexPointer> adjacentVertices;
            vcg::face::VVStarVF<MyMesh::FaceType>(vi.base(), adjacentVertices);

            for(auto adjV = adjacentVertices.begin(); adjV != adjacentVertices.end(); adjV++){
                double edge_length = Distance(vi->P(), adjV.operator*()->P());
                VertexPair vpair(vi.base(), (*adjV.base()));
                if(edge_length <= length && collapser.LinkConditions(vpair)){ //here should be another constraints which check if is collapsable, but not done
                    counter++;
                }
            }

            if (counter >= 2) {
                isFixed[vi] = true;
            } else {
                isFixed[vi] = false;
            }
        }
    }
};


#endif //MCFSKET_COLLAPSER_H
