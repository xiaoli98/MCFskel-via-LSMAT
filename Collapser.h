//
// Created by boy on 04/12/22.
//

#ifndef MCFSKET_COLLAPSER_H
#define MCFSKET_COLLAPSER_H

#define COLLAPESER_DEBUG 1

#include <map>
#include "utils.hpp"

class Collapser{
    typedef MyMesh::PerVertexAttributeHandle<double> omega;
    typedef vector<tuple<MyMesh::VertexType*, MyMesh::VertexType*, double>> laplacian_triple;
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

    //some utility functions
    bool isInCommon(const MyMesh::VertexType *v, const MyMesh::FaceType *fb){
        for(int i = 0; i < 3; i++){
            if(fb->V(i) == v) return true;
        }
        return false;
    }
    int edge_index(const MyMesh::FaceType* f, const MyMesh::VertexType *v){
        for(int i = 0; i < 3; i++){
            if(f->V(i) == v) return i;
        }
        return -1;
    }

    int third_vert(const MyMesh::FaceType *f, const MyMesh::VertexType *v1, const MyMesh::VertexType *v2){
        for(int i = 0; i < 3; i++){
            if(f->V(i) != v1 && f->V(i) != v2) return i;
        }
    }
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
        detect_degeneracies();
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
    void update_topology(){
        //todo
        // Pos is not suitable for mesoskeletons
        // this is not manifold, so the navigation is not good
        // a possible idea is to use only VF an FF adjacency

        edge_collapse();
        edge_split();
        cout << "topology updated, face and vertex cleaned"<<endl;
    }

    void edge_collapse(){
        int count=0;
        //edge collap
        MyMesh::VertexType *V0, *V1;
        vector<triVert> newFaces;
        double edge_length;
        tri::UpdateFlags<MyMesh>::FaceClearV(*m);
        for(auto fit = m->face.begin(); fit != m->face.end(); fit++){
            if(!fit->IsV() && !fit->IsD()){
                fit->SetV();
                for(int i = 0; i < fit->VN(); i++){
                    V0 = fit->V0(i);
                    V1 = fit->V1(i);
                    edge_length = Distance(V0->P(), V1->P());
                    if(edge_length < edge_threshold()){
#if COLLAPESER_DEBUG
                        cout << "edge_length:" << edge_length << " threshold: " << edge_threshold()<<endl;
                        cout << "V0: "<<V0<<endl;
#endif
                        V1->P() = (V1->P() + V0->P()) / 2;

                        auto mat0 = vert_mat[V0];
                        auto mat1 = vert_mat[V1];
                        //find the closest medial axis to V1
                        double dist1 = Distance(mat0, V1->P());
                        double dist2 = Distance(mat1, V1->P());
#if COLLAPESER_DEBUG
                        printf("V1 index: %lu\n", tri::Index(*m, V1));
#endif
                        if (dist1 < dist2)
                            vert_mat[V1] = vert_mat[V0];
                        else
                            vert_mat[V1] = vert_mat[V1];


                        //todo
                        // replace pos with FFadj
                        vcg::face::Pos<MyMesh::FaceType> p(fit.base(), V0);
                        vector<triVert> temp;

                        MyMesh::VertexType *anchorVert = p.VFlip();
                        MyMesh::FaceType *start_face = p.F();
                        vector<MyMesh::VertexType*> to_be_connected;
                        do{
                            if(!p.F()->IsD()){
                                tri::Allocator<MyMesh>::DeleteFace(*m, *p.F());
                            }
                            p.FlipE();
                            if (p.VFlip() != anchorVert) to_be_connected.emplace_back(p.VFlip());
                            p.FlipF();
                        }while(p.F() != start_face);
#if COLLAPESER_DEBUG
                        cout <<"to be connected size:"<< to_be_connected.size()<<endl;
#endif
                        temp.reserve(to_be_connected.size()-1);
                        for(int j = 0; j < to_be_connected.size()-1; j++){
                            temp.emplace_back(triVert(anchorVert, to_be_connected[j], to_be_connected[j+1]));
                        }
//        p.V()->SetD();
                        tri::Allocator<MyMesh>::DeleteVertex(*m, *p.V());


                        newFaces.insert(newFaces.end(), temp.begin(), temp.end());
                        count++;
                    }
                }
            }
        }
        // create new added faces
        for(auto it = newFaces.begin(); it != newFaces.end(); it++){
            tri::Allocator<MyMesh>::AddFace(*m,get<0>(*it), get<1>(*it), get<2>(*it))->C() = Color4b::Red;
        }

        tri::UpdateTopology<MyMesh>::VertexFace(*m);
        tri::UpdateTopology<MyMesh>::FaceFace(*m);
        tri::Allocator<MyMesh>::CompactFaceVector(*m);
        tri::Allocator<MyMesh>::CompactVertexVector(*m);


        printf("collapsed %d edges\n", count);
    }

    void edge_split(){
        int counts = 0;
        tri::UpdateFlags<MyMesh>::FaceClearV(*m);
        tri::UpdateFlags<MyMesh>::VertexClearV(*m);
        double angle_threshold = 110*(3.14/180); //110° in radiant
        MyMesh::VertexPointer tobeSplittedVert;
        MyMesh::FacePointer tobeSplittedFace;
        tri::Allocator<MyMesh>::PointerUpdater<MyMesh::VertexPointer> vpu;

        for (auto fit = m->face.begin(); fit != m->face.end(); fit++){
            fit->SetV();
            for(int i = 0; i < 3; i++){
                MyMesh::FacePointer faceA = fit.base();
                MyMesh::FaceType *faceB = fit->FFp(i);
                if(faceB->IsV()) continue; //the adjacent face is visited, so this edge is already processed

                MyMesh::VertexPointer from_vert = faceA->V(i);
                MyMesh::VertexPointer to_vert = faceA->V((i+1)%3);
                MyMesh::VertexPointer left_vert, right_vert;
                left_vert = faceA->V(third_vert(faceA, from_vert, to_vert));
                right_vert = faceB->V(third_vert(faceB, from_vert, to_vert));
                //edge too short to be split, it should be collapsed
                if (Distance(from_vert->P(), to_vert->P()) < zero_TH) continue;

                assert(isInCommon(to_vert, faceB));
                double alpha0 = Angle(faceA->V(faceA->Prev(i))->cP()-from_vert->cP(), faceA->V(faceA->Next(i))->cP()-from_vert->cP());
                int j = edge_index(faceB, from_vert);
                double alpha1 = Angle(faceB->V(faceB->Prev(j))->cP()-from_vert->cP(), faceB->V(faceB->Next(j))->cP()-from_vert->cP());
                if(alpha0 < angle_threshold || alpha1 < angle_threshold) continue;
                printf("angle0: %f\t", alpha0);
                printf("angle1: %f\n", alpha1);
                if(alpha0 > alpha1) {
                    tobeSplittedVert = left_vert;
                    tobeSplittedFace = faceA;
                }
                else {
                    tobeSplittedVert = right_vert;
                    tobeSplittedFace = faceB;
                }

                MyMesh::CoordType projector = (to_vert->P()-from_vert->P()).normalized();
                MyMesh::CoordType projectee = tobeSplittedVert->P() - from_vert->P();
                double t = projector.dot(projectee);
                MyMesh::VertexType *newV = &*tri::Allocator<MyMesh>::AddVertices(*m, 1, vpu);
                MyMesh::CoordType newCoord = from_vert->P() + projector * t;
                newV->P() = newCoord;

                tri::Allocator<MyMesh>::DeleteFace(*m, *tobeSplittedFace);
                if(vpu.NeedUpdate()) {
                    vpu.Update(tobeSplittedVert);
                    vpu.Update(from_vert);
                    vpu.Update(to_vert);
                }

                printf("from_vert idx: %lu \n", tri::Index(*m, from_vert));
                printf("newV idx: %lu \n", tri::Index(*m, newV));
                printf("tobeSplittedVert idx: %lu \n", tri::Index(*m, tobeSplittedVert));
                printf("vert.front idx: %lu \n", tri::Index(*m, m->vert.front()));
                printf("vert.back idx: %lu \n", tri::Index(*m, m->vert.back()));
                tri::Allocator<MyMesh>::AddFace(*m, from_vert, newV, tobeSplittedVert)->SetV();
                tri::Allocator<MyMesh>::AddFace(*m,  tobeSplittedVert, to_vert, newV)->SetV();

                MyMesh::CoordType mat0 = vert_mat[from_vert];
                MyMesh::CoordType mat1 = vert_mat[to_vert];
                MyMesh::CoordType mat_projector = (mat1 - mat0).normalized();
                vert_mat[newV] = mat0 + mat_projector * t;
                isSplitted[newV] = true;
                counts++;
            }
        }

//        for(auto fit = m->face.begin(); fit != m->face.end(); fit++){
//            if(!fit->IsV() && !fit->IsD()){
//                fit->SetV();
//                for(int i = 0; i < 3; i++){
//
//                }
//                printf("------------------------------------------------------------------------------------------\n");
//                for(int i = 0; i < fit->VN(); i++) {
//                    V = fit->V0(i);
//
//                    if(V->IsV()) continue;
//                    V->SetV();
//                    vcg::face::Pos<MyMesh::FaceType> p(fit.base(), V);
//                    MyMesh::FaceType *FaceA = V->VFp();
//                    MyMesh::FaceType *FaceB;
//                    MyMesh::VertexType *vA, *vB, *vC, *vD;
//                    vA = p.V();
//                    vB = p.VFlip();
//
//                    //the first halfedge and the relative opposite angle
//                    p.FlipE();
//                    p.FlipV();
//                    vC = p.V();
//                    MyMesh::VertexType *tempA, *tempB, *tempC;
//                    tempA = fit->V0(i);
//                    tempB = fit->V1(i);
//                    tempC = fit->V2(i);
//#if COLLAPESER_DEBUG
//                    double AB = Distance(vA->P(), vB->P());
//                    double AC = Distance(vA->P(), vC->P());
//                    double BC = Distance(vB->P(), vC->P());
////                    printf("AB: %f \t AC: %f \t BC: %f \n", AB, AC, BC);
//                    double X = (BC*BC + AC*AC - AB*AB) / (2 * AC * BC);
////                    printf("X: %f\n", X);
//                    double angleC = acos(X);
////                    printf("angle in C: %f\n", angleC);
//
//#endif
//                    double alpha0 = p.AngleRad();
////                    double alpha0 = angleC;
//                    p.FlipV();
//                    p.FlipE();//return on the starting edge
//
//                    //go to the other face in order to calculate the other opposite angle (angle of the second halfedge)
//                    p.FlipF();
//                    FaceB = p.F();
//                    p.FlipE();
//                    p.FlipV();
//                    vD = p.V();
//                    double alpha1 = p.AngleRad();
//                    p.FlipV();
//                    p.FlipE();
//#if COLLAPESER_DEBUG
////                    printf("angle threshold: %f\n", angle_threshold);
////                    printf("alpha0: %f \t alpha1: %f \n", alpha0, alpha1);
//#endif
//                    if(alpha0 < angle_threshold || alpha1 < angle_threshold) continue;
//                    double a,b,c;
//                    c = Distance(vA->P(), vB->P());
//                    b = Distance(vA->P(), vC->P());
//                    a = Distance(vB->P(), vC->P());
//
//                    if(a < zero_TH || b < zero_TH || c < zero_TH){ //edge too short, need collapse instead of split
//                        continue;
//                    }
//                    if(alpha0 > alpha1) {
//                        tobeSplittedVert = vC;
//                        tobeSplittedFace = FaceA;
//                    }
//                    else {
//                        tobeSplittedVert = vD;
//                        tobeSplittedFace = FaceB;
//                    }
//
//                    MyMesh::CoordType projector = (vB->P()-vA->P()).normalized();
//                    MyMesh::CoordType projectee = tobeSplittedVert->P() - vA->P();
//                    double t = projector.dot(projectee);
//
//                    MyMesh::VertexType *newV = &*tri::Allocator<MyMesh>::AddVertex(*m, vA->P() + projector * t);
//                    printf("-----------------------------------------------------------------\n");
//                    int abc = -1;
//                    for(auto vit = m->vert.begin(); vit != m->vert.end()-1; vit++){
//                        if(vert_idx[vit] > abc) {
//                            abc = vert_idx[vit];
//                        }
//                        else {
//                            printf("false\n");
//                            printf("%d %d\n", vert_idx[vit], abc);
//                        }
//                    }
//                    printf("-----------------------------------------------------------------\n");
//
//                    printf("vA_idx: %d\n", vert_idx[vA]);
//                    printf("vB_idx: %d\n", vert_idx[vB]);
//                    printf("tobeSplittedVert_idx: %d\n", vert_idx[tobeSplittedVert]);
//
//                    tri::Allocator<MyMesh>::DeleteFace(*m, *tobeSplittedFace);
//                    tri::Allocator<MyMesh>::AddFace(*m, vA, newV, tobeSplittedVert)->SetV();
//                    tri::Allocator<MyMesh>::AddFace(*m,  tobeSplittedVert, vB, newV)->SetV();
//
//                    MyMesh::CoordType mat0 = vert_mat[vA];
//                    MyMesh::CoordType mat1 = vert_mat[vB];
//                    MyMesh::CoordType mat_projector = (mat1 - mat0).normalized();
//                    vert_mat[newV] = mat0 + mat_projector * t;
//                    isSplitted[newV] = true;
//                    counts++;
//                    break;
//                }
//            }
//        }

        printf("splitted %d edges\n", counts);
        tri::UpdateTopology<MyMesh>::VertexFace(*m);
        tri::UpdateTopology<MyMesh>::FaceFace(*m);
        tri::Allocator<MyMesh>::CompactFaceVector(*m);
        tri::Allocator<MyMesh>::CompactVertexVector(*m);
//        MyMesh::FaceType *faceA = p.F();
//        MyMesh::FaceType *faceB = p.FFlip();
//        MyMesh::VertexType *A, *B, *C, *D;
//        p.FlipE();
//        A = p.VFlip();
//        p.FlipE();
//
//        p.FlipF();
//        p.FlipE();
//        C = p.VFlip();
//        p.FlipE();
//        B = p.V();
//        D = p.VFlip();
//
//        tri::Allocator<MyMesh>::DeleteFace(*m, *faceA);
//        tri::Allocator<MyMesh>::DeleteFace(*m, *faceB);
//        faceA->SetD();
//        faceB->SetD();
//        //todo
//        // try without SetV
//        tri::Allocator<MyMesh>::AddFace(*m, newV, A, B)->SetV();
//        tri::Allocator<MyMesh>::AddFace(*m, newV, A, D)->SetV();
//        tri::Allocator<MyMesh>::AddFace(*m, newV, C, B)->SetV();
//        tri::Allocator<MyMesh>::AddFace(*m, newV, C, D)->SetV();
    }

    void detect_degeneracies(){
        double length = edge_threshold()/10;
        for (auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if (isFixed[vi] || vi->IsD()) continue;

            int counter = 0;
            vector<MyMesh::VertexType*> adjacentsVertices;
            vcg::face::VVStarVF<MyMesh::FaceType>(vi.base(), adjacentsVertices);
            for(auto adjV = adjacentsVertices.begin(); adjV != adjacentsVertices.end(); adjV++){
                double edge_length = Distance(vi->P(), adjV.operator*()->P());
                if(edge_length <= length){ //here should be another constraints which check if is collapsable, but not done
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
