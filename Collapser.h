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

    MyMesh::PerVertexAttributeHandle<int> vert_idx;
    MyMesh::PerVertexAttributeHandle<MyMesh::CoordType> vert_mat;
    MyMesh::PerVertexAttributeHandle<bool> isFixed;
    MyMesh::PerVertexAttributeHandle<bool> isSplitted;
    MyMesh::PerMeshAttributeHandle<laplacian_triple> laplacian;
    MyMesh::PerMeshAttributeHandle<double> edge_threshold;
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

        vert_idx = tri::Allocator<MyMesh>::GetPerVertexAttribute<int>(*m, string("map_vert_idx"));
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
//#if COLLAPESER_DEBUG
//            printf("row: %d  col %d  ==>%f\n", vert_idx[get<0>(*l)], vert_idx[get<1>(*l)], get<2>(*l));
//#endif
            triplets.emplace_back(Eigen::Triplet<double>(vert_idx[get<0>(*l)], vert_idx[get<1>(*l)], get<2>(*l)));
        }

        for(auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
//#if COLLAPESER_DEBUG
//            printf("row: %d  col %d  ==>%f\n", vert_idx[*vi]+ncols, vert_idx[*vi], omega_H[*vi]);
//#endif
            triplets.emplace_back(Eigen::Triplet<double>(vert_idx[*vi]+ncols, vert_idx[*vi], omega_H[*vi]));
        }
        for(auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
//#if COLLAPESER_DEBUG
//            printf("row: %d  col %d  ==>%f\n", vert_idx[*vi]+2*ncols, vert_idx[*vi], omega_H[*vi]);
//#endif
            triplets.emplace_back(Eigen::Triplet<double>(vert_idx[*vi]+2*ncols, vert_idx[*vi], omega_M[*vi]));
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
        //todo
        // there is some col index that goes out of the bound
        LHS.setFromTriplets(triplets.begin(), triplets.end());
    }

    void createRHS(){
        for(auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
            MyMesh::CoordType c = (float)omega_H[vi] * vi->P();
            RHS.row(ncols + vert_idx[vi]) = Eigen::Vector3d(c.X(), c.Y(), c.Z());
        }

        for(auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
            MyMesh::CoordType c = (float) omega_M[vi] * vert_mat[vi];
            RHS.row(2*ncols + vert_idx[vi]) = Eigen::Vector3d(c.X(), c.Y(), c.Z());
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
            Eigen::Vector3d p = X.row(vert_idx[vi]);
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
                        printf("V1 index: %d\n", vert_idx[V1]);
#endif
                        if (dist1 < dist2)
                            vert_mat[V1] = vert_mat[V0];
                        else
                            vert_mat[V1] = vert_mat[V1];
                        vcg::face::Pos<MyMesh::FaceType> p(fit.base(), V0);
                        vector<triVert> temp;
                        edge_collapse(p, temp);
                        newFaces.insert(newFaces.end(), temp.begin(), temp.end());
                        count++;
                    }
                }
            }
        }
        // create new added faces
//        cout << "face number before adding new faces:"<< m->FN()<<endl;
//        cout << "f container size: "<<m->face.size()<<endl;
        for(auto it = newFaces.begin(); it != newFaces.end(); it++){
            tri::Allocator<MyMesh>::AddFace(*m,get<0>(*it), get<1>(*it), get<2>(*it))->C() = Color4b::Red;
        }
//        cout << "face number after added new faces:"<< m->FN()<<endl;
//        cout << "f container size: "<<m->face.size()<<endl;
        printf("collapsed %d edges\n", count);

        //cleaning and updating
        tri::UpdateTopology<MyMesh>::VertexFace(*m);
        tri::UpdateTopology<MyMesh>::FaceFace(*m);
        tri::Allocator<MyMesh>::CompactFaceVector(*m);
        tri::Allocator<MyMesh>::CompactVertexVector(*m);
        cout << "topology updated, face and vertex cleaned"<<endl;

        //vertex split
        int counts = 0;
        tri::UpdateFlags<MyMesh>::FaceClearV(*m);
        double angle_threshold = 110*(3.14/180); //110° in radiant
        MyMesh::VertexType *tobeSplitted;
        MyMesh::VertexType *V;
        for(auto fit = m->face.begin(); fit != m->face.end(); fit++){
            if(!fit->IsV() && !fit->IsD()){
                for(int i = 0; i < fit->VN(); i++) {
                    V = fit->V0(i);
                    vcg::face::Pos<MyMesh::FaceType> p(V->VFp(), V);
                    MyMesh::VertexType *A = p.V();
                    MyMesh::VertexType *B = p.VFlip();

                    p.FlipE();
                    p.FlipV();
                    MyMesh::VertexType *C = p.V();
                    double alpha0 = p.AngleRad();
                    p.FlipV();
                    p.FlipE();//return on the starting edge

                    //go to the other face in order to calculate the other opposite angle (angle of the second halfedge)
                    p.FlipF();
                    p.FlipE();
                    p.FlipV();
                    MyMesh::VertexType  *D = p.V();
                    double alpha1 = p.AngleRad();
                    p.FlipV();
                    p.FlipE();
                    if(alpha0 < angle_threshold || alpha1 < angle_threshold) continue;
                    double a,b,c;
                    c = Distance(A->P(), B->P());
                    b = Distance(A->P(), C->P());
                    a = Distance(B->P(), C->P());

                    if(a < zero_TH || b < zero_TH || c < zero_TH){ //edge too short
                        continue;
                    }
                    if(alpha0 > alpha1)
                        tobeSplitted = C;
                    else
                        tobeSplitted = D;

                    MyMesh::CoordType projector = (B->P()-A->P()).normalized();
                    MyMesh::CoordType projectee = tobeSplitted->P() - A->P();
                    double t = projector.dot(projectee);

                    MyMesh::VertexType *newV = &*tri::Allocator<MyMesh>::AddVertex(*m, A->P() + projector * t);
                    vertex_split(p, newV);

                    MyMesh::CoordType mat0 = vert_mat[A];
                    MyMesh::CoordType mat1 = vert_mat[B];
                    MyMesh::CoordType mat_projector = (mat1 - mat0).normalized();
                    vert_mat[newV] = mat0 + mat_projector * t;
                    isSplitted[newV] = true;
                    counts++;
                }
            }
        }

        printf("splitted %d verteces\n", counts);
        tri::UpdateTopology<MyMesh>::VertexFace(*m);
        tri::UpdateTopology<MyMesh>::FaceFace(*m);
        tri::Allocator<MyMesh>::CompactFaceVector(*m);
        tri::Allocator<MyMesh>::CompactVertexVector(*m);
        cout << "topology updated, face and vertex cleaned"<<endl;
    }

    void edge_collapse(vcg::face::Pos<MyMesh::FaceType> &p, vector<triVert> &newFaceVerts){
        MyMesh::VertexType *anchorVert = p.VFlip();
        MyMesh::FaceType *start_face = p.F();
        vector<MyMesh::VertexType*> to_be_connected;
        do{
            if(!p.F()->IsD()){
//                p.F()->SetD();
                tri::Allocator<MyMesh>::DeleteFace(*m, *p.F());
            }
            p.FlipE();
            if (p.VFlip() != anchorVert) to_be_connected.emplace_back(p.VFlip());
            p.FlipF();
        }while(p.F() != start_face);
#if COLLAPESER_DEBUG
        cout <<"to be connected size:"<< to_be_connected.size()<<endl;
#endif
        newFaceVerts.reserve(to_be_connected.size()-1);
        for(int i = 0; i < to_be_connected.size()-1; i++){
            newFaceVerts.emplace_back(triVert(anchorVert, to_be_connected[i], to_be_connected[i+1]));
        }
//        p.V()->SetD();
        tri::Allocator<MyMesh>::DeleteVertex(*m, *p.V());
    }

    void vertex_split(vcg::face::Pos<MyMesh::FaceType> &p, MyMesh::VertexType *newV){
        MyMesh::FaceType *faceA = p.F();
        MyMesh::FaceType *faceB = p.FFlip();
        MyMesh::VertexType *A, *B, *C, *D;
        p.FlipE();
        A = p.VFlip();
        p.FlipE();

        p.FlipF();
        p.FlipE();
        C = p.VFlip();
        p.FlipE();
        B = p.V();
        D = p.VFlip();

        tri::Allocator<MyMesh>::DeleteFace(*m, *faceA);
        tri::Allocator<MyMesh>::DeleteFace(*m, *faceB);
//        faceA->SetD();
//        faceB->SetD();

        tri::Allocator<MyMesh>::AddFace(*m, newV, A, B);
        tri::Allocator<MyMesh>::AddFace(*m, newV, A, D);
        tri::Allocator<MyMesh>::AddFace(*m, newV, C, B);
        tri::Allocator<MyMesh>::AddFace(*m, newV, C, D);
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
