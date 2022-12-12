//
// Created by boy on 04/12/22.
//

#ifndef MCFSKET_COLLAPSER_H
#define MCFSKET_COLLAPSER_H

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

    void compute(string filename){
        createLHS();
        createRHS();
        solveLS();
        MyMesh::PerVertexAttributeHandle<MyMesh::CoordType> sol = tri::Allocator<MyMesh>::GetPerVertexAttribute<MyMesh::CoordType>(*m, string("solution"));

        for (auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            vi->P() = sol[vi];
        }
        cout << "meso skel created"<<endl;
        tri::io::ExporterOFF<MyMesh>::Save(*m, filename.c_str(), tri::io::Mask::IOM_FACECOLOR);

        update_omega();

        // todo
        // update contraints
        // update topology
        // detect degeneracies
    }

    void createLHS(){
        laplacian = tri::Allocator<MyMesh>::GetPerMeshAttribute<laplacian_triple>(*m, string("laplacian"));
        nrows = 3 * m->VN();
        ncols = m->VN();
        LHS.resize(nrows, ncols);
        RHS = Eigen::MatrixXd::Zero(nrows, 3);
        X = Eigen::MatrixXd::Zero(ncols, 3);

        vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(ncols*9);

        for (auto l = laplacian().begin(); l != laplacian().end(); l++){
//            printf("row: %d  col %d  ==>%f\n", vert_idx[get<0>(*l)], vert_idx[get<1>(*l)], get<2>(*l));
            triplets.emplace_back(Eigen::Triplet<double>(vert_idx[get<0>(*l)], vert_idx[get<1>(*l)], get<2>(*l)));
        }

        for(auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            triplets.emplace_back(Eigen::Triplet<double>(vert_idx[*vi]+ncols, vert_idx[*vi], omega_H[*vi]));
        }
        for(auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            triplets.emplace_back(Eigen::Triplet<double>(vert_idx[*vi]+2*ncols, vert_idx[*vi], omega_M[*vi]));
        }

        LHS.setFromTriplets(triplets.begin(), triplets.end());
    }

    void createRHS(){
        for(auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            MyMesh::CoordType c = (float)omega_H[vi] * vi->P();
            RHS.row(ncols + vert_idx[vi]) = Eigen::Vector3d(c.X(), c.Y(), c.Z());
        }

        for(auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
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
            return;
        }

        MyMesh::PerVertexAttributeHandle<MyMesh::CoordType> sol = tri::Allocator<MyMesh>::GetPerVertexAttribute<MyMesh::CoordType>(*m, string("solution"));
        for(auto vi=m->vert.begin(); vi != m->vert.end(); vi++){
            Eigen::Vector3d p = X.row(vert_idx[vi]);
            sol[vi] = MyMesh::CoordType(p[0], p[1], p[2]);
        }
    }

    void update_omega(){
        for (auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
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


    void update_topology(){
        int count=0;
        //edge collap
        for(auto eit = m->edge.begin(); eit != m->edge.end(); eit++){
            auto v0 = eit->V0(0);
            auto v1 = eit->V1(0);
            double edge_length = Distance(v0->P(), v1->P());
            if(edge_length < edge_threshold()){
                if(!eit->IsD()) { // and collapsable
                    v1->P() = (v1->P() + v0->P()) / 2;

                    auto mat0 = vert_mat[v0];
                    auto mat1 = vert_mat[v1];
                    double dist1 = Distance(mat0, v1->P());
                    double dist2 = Distance(mat1, v1->P());

                    if (dist1 < dist2)
                        vert_mat[v1] = vert_mat[v0];
                    else
                        vert_mat[v1] = vert_mat[v1];
                    //todo collapse the edge
                }
            }
        }

        //vertex split
        double angle_threshold = 110*(3.14/180); //110° in radiant
        for (auto eit = m->edge.begin(); eit != m->edge.end(); eit++){
            MyMesh::FaceType *edgeFP = eit->EFp();
            vcg::face::Pos<MyMesh::FaceType> p(edgeFP, eit->EFi());
            MyMesh::VertexType *tobeSplitted;

            //calculate the 2 opposite angles of the edge

            //angle of the first halfedge
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

            MyMesh::VertexType newV = *tri::Allocator<MyMesh>::AddVertex(*m, A->P() + projector * t);
            //todo split the edge on newV

            MyMesh::CoordType mat0 = vert_mat[A];
            MyMesh::CoordType mat1 = vert_mat[B];
            MyMesh::CoordType mat_projector = (mat1 - mat0).normalized();
            vert_mat[newV] = mat0 + mat_projector * t;
            isSplitted[newV] = true;
        }
    }

    void detect_degeneracies(){
        double length = edge_threshold()/10;
        for (auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if (isFixed[vi]) continue;

            bool fix = false;
            int counter = 0;
            MyFace* start = vi->VFp();
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
