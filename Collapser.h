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

//    unordered_map<MyMesh::VertexType*, Point3d> map_vert_mat;
    MyMesh::PerVertexAttributeHandle<int> vert_idx;
    MyMesh::PerVertexAttributeHandle<MyMesh::CoordType> vert_mat;
    MyMesh::PerMeshAttributeHandle<laplacian_triple> laplacian;
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
    }

    void compute(){
        createLHS();
        createRHS();
        solveLS();
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


};


#endif //MCFSKET_COLLAPSER_H
