//
// Created by boy on 27/11/22.
//

#ifndef MCF_SKET_LAPLACEHELPER_HPP
#define MCF_SKET_LAPLACEHELPER_HPP

#include <unordered_map>
#include <math.h>
#include "utils.hpp"

class LaplaceHelper{
private:
    MyMesh *m;
    MyMesh::VertContainer vert;
    MyMesh::EdgeContainer edges;

    LaplaceHelper(MyMesh *m){
        this->m = m;
        vert = this->m->vert;
        edges = this->m->edge;
    }

public:
    double cotan(const Point3d a, const Point3d b){
        return (a.dot(b)) / cross(a, b).Norm();
    }
    Point3d cross(const Point3d a, const Point3d b){
        Point3d p;
        p.X() = a.Y() * b.Z() - a.Z() * b.Y();
        p.Y() = - (a.X() * b.Z() - a.Z() * b.X());
        p.Z() = a.X() * b.Y() - a.Y() * b.X();
        return p;
    }


    void get_vertex_laplace(){
        int n_vert = vert.size();
        vector<vector<double>> edge_weights;
        edge_weights.resize(n_vert);
        int vert_idx;
        MyMesh::VertexIterator vi;
        for(vi = vert.begin(), vert_idx = 0; vi != vert.end(); vi++, vert_idx++)
        {
            MyFace* start = vi.base()->VFp();
            vcg::face::Pos<MyFace> p(start,0,vi.base());
            double angle_sum = 0;
            Point3d p_sum = Point3d (0,0,0);
            double sin_alpha, cos_alpha, alpha, cot_alpha;
            double sin_beta, cos_beta, beta,cot_beta;
            double sum=0, sum_area=0;
            do
            {
                //todo
                //sumup area
                //at the end divide by 3

                //IN THE PREV EDGE
                //go to the opposite vert on the same edge
                p.FlipV();
                alpha = p.AngleRad();
                sincos(alpha, &sin_alpha, &cos_alpha);
                cot_alpha = sin_alpha/cos_alpha;
                p.FlipV();//return to the original vert

                //go to the NEXT EDGE of the adjacent face
                p.FlipE();
                p.FlipF();
                p.FlipE();

                //go to the opposite vert on the same edge
                p.FlipV();
                beta = p.AngleRad();
                sincos(beta, &sin_beta, &cos_beta);
                p.FlipV();
                cot_beta = sin_beta/cos_beta;
                sum += cot_alpha + cot_beta;

                //NEXT EDGE becomes CURRENT EDGE
                p.FlipE();
            }while(p.f!=start);
            //(2 * sum) / area of the vertex
        }

    }
};
#endif //MCF_SKET_LAPLACEHELPER_HPP
