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
    double compute_area(Point3d p, Point3d q, Point3d r, double angle){
        double pq = Distance(p, q);
        double pr = Distance(p, r);
        return 0.5 * pq * pr * sin(angle);
    }

    vector<vector<tuple<int, int, double>>> get_vertex_laplace(){
        int n_vert = vert.size();
        vector<vector<tuple<int, int, double>>> edge_weights;
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
            double w = 0, sum_area=0, sum = 0;

            //every cycle processes an edge
            int edge_idx = 0;
            do
            {
                double area = 0, angle;
                Point3d P, Q ,R; //3 vertices to be save for area computation
                P = p.V()->P();
                //IN THE PREV EDGE
                //go to the opposite vert on the same edge
                p.FlipV();
                alpha = p.AngleRad();
                sincos(alpha, &sin_alpha, &cos_alpha);
                cot_alpha = sin_alpha/cos_alpha;
                Q = p.V()->P();
                p.FlipV();//return to the original vert

                //go to the NEXT EDGE and the other vertex of the edge
                p.FlipE();
                p.FlipV();
                R = p.V()->P();
                p.FlipV();
                angle = p.AngleRad();
                //go to the adjacent face
                p.FlipF();
                p.FlipE();

                //go to the opposite vert on the same edge
                p.FlipV();
                beta = p.AngleRad();
                sincos(beta, &sin_beta, &cos_beta);
                p.FlipV();
                cot_beta = sin_beta/cos_beta;
                w = cot_alpha + cot_beta;

                area = compute_area(P, Q, R, angle);
                sum += w;
                sum_area += area;
                //inserting <i, j, w>
                edge_weights[vert_idx].emplace_back(tuple<int, int, double>(vert_idx, edge_idx, w));

                //NEXT EDGE becomes CURRENT EDGE
                p.FlipE();
            }while(p.f!=start);
            //inserting on the diagonal
            edge_weights[vert_idx].emplace_back(tuple<int, int, double>(vert_idx, vert_idx, -sum * sum_area/3));
        }
        return edge_weights;

    }
};
#endif //MCF_SKET_LAPLACEHELPER_HPP
