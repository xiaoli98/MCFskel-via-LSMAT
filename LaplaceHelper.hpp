//
// Created by boy on 27/11/22.
//

#ifndef MCF_SKET_LAPLACEHELPER_HPP
#define MCF_SKET_LAPLACEHELPER_HPP

#include <unordered_map>
#include <math.h>
#include <tuple>
#include "utils.hpp"


class LaplaceHelper{
private:
    MyMesh *m;
    MyMesh::VertContainer vert;
    MyMesh::EdgeContainer edges;

    vector<vector<tuple<int, int, double>>> laplacian;
    vector<tuple<int, int, double>> laplacian_weights;
    unordered_map<MyMesh::VertexType* ,int> map_vert_idx;
public:
    LaplaceHelper(MyMesh *m){
        this->m = m;
        vert = this->m->vert;
        edges = this->m->edge;

    }

    const vector<vector<tuple<int, int, double>>> &getLaplacian() const {
        return laplacian;
    }

    const vector<tuple<int, int, double>> &getLaplacianWeights() const {
        return laplacian_weights;
    }

    double compute_area(Point3d p, Point3d q, Point3d r, double angle){
        double pq = Distance(p, q);
        double pr = Distance(p, r);
        return 0.5 * pq * pr * sin(angle);
    }

    void compute_laplace() {
        int n_vert = vert.size();
        laplacian.resize(n_vert);
        laplacian_weights.resize(n_vert);
        int vert_idx;
        tri::UpdateTopology<MyMesh>::VertexFace(*m);
        tri::UpdateTopology<MyMesh>::FaceFace(*m);
        tri::UpdateFlags<MyMesh>::VertexClearV(*m);
//        MyMesh::VertexIterator vi;
        MyMesh::FaceIterator fi;
        MyMesh::VertexPointer vp;
        int visited_vert = 0;

        ////creating mapping of vertex and the index
        //todo
        //  non ben capito ma i puntatori nei container di vertex sono diversi dai puntatori dei vertici nelle facce
        //  quindi inizializzo la mappa con usando i puntatori dei vertici delle facce
        int idx = 0;
        for(fi = m->face.begin(); fi != m->face.end(); fi++){
            for (int i = 0; i < fi->VN(); i++) {
                if (! fi->V(i)->IsD() && !fi->V(i)->IsV()) {
                    fi->V(i)->SetV();
                    map_vert_idx.insert(make_pair(fi->V(i), idx++));
                }
            }
        }

        tri::UpdateFlags<MyMesh>::VertexClearV(*m);
        for (fi = m->face.begin(); fi != m->face.end(); fi++) {
            for (int i = 0; i < fi->VN(); i++) {
                vp = fi->V(i);
                if (!vp->IsD() && !vp->IsV()) {
//                    map_vert_idx.insert(make_pair(vp,idx++));
//                    cout << "vp:"<<vp<<" ";
                    PRINTP(vp->P())
                    visited_vert++;
                    vp->SetV();
                    vcg::face::Pos<MyFace> p(fi.base(), i, vp);
                    double angle_sum = 0;
                    Point3d p_sum = Point3d(0, 0, 0);
                    double sin_alpha, cos_alpha, alpha, cot_alpha;
                    double sin_beta, cos_beta, beta, cot_beta;
                    double w = 0, sum_area = 0, sum = 0;
                    //every cycle processes an edge
                    int edge_idx = 0;
                    do {
                        double area = 0, angle;
                        Point3d P, Q, R; //3 vertices to be save for area computation
                        P = p.V()->P();
                        vert_idx = map_vert_idx[p.V()];
//                        cout << "vert_idx: "<<vert_idx<<endl;
//                        cout << "BEFORE FlipV: "<< p.V()<<endl;
                        //IN THE PREV EDGE
                        //go to the opposite vert on the same edge
                        p.FlipV();
//                        cout << "AFTER FlipV: "<< p.V()<<endl;
                        alpha = p.AngleRad();
                        sincos(alpha, &sin_alpha, &cos_alpha);
                        cot_alpha = sin_alpha / cos_alpha;
                        Q = p.V()->P();
                        p.FlipV();//return to the original vert

                        //go to the NEXT EDGE and the other vertex of the edge
                        p.FlipE();
//                        cout << "AFTER FlipE: "<<endl;
                        p.FlipV();
//                        cout << "AFTER FlipV: "<< p.V()<<endl;
                        R = p.V()->P();
                        edge_idx = map_vert_idx[p.V()];
//                        cout << "edge index: "<<edge_idx<<endl;
                        p.FlipV();
//                        cout << "should be the same of the starting one:"<< p.V()<<endl;
                        angle = p.AngleRad();
                        //go to the adjacent face
                        p.FlipF();
                        p.FlipE();

                        //go to the opposite vert on the same edge
                        p.FlipV();
//                        cout << "AFTER FlipF, FlipE and FlipV (should be different form the others): " << p.V()<<endl;
                        beta = p.AngleRad();
                        sincos(beta, &sin_beta, &cos_beta);
                        p.FlipV();
                        cot_beta = sin_beta / cos_beta;
                        w = (cot_alpha + cot_beta) * 0.5;

                        area = compute_area(P, Q, R, angle);
                        sum += w;
                        sum_area += area;
                        //inserting <i, j, w>
                        laplacian[vert_idx].emplace_back(
                                tuple<int, int, double>(make_tuple(vert_idx, edge_idx, w)));

                        //NEXT EDGE becomes CURRENT EDGE
                        p.FlipE();
//                        cout << p.f<<endl;
//                        cout << fi.base()<<endl;
//                        cout << "EDGE PROCESSED "<<endl;
                    } while (p.f != fi.base());
                    //inserting on the diagonal
                    laplacian[vert_idx].emplace_back(tuple<int, int, double>(make_tuple(vert_idx, vert_idx, -sum)));
                    laplacian_weights[vert_idx] = tuple<int, int, double>(
                            make_tuple(vert_idx, vert_idx, 1 / (sum_area / 3)));
                }
            }
        }
    }


    void print_lapacian(){
        cout << "laplacian size:"<<laplacian.size() <<endl;
        for(auto item = laplacian.begin(); item != laplacian.end(); item ++){
            cout << "triples size: "<< item->size()<<endl;
            for(auto triple = item->begin(); triple != item->end(); triple++){
                printf("row: %d col:%d value:%f\n", get<0>(*triple), get<1>(*triple), get<2>(*triple));
            }
        }
    }
};
#endif //MCF_SKET_LAPLACEHELPER_HPP
