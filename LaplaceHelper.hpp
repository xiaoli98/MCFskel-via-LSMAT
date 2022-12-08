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
    typedef vector<tuple<MyMesh::VertexType*, MyMesh::VertexType*, double>> laplacian_triple;
private:
    MyMesh *m;
    MyMesh::VertContainer *vert;
    MyMesh::FaceContainer *faces;
    laplacian_triple laplacian;
    laplacian_triple laplacian_weights;

    MyMesh::PerVertexAttributeHandle<int> vert_idx;
    MyMesh::PerMeshAttributeHandle<laplacian_triple> laplacian_handle;
public:
    LaplaceHelper(MyMesh *m){
        this->m = m;
        vert = &this->m->vert;
        faces = &this->m->face;
        vert_idx = tri::Allocator<MyMesh>::GetPerVertexAttribute<int>(*m, string("map_vert_idx"));
    }

    const vector<tuple<MyMesh::VertexType*, MyMesh::VertexType*, double>> &getLaplacian() const {
        return laplacian;
    }

    const vector<tuple<MyMesh::VertexType*, MyMesh::VertexType*, double>> &getLaplacianWeights() const {
        return laplacian_weights;
    }

    double compute_area(Point3d p, Point3d q, Point3d r, double angle){
        double pq = Distance(p, q);
        double pr = Distance(p, r);
        return 0.5 * pq * pr * sin(angle);
    }

    void compute_laplace() {
        int n_vert = vert->size();
        laplacian.reserve(6 * n_vert);
        cout << "reserving: "<<6* n_vert<<endl;
//        laplacian_weights.resize(n_vert);
        tri::UpdateTopology<MyMesh>::VertexFace(*m);
        tri::UpdateTopology<MyMesh>::FaceFace(*m);
        tri::UpdateFlags<MyMesh>::VertexClearV(*m);
//        MyMesh::VertexIterator vi;
        MyMesh::FaceIterator fi;
        int visited_vert = 0;

        ///creating mapping of vertex and the index
        //todo -> risolto, la copia di vettori e' il male
        //  non ben capito ma i puntatori nei container di vertex sono diversi dai puntatori dei vertici nelle facce
        //  quindi inizializzo la mappa con usando i puntatori dei vertici delle facce
        /// TROPPO PYTHON HA ROVINATO LA NUOVA GENERAZIONE

        MyMesh::VertexType *vert0_p, *vert1_p;

        tri::UpdateFlags<MyMesh>::VertexClearV(*m);

        for (auto vit = vert->begin(); vit != vert->end(); vit++){
            MyFace* start = vit->VFp();
            vcg::face::Pos<MyFace> p(start, vit.base());
            double angle_sum = 0;
            Point3d p_sum = Point3d (0,0,0);
            double sin_alpha, cos_alpha, alpha, cot_alpha;
            double sin_beta, cos_beta, beta,cot_beta;
            double w = 0, sum_area=0, sum = 0;
            do
            {
                double area = 0, angle;
                Point3d P, Q, R; //3 vertices to be save for area computation
                P = p.V()->P();
                vert0_p = p.V();

                //IN THE PREV EDGE
                //go to the opposite vert on the same edge
                p.FlipV();
                alpha = p.AngleRad();
                sincos(alpha, &sin_alpha, &cos_alpha);
                cot_alpha = cos_alpha / sin_alpha;
                Q = p.V()->P();
                p.FlipV();//return to the original vert

                //go to the NEXT EDGE and the other vertex of the edge
                p.FlipE();
                p.FlipV();
                R = p.V()->P();
                vert1_p = p.V();

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
                cot_beta = cos_beta / sin_beta;
                w = (cot_alpha + cot_beta) * 0.5;

                area = compute_area(P, Q, R, angle);
                sum += w;
                sum_area += area;
                //inserting <i, j, w>
                laplacian.emplace_back(make_tuple(vert0_p, vert1_p, w));

                //NEXT EDGE becomes CURRENT EDGE
                p.FlipE();
            }while(p.f!=start);
            laplacian.emplace_back(make_tuple(vert0_p, vert0_p, -sum));
            laplacian_weights.emplace_back(make_tuple(vert0_p, vert0_p, 1 / (sum_area / 3)));
        }

        laplacian_handle = tri::Allocator<MyMesh>::GetPerMeshAttribute<laplacian_triple>(*m, string("laplacian"));
    }

    void print_lapacian(){
        cout << "laplacian size:"<<laplacian.size() <<endl;
        for(auto item = laplacian.begin(); item != laplacian.end(); item++){
            printf("row: %d col:%d value:%f\n", vert_idx[get<0>(*item)], vert_idx[get<1>(*item)], get<2>(*item));
        }
    }
};
#endif //MCF_SKET_LAPLACEHELPER_HPP
