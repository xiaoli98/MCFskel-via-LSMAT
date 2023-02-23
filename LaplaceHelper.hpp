//
// Created by boy on 27/11/22.
//

#ifndef MCF_SKET_LAPLACEHELPER_HPP
#define MCF_SKET_LAPLACEHELPER_HPP

#include <unordered_map>
#include <math.h>
#include <tuple>
#include <vcg/complex/algorithms/smooth.h>
#include "utils.hpp"


class LaplaceHelper{
    typedef vector<tuple<MyMesh::VertexType*, MyMesh::VertexType*, double>> laplacian_triple;
private:
    MyMesh *m;
    laplacian_triple laplacian;
    tri::Smooth<MyMesh> smoother;

    MyMesh::PerMeshAttributeHandle<laplacian_triple> laplacian_handle;
public:
    LaplaceHelper(MyMesh *m){
        this->m = m;
    }

    const vector<tuple<MyMesh::VertexType*, MyMesh::VertexType*, double>> &getLaplacian() const {
        return laplacian;
    }

//    double compute_area(Point3d p, Point3d q, Point3d r, double angle){
//        double pq = Distance(p, q);
//        double pr = Distance(p, r);
//        return 0.5 * pq * pr * sin(angle);
//    }

    struct EdgeComp{
        bool operator()(array<MyMesh::VertexPointer,2> e1, array<MyMesh::VertexPointer,2> e2){
            if(e1[0] != e2[0]) return e1[0] < e2[0];
            else return e1[1] < e2[1];
        }
    };

    typedef map<array<MyMesh::VertexPointer,2>, double, EdgeComp> edgeWeight;
    edgeWeight laplacian_weights(){
        edgeWeight allEdges;
        double angle = 0, weight;
        for(auto fit = m->face.begin(); fit != m->face.end(); fit++){
            for(int i = 0; i < fit->VN(); i++) {
                array<MyMesh::VertexPointer, 2> vertPair = {fit->V0(i), fit->V1(i)}; //aka edge
                if(vertPair[0] > vertPair[1]) swap(vertPair[0], vertPair[1]);
                MyMesh::VertexPointer oppositeVert = fit->V2(i);
                angle = Angle(vertPair[0]->cP() - oppositeVert->cP(), vertPair[1]->cP() - oppositeVert->cP());
                weight = tan((M_PI * 0.5) - angle);
                if(weight > 100000) {
                    printf("strange weight\n");
                    printf("angle: %f\n", angle);
                    Point3d p1 = vertPair[0]->cP() - oppositeVert->cP();
                    Point3d p2 = vertPair[1]->cP() - oppositeVert->cP();
                    PRINTP(p1)
                    PRINTP(p2)
                    printf("norm mult: %f\n", p1.Norm() * p2.Norm());
                    double t = (p1*p2)/p1.Norm()*p2.Norm();
                    printf("t: %f\n", t);
                    if(t > 1){
                        printf("not good t\t");
                        printf("acos(t): %f\n", acos(1));
                    }
                    else {
                        printf("acos(t): %f\n", acos(t));
                    }
                    printf("ciao\n");

                }
                if(allEdges.find(vertPair) != allEdges.end()) {
                    allEdges.insert(make_pair(vertPair, weight));
                }
                else{
                    allEdges[vertPair] += weight;
                }
            }
        }
        return allEdges;
    }

#define MY_LAPLACIAN 0

    typedef tri::Smooth<MyMesh>::LaplacianInfo LaplacianInfo;
    void compute_laplace() {
        int n_vert = m->vert.size();
        laplacian.clear();
        laplacian.reserve(6 * n_vert);
        tri::UpdateTopology<MyMesh>::VertexFace(*m);
#if MY_LAPLACIAN == 0
        LaplacianInfo lpz(MyMesh::CoordType (0, 0, 0), 0);
        SimpleTempData<typename MyMesh::VertContainer, LaplacianInfo> TD(m->vert, lpz);
        smoother.AccumulateLaplacianInfo(*m, TD, true);
        vector<MyMesh::VertexPointer> adjacentVertices;
        double w, sum_w;
        for(auto vit = m->vert.begin(); vit != m->vert.end(); vit++){
            sum_w = 0;
            vcg::face::VVStarVF<MyMesh::FaceType>(vit.base(), adjacentVertices);
            for(auto adjv = adjacentVertices.begin(); adjv != adjacentVertices.end(); adjv++){
//                i = tri::Index(*m, vit.base());
//                j = tri::Index(*m, *(adjv.base()));
//                assert(TD[vit.base()].cnt >= 0);
                w = TD[vit.base()].cnt < 2 ? 1 : TD[vit.base()].cnt;
                //todo
                // laplacian weights are strange
                if(w > 10000) {
                    printf("WARNING: w value too high (%f)\n", w);
                }
                sum_w += w;
                laplacian.emplace_back(make_tuple(vit.base(), *(adjv.base()), w));
            }
            laplacian.emplace_back(make_tuple(vit.base(), vit.base(), -sum_w));
        }
#else
        edgeWeight allEdges = laplacian_weights();
        map<MyMesh::VertexPointer, double> weight_sum;
        for(auto v = m->vert.begin(); v != m->vert.end(); v++)
            weight_sum[v.base()] = 0;

        for(auto const& item: allEdges){
            array<MyMesh::VertexPointer,2> edge = item.first;
            double weight = item.second;
            weight_sum[edge[0]] -= weight;
            weight_sum[edge[1]] -= weight;
            laplacian.emplace_back(make_tuple(edge[0], edge[1], weight));
        }

        //first is the vertex pointer
        //second is the sum of weights
        for(auto const& item: weight_sum)
            laplacian.emplace_back(make_tuple(item.first, item.first, item.second));
#endif
        laplacian_handle = tri::Allocator<MyMesh>::GetPerMeshAttribute<laplacian_triple>(*m, string("laplacian"));
        laplacian_handle() = laplacian;

    }

    void print_lapacian(){
        cout << "laplacian size:"<<laplacian.size() <<endl;
        for(auto item = laplacian.begin(); item != laplacian.end(); item++){
            printf("row: %lu col:%lu value:%f\n", tri::Index(*m,get<0>(*item)), tri::Index(*m, get<1>(*item)), get<2>(*item));
        }
    }
};
#endif //MCF_SKET_LAPLACEHELPER_HPP
