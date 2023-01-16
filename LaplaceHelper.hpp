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
//    MyMesh::VertContainer *vert;
//    MyMesh::FaceContainer *faces;
    laplacian_triple laplacian;
    laplacian_triple laplacian_weights;
    tri::Smooth<MyMesh> smoother;

    MyMesh::PerMeshAttributeHandle<laplacian_triple> laplacian_handle;
public:
    LaplaceHelper(MyMesh *m){
        this->m = m;
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

    typedef tri::Smooth<MyMesh>::LaplacianInfo LaplacianInfo;
    void compute_laplace() {
        int n_vert = m->vert.size();
        laplacian.clear();
        laplacian.reserve(6 * n_vert);
        tri::UpdateTopology<MyMesh>::VertexFace(*m);
        
        LaplacianInfo lpz(MyMesh::CoordType (0, 0, 0), 0);
        SimpleTempData<typename MyMesh::VertContainer, LaplacianInfo> TD(m->vert, lpz);
        smoother.AccumulateLaplacianInfo(*m, TD, true);
        vector<MyMesh::VertexPointer> adjacentVertices;
        size_t i, j;
        double w, sum_w;
        for(auto vit = m->vert.begin(); vit != m->vert.end(); vit++){
            sum_w = 0;
            vcg::face::VVStarVF<MyMesh::FaceType>(vit.base(), adjacentVertices);
            for(auto adjv = adjacentVertices.begin(); adjv != adjacentVertices.end(); adjv++){
//                i = tri::Index(*m, vit.base());
//                j = tri::Index(*m, *(adjv.base()));
                w = TD[vit.base()].cnt * 0.5;
                sum_w += w;
                laplacian.emplace_back(make_tuple(vit.base(), *(adjv.base()), w));
            }
            laplacian.emplace_back(make_tuple(vit.base(), vit.base(), -sum_w));
        }

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
