//
// Created by boy on 01/11/22.
//

#ifndef MCF_SKET_MEANCURVATUREFLOW_HPP
#define MCF_SKET_MEANCURVATUREFLOW_HPP

#include "Collapser.h"
#include "LSMAT.hpp"
#include "LaplaceHelper.hpp"
#include "sphereShrinking.hpp"
#include "utils.hpp"

class MeanCurvatureFlow{
private:
    MyMesh *m;
    MyMesh::PerVertexAttributeHandle<int> vert_idx;
    MyMesh::PerVertexAttributeHandle<bool> isFixed;
    MyMesh::PerVertexAttributeHandle<bool> isSplitted;

public:
    MeanCurvatureFlow(MyMesh *m){
        this->m = m;

        vert_idx = tri::Allocator<MyMesh>::GetPerVertexAttribute<int>(*m, string("map_vert_idx"));
        isFixed = tri::Allocator<MyMesh>::GetPerVertexAttribute<bool>(*m, string("isFixed"));
        isSplitted = tri::Allocator<MyMesh>::GetPerVertexAttribute<bool>(*m, string("isSplitted"));
        int i = 0;
        for(auto vi = m->vert.begin(); vi != m->vert.end(); vi++, i++){
            vert_idx[vi] = i;
            isFixed[vi] = false;
            isSplitted[vi] = false;
        }
    }
    double SignedVolumeOfTriangle(Point3d p1, Point3d p2, Point3d p3){
        double v321 = p3.X()*p2.Y()*p1.Z();
        double v231 = p2.X()*p3.Y()*p1.Z();
        double v312 = p3.X()*p1.Y()*p2.Z();
        double v132 = p1.X()*p3.Y()*p2.Z();
        double v213 = p2.X()*p1.Y()*p3.Z();
        double v123 = p1.X()*p2.Y()*p3.Z();
        return (1.0f/6.0f)*(-v321 + v231 + v312 - v132 - v213 + v123);
    }

    double VolumeOfMesh(MyMesh *m){
        double sum = 0;
        for(auto f = m->face.begin(); f < m->face.end(); f++){
            sum += abs(SignedVolumeOfTriangle(f->V(0)->P(),f->V(1)->P(),f->V(2)->P()));
        }
        return sum;
    }

    void compute_skel(){
        VertexConstDataWrapper<MyMesh> ww(*m);
        KdTree<float> tree(ww);
        SphereShrinking ss = SphereShrinking(m, tree);
        ss.compute_ma_point();
        vector<Sphere3d> medial = ss.getMedialSpheres();

        MyMesh::PerVertexAttributeHandle<MyMesh::CoordType> vert_mat = tri::Allocator<MyMesh>::GetPerVertexAttribute<MyMesh::CoordType>(*m, string("vert_mat"));

        int i = 0;
        for (auto vi = m->vert.begin(); vi != m->vert.end(); vi++, i++){
            vert_mat[vi] = medial[i].Center();
        }

        cout << "SS done"<<endl;
        MyMesh *mesoSkel;
        LaplaceHelper laplaceHelper(m);
        Collapser collapser(m);
        int counter = 0;
        while(true) {
            laplaceHelper.compute_laplace();
            collapser.compute("skel_"+ to_string(counter++)+".off");
            cin.get();
        }
    }




};
#endif //MCF_SKET_MEANCURVATUREFLOW_HPP
