//
// Created by boy on 01/11/22.
//

#ifndef MCF_SKET_MEANCURVATUREFLOW_HPP
#define MCF_SKET_MEANCURVATUREFLOW_HPP

#include "Collapser.h"
#include "LSMAT.hpp"
#include "LaplaceHelper.hpp"
#include "sphereShrinking.hpp"
#include "Voronoi.hpp"
#include "utils.hpp"

class MeanCurvatureFlow{
private:
    MyMesh *m;
    MyMesh::PerVertexAttributeHandle<bool> isFixed;
    MyMesh::PerVertexAttributeHandle<bool> isSplitted;

public:
    MeanCurvatureFlow(MyMesh *m){
        this->m = m;

        isFixed = tri::Allocator<MyMesh>::GetPerVertexAttribute<bool>(*m, string("isFixed"));
        isSplitted = tri::Allocator<MyMesh>::GetPerVertexAttribute<bool>(*m, string("isSplitted"));
        int i = 0;
        for(auto vit = m->vert.begin(); vit != m->vert.end(); vit++, i++){
            isFixed[*vit] = false;
            isSplitted[*vit] = false;
        }
    }

//    double SignedVolumeOfTriangle(Point3d p1, Point3d p2, Point3d p3){
//        double v321 = p3.X()*p2.Y()*p1.Z();
//        double v231 = p2.X()*p3.Y()*p1.Z();
//        double v312 = p3.X()*p1.Y()*p2.Z();
//        double v132 = p1.X()*p3.Y()*p2.Z();
//        double v213 = p2.X()*p1.Y()*p3.Z();
//        double v123 = p1.X()*p2.Y()*p3.Z();
//        return (1.0f/6.0f)*(-v321 + v231 + v312 - v132 - v213 + v123);
//    }
//
//    double VolumeOfMesh(){
//        double sum = 0;
//        for(auto f = m->face.begin(); f < m->face.end(); f++){
//            sum += abs(SignedVolumeOfTriangle(f->V(0)->P(),f->V(1)->P(),f->V(2)->P()));
//        }
//        return sum;
//    }
#define VORONOI 0
    void compute_skel(){
#if VORONOI
        Voronoi_MAT vmat(m);
        vmat.computeVoronoiDiagram();
        vmat.searchVoronoiPoles();
        vmat.getMedialSpokeAngleAndRadii();
        vmat.setMedialAttribute();
        MyMesh::PerVertexAttributeHandle<MyMesh::CoordType> vert_mat = tri::Allocator<MyMesh>::GetPerVertexAttribute<MyMesh::CoordType>(*m, string("vert_mat"));
#else

        VertexConstDataWrapper<MyMesh> ww(*m);
        KdTree<double> tree(ww);
        SphereShrinking ss = SphereShrinking(m, tree);
        ss.compute_ma_point();
        vector<Sphere3d> medial = ss.getMedialSpheres();
        int i = 0;
        MyMesh::PerVertexAttributeHandle<MyMesh::CoordType> vert_mat = tri::Allocator<MyMesh>::GetPerVertexAttribute<MyMesh::CoordType>(*m, string("vert_mat"));
        for (auto vi = m->vert.begin(); vi != m->vert.end(); vi++, i++){
            vert_mat[vi] = medial[i].Center();
        }
#endif
        cout << "SS done"<<endl;

        LaplaceHelper laplaceHelper(m);
        Collapser collapser(m);

        int counter = 0;
        while(counter < 15) {
            cout << "-------------------------------------------------"<<counter<<"-------------------------------------------------"<<endl;
//            printf("bbox_diag: %f\n", m->bbox.Diag());
//            tri::UpdateTopology<MyMesh>::VertexFace(*m);
//            tri::UpdateTopology<MyMesh>::FaceFace(*m);
//            tri::Allocator<MyMesh>::CompactEveryVector(*m);
            laplaceHelper.compute_laplace();
//            laplaceHelper.print_lapacian();
            collapser.compute();

            cout << "FN: "<< m->FN() << " VN: "<<m->VN()<<endl;
            cout << "meso skel created"<< ": skel_"+to_string(counter)+".off" <<endl;
            tri::io::ExporterOFF<MyMesh>::Save(*m, ("skel_"+to_string(counter++)+".off").c_str(), tri::io::Mask::IOM_FACECOLOR);
        }
    }

};
#endif //MCF_SKET_MEANCURVATUREFLOW_HPP
