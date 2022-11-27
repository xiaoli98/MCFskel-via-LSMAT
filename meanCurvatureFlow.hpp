//
// Created by boy on 01/11/22.
//

#ifndef MCF_SKET_MEANCURVATUREFLOW_HPP
#define MCF_SKET_MEANCURVATUREFLOW_HPP

#include "LSMAT.hpp"
#include "sphereShrinking.hpp"
#include "utils.hpp"

class MeanCurvatureFlow{
private:
    MyMesh *m;

public:
    MeanCurvatureFlow(MyMesh *m){
        this->m = m;
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
    }

    void compute_skel(){
        VertexConstDataWrapper<MyMesh> ww(*m);
        KdTree<float> tree(ww);
        SphereShrinking ss = SphereShrinking(m, tree);
        ss.compute_ma_point();
        vector<Sphere3d> medial = ss.getMedialSpheres();

        MyMesh *mesoSkel;
        while(true) {
            if (VolumeOfMesh(mesoSkel) > 1) {
                //  update laplacian and weights
                    //v^{t+1}_i = v^t_i + h * Delta(v^t_i)
                        //Delta(v^t_i) = -d * H(v)_p * n(p) - H(p)* \hat(D)_v * n_p
                    // solve Eq. 4
                    // which minimize E = norm(LV^{t+1})^2 + w^2_H * \sum_i(norm(v^{t+1}_i - v^t_i)^2)
                            // or better E = E_smooth + E_velocity + E_medial
                            //
                //  perform mesh contraction
                //  update connectivity
                //  collapse the shortest edges
            }
        }
    }




};
#endif //MCF_SKET_MEANCURVATUREFLOW_HPP
