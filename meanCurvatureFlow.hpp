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
    MyMesh m;
public:
    MeanCurvatureFlow(){
        VolumeOfMesh(m)
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

    double VolumeOfMesh(MyMesh m){
        double sum = 0;
        for(auto f = m.face.begin(); f < m.face.end(); f++){
            sum += abs(SignedVolumeOfTriangle(f->V(0)->P(),f->V(1)->P(),f->V(2)->P()));
        }
    }

    void compute_skel(){

    }

};
#endif //MCF_SKET_MEANCURVATUREFLOW_HPP
