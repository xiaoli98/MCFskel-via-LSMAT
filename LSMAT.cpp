//
// Created by boy on 22/09/22.
//
#include <iostream>
#include <math.h>

#include <vcg/complex/complex.h>
#include <vcg/simplex/vertex/component.h>
#include <Eigen/Dense>
#include<wrap/io_trimesh/import_obj.h>
#include<wrap/io_trimesh/export_obj.h>
#include<wrap/io_trimesh/import_off.h>
#include<wrap/io_trimesh/export_off.h>

#define PRINTP(p) cout << p.X() <<" "<< p.Y() <<" "<< p.Z() <<endl;

class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public vcg::UsedTypes<vcg::Use<MyVertex>   ::AsVertexType,
        vcg::Use<MyEdge>     ::AsEdgeType,
        vcg::Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::VFAdj, vcg::vertex::BitFlags, vcg::vertex::Color4b>{};
class MyFace    : public vcg::Face<   MyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef, vcg::face::Color4b, vcg::face::BitFlags > {};
class MyEdge    : public vcg::Edge<   MyUsedTypes> {};

class MyMesh    : public vcg::tri::TriMesh<std::vector<MyVertex>, std::vector<MyFace> , std::vector<MyEdge>  > {};

using namespace vcg;
using namespace std;
//TODO
//sphere shringking algorithm
//  dato un mesh watertight con i suoi punti {(p_n, n_n)}_n
//  un punto P sul bordo del mesh deve soddisfare 2 proprieta':
//      1) la sfera passante per P deve essere vuota
//      2) il vettore distanza tra P e centro C e' parallela alla norma n del punto P
//  iniziamo una sfera grande che passa per P, la sfera deve contenere tutto il set di punti del mesh,
//  e il centro C deve stare lungo il segmente (p, -n)
//  se la distanza tra C e un punto F del mesh e' piu' corta del raggio R,
//  allora il centro e' ricalcolato con la tangente a (P, n) e passante per F.


//TODO
//optimization
//

class LSMAT{
private:

    float sigma;
    vector<Point3f*> point_list;
    vector<Point3f*> normal_list;
    vector<Sphere3f> sphere;
    vector<Sphere3f> sphere_old;

    LSMAT(MyMesh *m){
        tri::UpdateBounding<MyMesh>::Box(*m);
        tri::RequirePerVertexNormal(*m);
        tri::UpdateNormal<MyMesh>::PerVertexNormalized(*m);
        for(int i = 0; i < m->VN(); i++){
            point_list.emplace_back((Point3f*)&(m->vert[i].P()));
            normal_list.emplace_back((Point3f*)&(m->vert[i].N()));
        }
        this->sphere = vector<Sphere3f>(m->VN(), Sphere3f(m->bbox.Center(),m->bbox.Diag()));
        this->sphere_old.resize(m->VN());
    }

    static float ramp(float x){
        if (x > 0) return x;
        else return 0;
    }

    float e_medial(float omega_1, float omega_2, float e_maximal, float e_inscribed ){
        return omega_1*e_maximal + omega_2*e_inscribed;
    }

    float e_maximal(float radius, float radius_old, float epsilon){
        float energy = radius - (radius_old + epsilon);
        return sqrt(energy);

    }
    //gradient of e_inscribed is the derivative of Phi_blend()
    float e_inscribed(Sphere3f s, Sphere3f s_old, Point3f p, Point3f n, float h_blend, float h_support){
        float sum = 0;
        for (int i = 0; i < point_list.size(); i++){
            float _ramp = ramp(s_old.Center().operator-(p).Norm() - s_old.Radius());
            float _phi = phi(_ramp, h_support);
            sum += _phi;
        }
        return sum * pow(phi_blend(s, s_old, p, n, h_blend),2);
    }

    float e_pinning(Sphere3f s, Point3f p, float d_in){
        float _ramp = ramp(s.Center().operator-(p).Norm() - (s.Radius() + d_in));
        return pow(_ramp, 2);
    }

    float b(float x){
        if (x < 1){
            return pow(1-pow(x,2), 4);
        }
        else return 0;
    }

    float phi(Point3f p, float h){
        float x = p.Norm();
        return b(x/h);
    }
    float phi(float x, float h){
        return b(x/h);
    }

    float phi_plane(Sphere3f s, Point3f p, Point3f n){
        //RAMP(r - (p_n - c ) * n )
        return ramp(s.Radius() - (p.operator-(s.Center())).dot(n));
    }

    float phi_point(Sphere3f s, Point3f p){
        //RAMP(r-|| p_n - c ||_2)
        Point3f temp = p.operator-(s.Center());

        double sum_square = 0;
        for (int i= 0; i < 3; i++){
            sum_square += pow(temp.V(i), 2);
        }
        return ramp(s.Radius() - sqrt(sum_square));
    }

    float mix(float a, float b, float x){
        return x * a + (1 - x) * b;
    }

    float phi_blend(Sphere3f s, Sphere3f s_old, Point3f p, Point3f n, float h_blend){
        //c^(t-1) - n(c^(t-1)-p) * n
        Point3f projection = s_old.Center().operator-(n.operator*(s_old.Center().operator-(p).operator*(n)));
        return mix(pow(phi_plane(s, p, n),2), pow(phi_point(s, p),2), phi( projection.operator-(p), h_blend));
    }

public:
    pair<Point3f, float> compute(){
        //TODO trasform
    }

    void gauss_newton(Point3f p1, Point3f n1, Sphere3f s1, Sphere3f s1_old){
        Point3f p = p1;
        Point3f n = n1;
        Sphere3f s = s1;
        Sphere3f s_old = s1_old;

        Point3f new_centre;
        float new_radius;

        Eigen::Vector4f gradient_phi_plane(-n.X(), -n.Y(), -n.Z(), 1);
        float distancePC = p.operator-(s.Center()).Norm();
        Eigen::Vector4f gradient_phi_point((p.X()-s.Center().X())/distancePC,
                                           (p.Y()-s.Center().Y())/distancePC,
                                           (p.Z()-s.Center().Z())/distancePC,
                                           1);
        Eigen::Vector4f gradient_maximal(0,0,0,1);
        float distanceCP = s.Center().operator-(p).Norm();
        Eigen::Vector4f gradient_pinning((s.Center().X()-p.X())/distanceCP,
                                         (s.Center().Y()-p.Y())/distanceCP,
                                         (s.Center().Z()-p.Z())/distanceCP,
                                         -1);
        Eigen::Vector4f jacobian(gradient_maximal+2*gradient_phi_plane+2*gradient_phi_point+gradient_pinning);
        Eigen::Matrix<float, 4, 4> matrix = jacobian.transpose()*jacobian;

    }
};



//int main(int argc, char* argv[]){
//    MyMesh m;
//
//    tri::io::ImporterOFF<MyMesh>::Open(m,argv[1]);
//    tri::UpdateBounding<MyMesh>::Box(m);
//    tri::RequirePerVertexNormal(m);
//    tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
//
//    vector<Point3f> points;
//    vector<Point3f> normals;
//    for(int i = 0; i < m.VN(); i++){
//        points.emplace_back(m.vert[i].P());
//        normals.emplace_back(m.vert[i].N());
//    }
//
//
//    Point3f p;
//
//
//}