//
// Created by boy on 22/09/22.
//
#ifndef LSMAT_HPP
#define LSMAT_HPP

#define DEBUG_LSMAT 0

#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <random>
#include <string>

#include "utils.hpp"

using namespace vcg;
using namespace std;
//TODO
// come inizializzare
// cosa e' s^t-1 quando sono all'iterazione iniziale?
// instabilita' numerica


//TODO
//  implemantare la versione "classica"
//  debugging visuale
// vedere se le sfere sono corrette
// gnerare punti e trovare l'energia minima nella bounding box


//TODO
// inizializzare bene anche sfere vecchi 
class LSMAT{
public:
    LSMAT(MyMesh *m){
        this->m = m;
        tri::UpdateBounding<MyMesh>::Box(*m);
        tri::RequirePerVertexNormal(*m);
        tri::UpdateNormal<MyMesh>::PerVertexNormalized(*m);
        for(int i = 0; i < m->VN(); i++){
            point_list.emplace_back(m->vert[i].P());
            normal_list.emplace_back(m->vert[i].N());
        }
//        random_device rd;
//        mt19937 gen(rd());
//        uniform_real_distribution<> dis(0,1.0);

        double r = 0;
        Point3d center;
//        double offset = m->bbox.Diag()/10.0f;
//        cout <<"offset:" <<offset<<endl;
        for(int i = 0; i < m->VN(); i++){
//            r = ((double) rand() / (RAND_MAX)) + 1;
            r = m->bbox.Diag()/10;
            center = point_list[i].operator-(normal_list[i].operator*(r));
            this->spheres.emplace_back(Sphere3d(center, (r)));
            if(i == 0){
//                add_sphere(*m, center, r);
            }
        }
        this->spheres_old = vector<Sphere3d>(m->VN(), Sphere3d(Point3d(0,0,0), m->bbox.Diag()));
#if DEBUG_LSMAT
        cout << "***********INITIAL POINT_LIST**********"<<endl;
        for(auto & i : point_list){
            PRINTP(i)
        }

        cout << "***********INITIAL SPHERES**********"<<endl;
        for(auto & sphere : spheres){
            cout << "center: ";
            PRINTP(sphere.Center())
            cout << "radius: "<<sphere.Radius()<<endl;
        }
#endif
        cout <<"========================================================================================"<<endl;
    }


private:
    MyMesh *m;

    vector<Point3d> point_list;
    vector<Point3d> normal_list;
    vector<Sphere3d> spheres;
    vector<Sphere3d> spheres_old;
    double sigma = 0;
    double omega_1 = 0.001;     //modifying one will update the other omega_*
    double omega_2 = omega_1 / (0.007 * sigma + 0.02);
    double epsilon = 0.01;
    double h_blend = 0.74 * sigma + 0.49;
    double h_support = 0.74 * sigma + 0.49;
    double d_in = 0.75 * sigma;

    static double ramp(double x){
        if (x > 0) return x;
        else return 0;
    }
    static double b(double x){
        if (x < 1){
            return pow(1-pow(x,2), 4);
        }
        else return 0;
    }
    static double phi(Point3d p, double h){
        double x = p.Norm();
        return b(x/h);
    }
    static double phi(double x, double h){
        return b(x/h);
    }
    static double mix(double a, double b, double x){
        return x * a + (1 - x) * b;
    }
    double e_medial(Sphere3d s, Sphere3d s_old){
        double _e_maximal = e_maximal(s, s_old);
        double _e_inscribed = e_inscribed(s, s_old);
#if DEBUG_LSMAT
        cout <<"%%%%%%%%%%%%%% E_MEDIAL %%%%%%%%%%%%%%"<<endl;
        cout << "%\tomega_1: "<<omega_1<<endl;
        cout << "%\tomega_2: "<<omega_2<<endl;
        cout << "%\te_maximal: " << _e_maximal<<endl;
        cout << "%\te_inscribed: " << _e_inscribed<<endl;
        cout <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
#endif
        double result = omega_1*_e_maximal + omega_2*_e_inscribed;
        return result;
    }
    double e_maximal(Sphere3d s, Sphere3d s_old){
        double energy = s.Radius() - (s_old.Radius() + epsilon);
        return sqrt(pow(energy, 2));
    }
    double e_inscribed(Sphere3d s, Sphere3d s_old){
        double sum = 0;
        for (int i = 0; i < point_list.size(); i++){
            double _ramp = ramp(Distance(s_old.Center(), point_list[i]) - s_old.Radius());
            double _phi = phi(_ramp, h_support);
            Point3d projection = projection_c(s_old.Center(), point_list[i], normal_list[i]);
            sum += _phi * phi_blend_square(s, s_old, point_list[i], normal_list[i], projection);
        }
        return sum;
    }
    Eigen::VectorXd grad_e_inscribed(Sphere3d s, Sphere3d s_old){
        Eigen::Vector4d gradient;
        double _ramp, _phi1, _phi2;
        Point3d projection;

        double distancePC;
        for(int i = 0; i < point_list.size(); i++){
            _ramp = ramp(Distance(s_old.Center(), point_list[i]) - s_old.Radius());
            _phi1 = phi(_ramp, h_support);
            projection = projection_c(s_old.Center(), point_list[i], normal_list[i]);
            _phi2 = phi(projection, h_blend);

            Eigen::Vector4d gradient_phi_plane;
            double _phi_plane = phi_plane(s, point_list[i], normal_list[i]);
            if(_phi_plane > 0){
                gradient_phi_plane << 2*(-normal_list[i].X())*(s.Radius()-(point_list[i].operator-(s.Center()).dot(normal_list[i]))),
                        2*(-normal_list[i].Y())*(s.Radius()-(point_list[i].operator-(s.Center()).dot(normal_list[i]))),
                        2*(-normal_list[i].Z())*(s.Radius()-(point_list[i].operator-(s.Center()).dot(normal_list[i]))),
                        2*(s.Radius()-(point_list[i].operator-(s.Center()).dot(normal_list[i])));
            }
            else{
                gradient_phi_plane << 0, 0, 0, 0;
            }
            distancePC = point_list[i].operator-(s.Center()).Norm();
            Eigen::Vector4d gradient_phi_point;
            double _phi_point = phi_point(s, point_list[i]);
            if(_phi_point > 0 ){
                gradient_phi_point <<   2*((point_list[i].X()-s.Center().X())/distancePC)*(s.Radius() - distancePC),
                        2*((point_list[i].Y()-s.Center().Y())/distancePC)*(s.Radius() - distancePC),
                        2*((point_list[i].Z()-s.Center().Z())/distancePC)*(s.Radius() - distancePC),
                        2*(s.Radius() - distancePC);
            }
            else{
                gradient_phi_point << 0, 0, 0, 0;
            }
            gradient  += 2* ( (_phi2-1) * _phi_point * gradient_phi_point - _phi2 * _phi_plane * gradient_phi_plane);
        }
        return gradient;
    }

    double e_pinning(Sphere3d s, Point3d p){
        double _ramp = ramp(s.Center().operator-(p).Norm() - (s.Radius() + d_in));
        return pow(_ramp, 2);
    }
    double phi_plane(Sphere3d s, Point3d p, Point3d n){
        //RAMP(r - (p_n - c ) * n )
        return ramp(s.Radius() - (p.operator-(s.Center())).dot(n));
    }
    double phi_point(Sphere3d s, Point3d p){
        //RAMP(r-|| p_n - c ||_2)
//        Point3d temp = p.operator-(s.Center());
        return ramp(s.Radius() - Distance(p, s.Center()));
    }

    double phi_blend_square(Sphere3d s, Sphere3d s_old, Point3d p, Point3d n, Point3d projection){
        return mix(pow(phi_plane(s, p, n),2), pow(phi_point(s, p),2), phi( projection.operator-(p), h_blend));
    }

    Point3d projection_c(Point3d c, Point3d p, Point3d n){
        return c.operator-(n.operator*(c.operator-(p).operator*(n)));
    }


public:
#pragma region getter/setter
    double getSigma() const {
        return sigma;
    }

    void setSigma(double sigma) {
        LSMAT::sigma = sigma;
    }

    double getOmega1() const {
        return omega_1;
    }

    void setOmega1(double omega1) {
        omega_1 = omega1;
        setOmega2(omega_1/(0.007 * sigma + 0.02));
    }

    double getOmega2() const {
        return omega_2;
    }

    void setOmega2(double omega2) {
        omega_2 = omega2;
        setOmega1((0.007 * sigma + 0.02) * omega_2);
    }

    double getEpsilon() const {
        return epsilon;
    }

    void setEpsilon(double epsilon) {
        LSMAT::epsilon = epsilon;
    }

    double getHBlend() const {
        return h_blend;
    }

    void setHBlend(double hBlend) {
        h_blend = hBlend;
    }

    double getHSupport() const {
        return h_support;
    }

    void setHSupport(double hSupport) {
        h_support = hSupport;
    }

    double getDIn() const {
        return d_in;
    }

    void setDIn(double dIn) {
        d_in = dIn;
    }

#pragma endregion

    vector<Sphere3d> compute(int n_runs, double threshold = 0.01){
        double energy;
        double energy_old = 0;
        double residual;
        int j = 0;
//        vector<Sphere3d> sphere_list;
//        sphere_list.resize(point_list.size());
        Sphere3d newSphere;
        Point3d newPoint;
        Eigen::Vector4d vec;


        for(int i = 0; i<point_list.size(); i++){
            while(j++ < n_runs || residual > threshold){
#if DEBUG_LSMAT
                cout << "################# POINT "<<i<<" #################"<<endl;
#endif
                tuple<Eigen::VectorXd, double, double> t = gauss_newton(point_list[i], normal_list[i], spheres[i], spheres_old[i], energy_old);
                energy_old = energy;
                energy = get<1>(t);
                residual = get<2>(t);
                spheres_old[i] = spheres[i];
                vec = get<0>(t);
                newPoint.X() = vec[0];
                newPoint.Y() = vec[1];
                newPoint.Z() = vec[2];
                newSphere.Center() = newPoint;
                newSphere.Radius() = vec[3];
#if DEBUG_LSMAT
//                this->m->vert[i].C() = Color4b::Red;
//                tri::Allocator<MyMesh>::AddVertex(*m, newSphere.Center(), Color4b::Red);
//                tri::Allocator<MyMesh>::AddFace(*m, m->vert[i].P(), spheres[i].Center(), spheres_old[i].Center());
//                tri::Append<MyMesh,MyMesh>::Mesh(*m, *uniform_sphere_points(newSphere.Center(), newSphere.Radius()));
                string filename = "LSMAT_it_" + to_string(j) + "_point_" + to_string(i)+".off";
//                add_octahedron(*m, spheres[i].Center(), spheres[i].Radius(), filename);
//                add_octahedron(*m, newSphere.Center(), newSphere.Radius(),filename);
//                tri::io::ExporterOFF<MyMesh>::Save(*m, filename.c_str(), tri::io::Mask::IOM_FACECOLOR);
                cin.get();
#endif
                spheres[i] = newSphere;
            }
//            spheres_old = spheres;
        }
        return spheres;
    }

    ///
    /// \param p1 i-th point
    /// \param n1 normal associated with p1
    /// \param s1 sphere associated with p1
    /// \param s1_old old sphere associated with p1
    /// \param energy_old  old energy
    /// \return 1 - X_{n+1},
    ///         2 - new energy
    ///         3 - residual
    tuple<Eigen::VectorXd, double, double> gauss_newton(Point3d p1, Point3d n1, Sphere3d s1, Sphere3d s1_old, double energy_old){
        Point3d p = p1;
        Point3d n = n1;
        Sphere3d s = s1;
        Sphere3d s_old = s1_old;

        Eigen::Vector4d X_n (s.Center().X(), s.Center().Y(), s.Center().Z(), s.Radius());

        //gradient phi_blend_square
        Eigen::Vector4d gradient_e_inscribed = grad_e_inscribed(s, s_old);

        Eigen::Vector4d gradient_maximal(0,0,0,2*(s.Radius()-s_old.Radius()-epsilon));

        double distanceCP = Distance(s.Center(), p);
        Eigen::Vector4d gradient_pinning;
        if(e_pinning(s, p) > 0){
            gradient_pinning << 2*((s.Center().X()-p.X())/distanceCP) * (distanceCP - s.Radius() - d_in),
                                2*((s.Center().Y()-p.Y())/distanceCP) * (distanceCP - s.Radius() - d_in),
                                2*((s.Center().Z()-p.Z())/distanceCP) * (distanceCP - s.Radius() - d_in),
                                -2*(distanceCP - s.Radius() - d_in);
        }
        else{
            gradient_pinning << 0, 0, 0, 0;
        }
        Eigen::Vector4d jacobian = omega_1 * gradient_maximal + omega_2 * gradient_e_inscribed + gradient_pinning;
        Eigen::Matrix<double, 4, 4> matrix = jacobian * jacobian.transpose();
        Eigen::Vector4d vec = matrix.inverse() * jacobian;
#if DEBUG_LSMAT
        cout << "X_n: "<<endl << X_n<<endl;
//        cout << "phi_plane: " << phi_plane(s, p, n)<<endl;
//        cout << "gradient_phi_plane: " <<endl<< gradient_phi_plane<<endl;
//        cout << "distance PC: " << distancePC<<endl;
//        cout << "phi_point: " << phi_point(s, p)<<endl;
//        cout << "gradient_phi_point: " <<endl<< gradient_phi_point<<endl;
        cout << "gradient_e_inscribed: "<< endl << gradient_e_inscribed<<endl;
        cout << "gradient_maximal: " <<endl<< gradient_maximal<<endl;
        cout << "distance CP: " << distanceCP<<endl;
        cout << "e_pinning: "<< e_pinning(s,p)<<endl;
        cout << "gradient_pinning: " <<endl<<gradient_pinning <<endl;
        cout << "jacobian : " << endl<<jacobian<<endl;
        cout << "J^T J: " <<endl<<matrix <<endl;
        cout << "(J^T J)^-1"<<endl<<matrix.inverse()<<endl;
        cout << "(J^T J)^-1 J^T " << endl<<vec<<endl;
#endif

        // TODO
        // entries of Jacobian are the partial derivative of the residual on s
        // the problem should be linearized, but how?
        s_old = s;
        s.Center() = Point3d(vec[0], vec[1], vec[2]);
        s.Radius() = vec[3];
        double energy = e_medial(s, s_old);
        double residual = energy;

        vec = vec * residual;
#if DEBUG_LSMAT
//        cout << "p: ";
//        PRINTP(p)
//        cout << "n: ";
//        PRINTP(n)
//        cout << "sphere center: ";
//        PRINTP(s.Center());
//        cout <<"sphere radius: "<<s.Radius()<<endl;
//        cout << "sphere_old center: ";
//        PRINTP(s_old.Center());
//        cout <<"sphere_old radius: "<<s_old.Radius()<<endl;
        cout << "energy: "<<energy<<endl;
        cout <<"residual: " << residual<<endl;
        cout << "vec = "<<endl<<vec<<endl;
        cout <<"X_n - vec = "<< endl<<X_n - vec<<endl;
        cout <<"residual: "<<residual<<endl;

#endif
        return tuple<Eigen::VectorXd, double, double> {X_n - vec, energy,residual};
    }
};


#endif //LSMAT_HPP