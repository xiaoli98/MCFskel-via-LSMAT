//
// Created by boy on 22/09/22.
//
#ifndef LSMAT_HPP
#define LSMAT_HPP

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
class LSMAT{
public:
    LSMAT(MyMesh *m){
        tri::UpdateBounding<MyMesh>::Box(*m);
        tri::RequirePerVertexNormal(*m);
        tri::UpdateNormal<MyMesh>::PerVertexNormalized(*m);
        for(int i = 0; i < m->VN(); i++){
//            PRINTP(m->vert[i].P());
            point_list.emplace_back(m->vert[i].P());
            normal_list.emplace_back(m->vert[i].N());
        }
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(0,1.0);

        double r = 0;
        Point3d center;
        double offset = m->bbox.Diag()/10.0f;
        cout <<"offset:" <<offset<<endl;
        for(int i = 0; i < m->VN(); i++){
            r = ((double) rand() / (RAND_MAX)) + 1;
            center = point_list[i].operator-(normal_list[i]*r*offset);
            this->spheres.emplace_back(Sphere3d(center, m->bbox.Diag()*dis(gen)));
        }

        this->spheres_old.resize(m->VN());
#if DEBUG
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

    vector<Point3d> point_list;
    vector<Point3d> normal_list;
    vector<Sphere3d> spheres;
    vector<Sphere3d> spheres_old;
    double sigma = 0;
    double omega_1 = 0.1;     //modifying one will update the other omega_*
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
        return omega_1*e_maximal(s, s_old) + omega_2*e_inscribed(s, s_old);
    }
    double e_maximal(Sphere3d s, Sphere3d s_old){
        double energy = s.Radius() - (s_old.Radius() + epsilon);
        return sqrt(energy);
    }
    //gradient_s of e_inscribed is the derivative_s of Phi_blend()
    double e_inscribed(Sphere3d s, Sphere3d s_old){
        double sum = 0;
        for (int i = 0; i < point_list.size(); i++){
            double _ramp = ramp(s_old.Center().operator-(point_list[i]).Norm() - s_old.Radius());
            double _phi = phi(_ramp, h_support);

            sum += _phi *  pow(phi_blend(s, s_old, point_list[i], normal_list[i]),2);
        }
        return sum;
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
        Point3d temp = p.operator-(s.Center());
        return ramp(s.Radius() - temp.Norm());
    }

    double phi_blend(Sphere3d s, Sphere3d s_old, Point3d p, Point3d n){
        //c^(t-1) - n(c^(t-1)-p) * n
        Point3d projection = s_old.Center().operator-(n.operator*(s_old.Center().operator-(p).operator*(n)));
        return mix(pow(phi_plane(s, p, n),2), pow(phi_point(s, p),2), phi( projection.operator-(p), h_blend));
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
#if DEBUG
        for(int i = 0 ; i < point_list.size(); i++){
            cout<<"point" <<i<<": ";
            PRINTP(point_list[i])
            cout <<"normal:";
            PRINTP(normal_list[i])
        }
#endif
        double energy;
        double energy_old = 0;
        double residual;
        int j = 0;
//        vector<Sphere3d> sphere_list;
//        sphere_list.resize(point_list.size());
        Sphere3d newSphere;
        Point3d newPoint;
        Eigen::Vector4d vec;

        while(j++ < n_runs || residual > threshold){
            for(int i = 0; i<point_list.size(); i++){
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
                spheres[i] = newSphere;
            }
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
#if DEBUG
        cout << "p: ";
        PRINTP(p)
        cout << "n: ";
        PRINTP(n)
        cout << "sphere center: ";
        PRINTP(s.Center());
        cout <<"sphere radius: "<<s.Radius()<<endl;
        cout << "sphere_old center: ";
        PRINTP(s_old.Center());
        cout <<"sphere_old radius: "<<s_old.Radius()<<endl;
        cout << "X_n: "<<endl << X_n<<endl;
#endif

        Eigen::Vector4d gradient_phi_plane;
        if(phi_plane(s, p, n) > 0){
            gradient_phi_plane <<   2*(-n.X())*(s.Radius()-(p.operator-(s.Center()).dot(n))),
                                    2*(-n.Y())*(s.Radius()-(p.operator-(s.Center()).dot(n))),
                                    2*(-n.Z())*(s.Radius()-(p.operator-(s.Center()).dot(n))),
                                    2*(s.Radius()-(p.operator-(s.Center()).dot(n)));
        }
        else{
            gradient_phi_plane << 0, 0, 0, 0;
        }
#if DEBUG
        cout << "phi_plane: " << phi_plane(s, p, n)<<endl;
        cout << "gradient_phi_plane: " <<endl<< gradient_phi_plane<<endl;
#endif
        double distancePC = p.operator-(s.Center()).Norm();
        Eigen::Vector4d gradient_phi_point;
        if(phi_point(s, p) > 0 ){
            gradient_phi_point <<   2*((p.X()-s.Center().X())/distancePC)*(s.Radius() - distancePC),
                                    2*((p.Y()-s.Center().Y())/distancePC)*(s.Radius() - distancePC),
                                    2*((p.Z()-s.Center().Z())/distancePC)*(s.Radius() - distancePC),
                                    2*(s.Radius() - distancePC);
        }
        else{
            gradient_phi_point << 0, 0, 0, 0;
        }
#if DEBUG
        cout << "distance PC: " << distancePC<<endl;
        cout << "phi_point: " << phi_point(s, p)<<endl;
        cout << "gradient_phi_point: " <<endl<< gradient_phi_point<<endl;
#endif
        Eigen::Vector4d gradient_maximal(0,0,0,2*s.Radius()-s_old.Radius()-epsilon);
#if DEBUG
        cout << "gradient_maximal: " <<endl<< gradient_maximal<<endl;
#endif
        double distanceCP = s.Center().operator-(p).Norm();
        Eigen::Vector4d gradient_pinning;
        if(e_pinning(s, p)){
            gradient_pinning << 2*((s.Center().X()-p.X())/distanceCP) * (distanceCP - s.Radius() - d_in),
                                2*((s.Center().Y()-p.Y())/distanceCP) * (distanceCP - s.Radius() - d_in),
                                2*((s.Center().Z()-p.Z())/distanceCP) * (distanceCP - s.Radius() - d_in),
                                -2*(distanceCP - s.Radius() - d_in);
        }
        else{
            gradient_pinning << 0, 0, 0, 0;
        }
#if DEBUG
        cout << "distance CP: " << distanceCP<<endl;
        cout << "e_pinning: "<< e_pinning(s,p)<<endl;
        cout << "gradient_pinning: " <<endl<<gradient_pinning <<endl;
#endif
        Eigen::Vector4d jacobian = gradient_maximal+gradient_phi_plane+gradient_phi_point+gradient_pinning;
        Eigen::Matrix<double, 4, 4> matrix = jacobian * jacobian.transpose();
#if DEBUG
        cout << "jacobian : " << endl<<jacobian<<endl;
        cout << "J^T J: " <<endl<<matrix <<endl;
#endif

        Eigen::Vector4d vec = matrix.inverse() * jacobian;
#if DEBUG
        cout << "(J^T J)^-1 J^T " << endl<<vec<<endl;
#endif
        double energy = e_medial(s, s_old);
        double residual =energy_old - energy;

        vec = vec * residual;
#if DEBUG
        cout <<"X_n - vec = "<< endl<<X_n - vec<<endl;
        cout << "eneegy: "<<energy<<endl;
        cout <<"residual: "<<residual<<endl;
        cin.get();
#endif
        return tuple<Eigen::VectorXd, double, double> {X_n - vec, energy,residual};
    }
};


#endif //LSMAT_HPP