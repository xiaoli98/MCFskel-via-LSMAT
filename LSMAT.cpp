//
// Created by boy on 22/09/22.
//
#include <cmath>

#include <vcg/complex/complex.h>
#include <vcg/simplex/vertex/component.h>
#include <Eigen/Dense>
#include<wrap/io_trimesh/import_obj.h>
#include<wrap/io_trimesh/export_obj.h>
#include<wrap/io_trimesh/import_off.h>
#include<wrap/io_trimesh/export_off.h>

#define PRINTP(p) cout << p.X() <<" "<< p.Y() <<" "<< p.Z() <<endl;
#define PRINTPP(p) cout << p->X() <<"\t"<< p->Y() <<"\t"<< p->Z() <<endl;

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
public:
    LSMAT(MyMesh *m){
        tri::UpdateBounding<MyMesh>::Box(*m);
        tri::RequirePerVertexNormal(*m);
        tri::UpdateNormal<MyMesh>::PerVertexNormalized(*m);
        for(int i = 0; i < m->VN(); i++){
            point_list.emplace_back((Point3d*)&(m->vert[i].P()));
            normal_list.emplace_back((Point3d*)&(m->vert[i].N()));
        }
        this->spheres = vector<Sphere3d>(m->VN(), Sphere3d(m->bbox.Center(),m->bbox.Diag()));
        this->spheres_old.resize(m->VN());
    }


private:

    vector<Point3d*> point_list;
    vector<Point3d*> normal_list;
    vector<Sphere3d> spheres;
    vector<Sphere3d> spheres_old;
    double sigma = 0;
    double omega_1 = 1;     //modifying one will update the other omega_*
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
    //gradient of e_inscribed is the derivative of Phi_blend()
    double e_inscribed(Sphere3d s, Sphere3d s_old){
        double sum = 0;
        for (int i = 0; i < point_list.size(); i++){
            double _ramp = ramp(s_old.Center().operator-(*point_list[i]).Norm() - s_old.Radius());
            double _phi = phi(_ramp, h_support);

            sum += _phi*  pow(phi_blend(s, s_old, *point_list[i], *normal_list[i], h_blend),2);
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

        double sum_square = 0;
        for (int i= 0; i < 3; i++){
            sum_square += pow(temp.V(i), 2);
        }
        return ramp(s.Radius() - sqrt(sum_square));
    }

    double phi_blend(Sphere3d s, Sphere3d s_old, Point3d p, Point3d n, double h_blend){
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
        for(int i = 0 ; i < point_list.size(); i++){
            cout<<"point" <<i<<": ";
            PRINTPP(point_list[i])
        }
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
                tuple<Eigen::VectorXd, double, double> t = gauss_newton(*point_list[i], *normal_list[i], spheres[i], spheres_old[i], energy_old);
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

        Eigen::Vector4d gradient_phi_plane;
        if(phi_plane(s, p, n) > 0){
            gradient_phi_plane << -n.X(), -n.Y(), -n.Z(), 1;
        }
        else{
            gradient_phi_plane << 0, 0, 0, 0;
        }

        double distancePC = p.operator-(s.Center()).Norm();
        Eigen::Vector4d gradient_phi_point;
        if(phi_point(s, p) > 0 ){
            gradient_phi_point <<   (p.X()-s.Center().X())/distancePC,
                                    (p.Y()-s.Center().Y())/distancePC,
                                    (p.Z()-s.Center().Z())/distancePC,
                                    1;
        }
        else{
            gradient_phi_point << 0, 0, 0, 0;
        }
        Eigen::Vector4d gradient_maximal(0,0,0,1);
        double distanceCP = s.Center().operator-(p).Norm();
        Eigen::Vector4d gradient_pinning;
        if(e_pinning(s, p)){
            gradient_pinning << (s.Center().X()-p.X())/distanceCP,
                                (s.Center().Y()-p.Y())/distanceCP,
                                (s.Center().Z()-p.Z())/distanceCP,
                                -1;
        }
        else{
            gradient_pinning << 0, 0, 0, 0;
        }
        Eigen::Vector4d jacobian(2*gradient_maximal+2*gradient_phi_plane+2*gradient_phi_point+gradient_pinning);
        Eigen::Matrix<double, 4, 4> matrix = jacobian * jacobian.transpose();

        Eigen::Vector4d vec = matrix.inverse() * jacobian;
        double energy = e_medial(s, s_old);
        double residual = abs(energy_old - energy);

        vec = vec * residual;
        cout <<"X_n - vec = "<< X_n - vec<<endl;
        cout << "eneegy: "<<energy<<endl;
        cout <<"residual: "<<residual<<endl;
        return tuple<Eigen::VectorXd, double, double> {X_n - vec, energy,residual};
    }
};

int main(int argc, char* argv[]){
    MyMesh m;

    tri::io::ImporterOFF<MyMesh>::Open(m,argv[1]);
    tri::UpdateBounding<MyMesh>::Box(m);
    tri::RequirePerVertexNormal(m);
    tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);

//    vector<Point3d> points;
//    vector<Point3d> normals;
//    for(int i = 0; i < m.VN(); i++){
//        points.emplace_back(m.vert[i].P());
//        normals.emplace_back(m.vert[i].N());
//    }


    Point3d p;

    LSMAT lsmat(&m);
    for(int j = 0; j < 1; j++){
        cout << "*****************LSMAT iteration "<< j <<"******************"<<endl;
        string filename = "LSMAT_"+ to_string(j)+".off";
        vector<Sphere3d> spheres = lsmat.compute(1);

        for(int i = 0; i<m.VN(); i++){
            m.vert[i].P() = spheres[i].Center();
        }

        tri::io::ExporterOFF<MyMesh>::Save(m, filename.c_str(), tri::io::Mask::IOM_FACECOLOR);
    }

}