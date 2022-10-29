//
// Created by boy on 29/10/22.
//

#ifndef MCF_SKET_SPHERESHRINKING_HPP
#define MCF_SKET_SPHERESHRINKING_HPP


#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <random>
#include <string>

#include <vcg/complex/algorithms/create/platonic.h>

#include "utils.hpp"

using namespace vcg;
using namespace std;

class SphereShrinking {
private:
    vector<Point3d> point_list;
    vector<Point3d> normal_list;

    vector<Sphere3d> medial_spheres;

    MyMesh *m;
public:
    const vector<Sphere3d> &getMedialSpheres() const {
        return medial_spheres;
    }

    void setMedialSpheres(const vector<Sphere3d> &medialSpheres) {
        medial_spheres = medialSpheres;
    }

public:
    SphereShrinking(MyMesh* m){
        this->m = m;
        tri::UpdateBounding<MyMesh>::Box(*m);
        tri::RequirePerVertexNormal(*m);
        tri::UpdateNormal<MyMesh>::PerVertexNormalized(*m);
        for(int i = 0; i < m->VN(); i++){
            point_list.emplace_back(m->vert[i].P());
            normal_list.emplace_back(m->vert[i].N());
        }
    }
    double safe_acos(double value) {
        if (value<=-1.0) {
            return M_PI;
        } else if (value>=1.0) {
            return 0;
        } else {
            return acos(value);
        }
    }

    double compute_radius(Point3d p1, Point3d n1, Point3d p2){
        double distance = Distance(p1,p2);
#if DEBUG
        cout << "\t##################################RADIUS###########################################"<<endl;
        cout <<"\tp1: ";
        PRINTP(p1)
        cout <<"\tn1: ";
        PRINTP(n1)
        cout<<"\tp2: ";
        PRINTP(p2)
        cout << "\tdistance: "<< distance<<endl;
        cout << "\tnumerator of theta: "<< n1.operator*(p1.operator-(p2)) << endl;
#endif
        double theta = safe_acos((n1.operator*(p1.operator-(p2))) / distance);
#if DEBUG
        cout << "\ttheta: "<< theta<<endl;
        cout << "\tout: "<< distance / (2 * cos(theta))<<endl;
        cout << "\t##################################RADIUS###########################################"<<endl;
#endif
        return distance / (2 * cos(theta));
    }

    /// random point different from p
    /// \param p point
    /// \return random point from point_list
    Point3d randomPoint(Point3d p) {
        size_t i =0;
        do {
            i = rand() % (point_list.size() - 1);
#if DEBUG
            cout <<"i = "<<i<<endl;
            cout <<"p[i]: ";
            PRINTP(point_list[i])
#endif
        }while(point_list[i] == p);
        return point_list[i];
    }

    /// return the nearest neighbor of p in th point_list
    /// \param p the point which to find the neighbor
    /// \return the nearest neighbor point
    Point3d nearestNeighbors(Point3d p){
        auto min_dist = DBL_MAX;
        double dist = 0;
        Point3d min_dist_point;
        for(auto point = point_list.begin(); point < point_list.end() && *point != p; point++){
            dist = Distance(p, *point);
            if(dist < min_dist) {
                min_dist = dist;
                min_dist_point = *point;
            }
        }
#if DEBUG
        cout <<"NN of: ";
        cout << p.X() <<"\t"<< p.Y() <<"\t"<< p.Z();
        cout <<"\tis: ";
        PRINTP(min_dist_point)
#endif
        return min_dist_point;
    }

    void compute_ma_point(){
        Point3d p;
        Point3d n;
        Point3d p_tilde;
        double radius = 0;
        double radius_new;
        Point3d c;
        for(int i = 0; i < point_list.size(); i++){
            p = point_list[i];
            n = normal_list[i];
#if DEBUG
            cout << "p:";
            PRINTP(p)
            cout << "n:";
            PRINTP(n)
#endif
            if(p == *(point_list.begin()))
                p_tilde = randomPoint(p);
#if DEBUG
            cout <<"random point:";
            PRINTP(p_tilde);
#endif
            radius = compute_radius(p, n, p_tilde);
            radius_new = radius;

            string fileName = "Sphere_" + to_string(i) +"_init.off";
            MyMesh myMeshInit;
            create_sphere(myMeshInit, p_tilde, radius);
            tri::io::ExporterOFF<MyMesh>::Save(myMeshInit, fileName.c_str(), tri::io::Mask::IOM_FACECOLOR);
            int j = 0;
            do{
                radius = radius_new;
                c = p.operator-(n.operator*(radius));
                p_tilde = nearestNeighbors(c);
                radius_new = compute_radius(p, n, p_tilde);
#if DEBUG
                cout << "radius: "<<radius <<"\t"<<"radius_new:" << radius_new<<endl;

                string fileName = "Sphere_" + to_string(i) +"_"+ to_string(j)+".off";

                MyMesh myMesh;
                create_sphere(myMesh, c, radius_new);
                tri::io::ExporterOFF<MyMesh>::Save(myMesh, fileName.c_str(), tri::io::Mask::IOM_FACECOLOR);
                j++;
                cin.get();
#endif
            }while(radius != radius_new);
            medial_spheres.emplace_back(Sphere3d(c, radius));
        }
    }

    template<class MeshType>
    void create_sphere(MeshType &m, Point3d center, double radius){
        typedef typename MeshType::ScalarType ScalarType;
        typedef typename MeshType::CoordType CoordType;
        typedef typename MeshType::VertexPointer  VertexPointer;
        typedef typename MeshType::VertexIterator VertexIterator;
        typedef typename MeshType::FaceIterator   FaceIterator;
//        std::mt19937 generator ();
//        std::uniform_real_distribution<double> uniform01(0.0, 1.0);
//
//        double theta, phi, x, y, z;
        int n_vert = 12;
        int n_faces = 20;
        tri::Allocator<MyMesh>::AddVertices(m, n_vert);
        tri::Allocator<MyMesh>::AddFaces(m, n_faces);

        ScalarType L = ScalarType((math::Sqrt(5.0)+1.0)/2.0);
        CoordType vv[12]={
                CoordType ( 0, L, 1),
                CoordType ( 0, L,-1),
                CoordType ( 0,-L, 1),
                CoordType ( 0,-L,-1),

                CoordType ( L, 1, 0),
                CoordType ( L,-1, 0),
                CoordType (-L, 1, 0),
                CoordType (-L,-1, 0),

                CoordType ( 1, 0, L),
                CoordType (-1, 0, L),
                CoordType ( 1, 0,-L),
                CoordType (-1, 0,-L)
        };

        int ff[20][3]={
                {1,0,4},{0,1,6},{2,3,5},{3,2,7},
                {4,5,10},{5,4,8},{6,7,9},{7,6,11},
                {8,9,2},{9,8,0},{10,11,1},{11,10,3},
                {0,8,4},{0,6,9},{1,4,10},{1,11,6},
                {2,5,8},{2,9,7},{3,10,5},{3,7,11}
        };

        m.Clear();
        tri::Allocator<MeshType>::AddVertices(m,12);
        tri::Allocator<MeshType>::AddFaces(m,20);
        VertexPointer ivp[12];

        VertexIterator vi;
        int i;
        for(i=0,vi=m.vert.begin();vi!=m.vert.end();++i,++vi){
            (*vi).P() = radius * center.operator+(vv[i]);
            ivp[i]=&*vi;
        }

        FaceIterator fi;
        for(i=0,fi=m.face.begin();fi!=m.face.end();++i,++fi){
            (*fi).V(0)=ivp[ff[i][0]];
            (*fi).V(1)=ivp[ff[i][1]];
            (*fi).V(2)=ivp[ff[i][2]];
        }

        tri::Sphere(m, 4);
//        tri::Append<MyMesh, MyMesh>::MeshCopy(m, in);
//
//        CoordType coord;
//        VertexIterator vi;
//        int i;
//        for (i = 0, vi=in.vert.begin(); vi != in.vert.end(); i++, vi++){
//
//            theta = 2 * M_PI * uniform01(generator);
//            phi = M_PI * uniform01(generator);
//            x = center.X() + radius * sin(phi) * cos(theta);
//            y = center.Y() + radius * sin(phi) * sin(theta);
//            z = center.Z() + radius * cos(phi);
//            coord = CoordType(x,y,z);
//            (*vi).P() = coord;
//        }

    }
};

#endif //MCF_SKET_SPHERESHRINKING_HPP
