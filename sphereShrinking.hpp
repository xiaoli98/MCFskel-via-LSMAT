//
// Created by boy on 29/10/22.
//

#ifndef MCF_SKET_SPHERESHRINKING_HPP
#define MCF_SKET_SPHERESHRINKING_HPP


#include <cmath>
#include <iostream>
#include <cstdlib>
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

    KdTree<float> kdtree;

    MyMesh *m;
public:
    const vector<Sphere3d> &getMedialSpheres() const {
        return medial_spheres;
    }

    void setMedialSpheres(const vector<Sphere3d> &medialSpheres) {
        medial_spheres = medialSpheres;
    }

public:
    SphereShrinking(MyMesh *m, KdTree<float> kdtree) : kdtree(kdtree) {
        this->m = m;
        tri::UpdateBounding<MyMesh>::Box(*m);
        tri::RequirePerVertexNormal(*m);
        tri::UpdateNormal<MyMesh>::PerVertexNormalized(*m);
        for(int i = 0; i < m->VN(); i++){
            point_list.emplace_back(m->vert[i].P());
            normal_list.emplace_back(m->vert[i].N());
        }
    }

    double compute_radius(Point3d p1, Point3d n1, Point3d p2){
        double distance = Distance(p1,p2);
        double cos_theta = n1.operator*(p1.operator-(p2))/distance;
        return distance / (2 * cos_theta);
    }

    /// initial point from p
    /// \param p point
    /// \return random point from point_list
    Point3d init_point(Point3d p, Point3d n) {
        Point3d center = p.operator-(n.operator*(m->bbox.Diag()/2));
        return center;
    }

    /// return the nearest neighbor of c in th point_list without p
    /// \param c the point which to find the neighbor
    /// \param p the one outside the set
    /// \return the nearest neighbor point
    Point3d nearestNeighbor(Point3d c, Point3d p){
        auto min_dist = DBL_MAX;
        double dist = 0;
        Point3d min_dist_point;
        for(auto point = point_list.begin(); point < point_list.end(); point++){
            if(point->operator==(p)) continue;
            dist = Distance(c, *point);
            if(dist < min_dist) {
                min_dist = dist;
                min_dist_point = *point;
            }
        }

        return min_dist_point;
    }

    /// calculate the knn using kd-tree
    /// \param c the point which to find the neighbor
    /// \param p the one outside the set
    /// \return the nearest point to c
    Point3d nearestNeighbor_kdtree(Point3d c, Point3d p){
        unsigned int idx;
        float dist;
        KdTree<float>::PriorityQueue queue;
        this->kdtree.doQueryK(c, 2, queue);
        for(int i = 0; i < queue.getNofElements(); i++){
            int idx = queue.getIndex(i);
            if(point_list[idx].operator!=(p))
                return point_list[idx];
        }
    }

    void compute_ma_point(){
        Point3d p;
        Point3d n;
        Point3d p_tilde;
        double radius = 0;
        double radius_new;
        double radius_init;
        Point3d c;

        for(int i = 0; i < point_list.size(); i++){
            cout<<"Point "<<i<<endl;
            p = point_list[i];
            n = normal_list[i];

            //could be improved
            //the original paper says to use the last p_tilde,
            // but if the point cloud is not ordered it brings to negative radius
            p_tilde = init_point(p, n);

            radius_init = compute_radius(p, n,p_tilde);
            radius_new = radius_init;
//#if DEBUG
//            cout << "initial radius:" <<radius_init<<endl;
//#endif
            int z = 0;
            while(true){
                radius = radius_new;
                c = p.operator-(n.operator*(radius));
                p_tilde = nearestNeighbor_kdtree(c, p);
                radius_new = compute_radius(p, n,p_tilde);
#if DEBUG
                if(i == 0)
//                    add_octahedron(*m, c, radius_new, "SS_"+ to_string(z++)+".off");
                    add_sphere(*m, c, radius_new);
#endif
                if(radius - radius_new < 0.0001)
                    break;
            };
//#if DEBUG
//            MyMesh m1;
//            add_sphere(m1, c, radius_new, Color4b::Red);
//            tri::Sphere(m1);
//            tri::Append<MyMesh, MyMesh>::Mesh(m1, *m);
//            string filename = "sphere_of_" + to_string(i) + ".off";
//            tri::io::ExporterOFF<MyMesh>::Save(m1, filename.c_str(), tri::io::Mask::IOM_FACECOLOR);
//#endif
            medial_spheres.emplace_back(Sphere3d(c, radius));
        }
    }

};

#endif //MCF_SKET_SPHERESHRINKING_HPP
