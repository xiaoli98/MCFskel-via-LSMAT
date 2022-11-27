//
// Created by boy on 29/10/22.
//

#ifndef MCF_SKET_UTILS_HPP
#define MCF_SKET_UTILS_HPP

#include <vcg/complex/complex.h>
#include <vcg/simplex/vertex/component.h>
#include <Eigen/Dense>
#include<wrap/io_trimesh/import_obj.h>
#include<wrap/io_trimesh/export_obj.h>
#include<wrap/io_trimesh/import_off.h>
#include<wrap/io_trimesh/export_off.h>

#define PRINTP(p) cout << p.X() <<"\t"<< p.Y() <<"\t"<< p.Z() <<endl;
#define PRINTPP(p) cout << p->X() <<"\t"<< p->Y() <<"\t"<< p->Z() <<endl;
#define DEBUG 1
#define DEBUG_FUNCTION 0

#define LS_MAT 0
#define SPHERE_SHRINK 1

#ifndef MY_MESH_STRUCTURE
#define MY_MESH_STRUCTURE

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

void add_octahedron(MyMesh &m, Point3d p, double radius, string filename){
    MyMesh newMesh;
//    radius = radius/2;
    PRINTP(p)
    cout <<"radius: "<<radius<<endl;
    Point3d vv[6] = {
            Point3d(p.X(), p.Y()+radius, p.Z()),
            Point3d(p.X()+radius, p.Y(), p.Z()),
            Point3d(p.X(), p.Y(), p.Z()+radius),
            Point3d(p.X()-radius, p.Y(), p.Z()),
            Point3d(p.X(), p.Y(), p.Z()-radius),
            Point3d(p.X(), p.Y()-radius, p.Z()),
    };
    int ff[8][3] = {
            {0, 2,1},
            {0, 3, 2},
            {0, 4, 3},
            {0, 1, 4},
            {1, 2, 5},
            {2, 3, 5},
            {3, 4, 5},
            {4, 1, 5}
    };
    tri::Allocator<MyMesh>::AddVertices(newMesh, 6);
    tri::Allocator<MyMesh>::AddFaces(newMesh, 8);
    MyMesh::VertexPointer ivp[6];

    MyMesh::VertexIterator vi;
    int i;
    for (i = 0, vi = newMesh.vert.begin(); vi != newMesh.vert.end(); ++i, ++vi) {
        (*vi).P() = vv[i];
        ivp[i] = &*vi;
    }

    MyMesh::FaceIterator fi;
    for (i = 0, fi = newMesh.face.begin(); fi != newMesh.face.end(); ++i, ++fi) {
        (*fi).V(0) = ivp[ff[i][0]];
        (*fi).V(1) = ivp[ff[i][1]];
        (*fi).V(2) = ivp[ff[i][2]];
    }
//    tri::Sphere(newMesh, 3);
//    MyMesh mesh1;
//    tri::Sphere(mesh1);
//    tri::RequirePerVertexNormal(mesh1);
//    tri::UpdateNormal<MyMesh>::PerVertexNormalized(mesh1);
//    Point3d mesh1_centre = mesh1.vert[0].P().operator-(mesh1.vert[0].N()*radius);
//    Point3d direction = p.operator-(mesh1_centre).normalized();
//    double dist = Distance(p, mesh1_centre);
//    tri::UpdatePosition<MyMesh>::Translate(mesh1, dist*direction);
//    tri::UpdatePosition<MyMesh>::Scale(mesh1, radius);
//    tri::Append<MyMesh, MyMesh>::Mesh(m, mesh1);

//    tri::Sphere(newMesh,1);

    for(int i = 0 ; i < 3; ++i)
    {
        MyMesh newM;
        for(MyMesh::FaceIterator fi=newMesh.face.begin();fi!=newMesh.face.end();++fi)
        {
            MyMesh::CoordType me01 =  (fi->P(0)+fi->P(1))/2.0;
            MyMesh::CoordType me12 =  (fi->P(1)+fi->P(2))/2.0;
            MyMesh::CoordType me20 =  (fi->P(2)+fi->P(0))/2.0;
            tri::Allocator<MyMesh>::AddFace(newM,me01,me12,me20);
            tri::Allocator<MyMesh>::AddFace(newM,fi->P(0),me01,me20);
            tri::Allocator<MyMesh>::AddFace(newM,fi->P(1),me12,me01);
            tri::Allocator<MyMesh>::AddFace(newM,fi->P(2),me20,me12);
        }
        tri::Clean<MyMesh>::RemoveDuplicateVertex(newM);
        tri::Append<MyMesh,MyMesh>::MeshCopy(newMesh,newM);

//        for(MyMesh::VertexIterator vi = newMesh.vert.begin(); vi != newMesh.vert.end(); ++vi)
//            vi->P().Normalize();
    }

    tri::Append<MyMesh, MyMesh>::Mesh(m, newMesh);
    tri::io::ExporterOFF<MyMesh>::Save(m, filename.c_str(), tri::io::Mask::IOM_FACECOLOR);
}

template<class MeshType>
void add_sphere(MeshType &m, Point3d center, double radius) {
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::VertexPointer VertexPointer;
    typedef typename MeshType::VertexIterator VertexIterator;
    typedef typename MeshType::FaceIterator FaceIterator;

    MeshType newMesh;

    tri::Sphere(newMesh,2);
    double scale_factor = radius;
    tri::UpdatePosition<MyMesh>::Scale(newMesh, scale_factor);
    tri::UpdatePosition<MyMesh>::Translate(newMesh, center);
    tri::Append<MeshType,MeshType>::Mesh(m, newMesh);
//    tri::io::ExporterOFF<MyMesh>::Save(m, filename.c_str(), tri::io::Mask::IOM_FACECOLOR);
}
#endif //MY_MESH_STRUCTURE

#endif //MCF_SKET_UTILS_HPP
