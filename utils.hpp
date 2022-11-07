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

#define LS_MAT 1
#define SPHERE_SHRINK 0

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
            {0, 1, 2},
            {0, 2, 3},
            {0, 3, 4},
            {0, 4, 1},
            {1, 5, 2},
            {2, 5, 3},
            {3, 5, 4},
            {4, 5, 1}
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
    ScalarType L = ScalarType((math::Sqrt(5.0) + 1.0) / 2.0);
    CoordType vv[12] = {
            CoordType(0, L, 1),
            CoordType(0, L, -1),
            CoordType(0, -L, 1),
            CoordType(0, -L, -1),

            CoordType(L, 1, 0),
            CoordType(L, -1, 0),
            CoordType(-L, 1, 0),
            CoordType(-L, -1, 0),

            CoordType(1, 0, L),
            CoordType(-1, 0, L),
            CoordType(1, 0, -L),
            CoordType(-1, 0, -L)
    };

    int ff[20][3] = {
            {1,  0,  4},
            {0,  1,  6},
            {2,  3,  5},
            {3,  2,  7},
            {4,  5,  10},
            {5,  4,  8},
            {6,  7,  9},
            {7,  6,  11},
            {8,  9,  2},
            {9,  8,  0},
            {10, 11, 1},
            {11, 10, 3},
            {0,  8,  4},
            {0,  6,  9},
            {1,  4,  10},
            {1,  11, 6},
            {2,  5,  8},
            {2,  9,  7},
            {3,  10, 5},
            {3,  7,  11}
    };

    newMesh.Clear();
    tri::Allocator<MeshType>::AddVertices(newMesh, 12);
    tri::Allocator<MeshType>::AddFaces(newMesh, 20);
    VertexPointer ivp[12];

    VertexIterator vi;
    int i;
    for (i = 0, vi = newMesh.vert.begin(); vi != newMesh.vert.end(); ++i, ++vi) {
        (*vi).P() = vv[i];
        ivp[i] = &*vi;
    }

    FaceIterator fi;
    for (i = 0, fi = newMesh.face.begin(); fi != newMesh.face.end(); ++i, ++fi) {
        (*fi).V(0) = ivp[ff[i][0]];
        (*fi).V(1) = ivp[ff[i][1]];
        (*fi).V(2) = ivp[ff[i][2]];
        (*fi).C().ColorRamp(0, 20, i);
    }

    tri::Sphere(newMesh, 1);
    tri::UpdatePosition<MyMesh>::Translate(newMesh, center);
    double scale_factor = radius / (2 * sin(2*M_PI/5));

    tri::UpdatePosition<MyMesh>::Scale(newMesh, scale_factor);
    tri::Append<MeshType,MeshType>::Mesh(m, newMesh);
}
#endif //MY_MESH_STRUCTURE

#endif //MCF_SKET_UTILS_HPP
