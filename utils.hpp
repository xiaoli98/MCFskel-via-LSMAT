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
#define DEBUG 0

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

#endif //MY_MESH_STRUCTURE

#endif //MCF_SKET_UTILS_HPP
