//
// Created by boy on 04/12/22.
//

#ifndef MCFSKET_COLLAPSER_H
#define MCFSKET_COLLAPSER_H

#define COLLAPESER_DEBUG 1

#include <map>
#include "utils.hpp"

class Collapser{
    typedef MyMesh::PerVertexAttributeHandle<double> omega;
    typedef vector<tuple<MyMesh::VertexType*, MyMesh::VertexType*, double>> laplacian_triple;
private:
    MyMesh *m;
    double omega_L_0;
    omega omega_L;
    double omega_H_0;
    omega omega_H;
    double omega_M_0;
    omega omega_M;
    double zero_TH;

    Eigen::SparseMatrix<double> LHS;
    Eigen::MatrixXd RHS;
    Eigen::MatrixXd X;

    int nrows, ncols;

    MyMesh::PerVertexAttributeHandle<MyMesh::CoordType> vert_mat;
    MyMesh::PerVertexAttributeHandle<bool> isFixed;
    MyMesh::PerVertexAttributeHandle<bool> isSplitted;
    MyMesh::PerMeshAttributeHandle<laplacian_triple> laplacian;
    MyMesh::PerMeshAttributeHandle<double> edge_threshold;

    //some utility functions
//    bool isInCommon(const MyMesh::VertexType *v, const MyMesh::FaceType *fb){
//        for(int i = 0; i < 3; i++){
//            if(fb->V(i) == v) return true;
//        }
//        return false;
//    }
//    int edge_index(const MyMesh::FaceType* f, const MyMesh::VertexType *v){
//        for(int i = 0; i < 3; i++){
//            if(f->V(i) == v) return i;
//        }
//        return -1;
//    }

//    int third_vert(const MyMesh::FaceType *f, const MyMesh::VertexType *v1, const MyMesh::VertexType *v2){
//        for(int i = 0; i < 3; i++){
//            if(f->V(i) != v1 && f->V(i) != v2) return i;
//        }
//    }
public:
    Collapser(MyMesh* m){
        this->m = m;
        omega_L = tri::Allocator<MyMesh>::GetPerVertexAttribute<double>(*m, string("omega_L"));
        omega_L_0 = 1;
        omega_H = tri::Allocator<MyMesh>::GetPerVertexAttribute<double>(*m, string("omega_H"));
        omega_H_0 = 20;
        omega_M = tri::Allocator<MyMesh>::GetPerVertexAttribute<double>(*m, string("omega_M"));
        omega_M_0 = 40;
        for (auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            omega_L[vi] = omega_L_0;
            omega_H[vi] = omega_H_0;
            omega_M[vi] = omega_M_0;
        }
        zero_TH = 1e-7;

        vert_mat = tri::Allocator<MyMesh>::GetPerVertexAttribute<MyMesh::CoordType>(*m, string("vert_mat"));
        isFixed = tri::Allocator<MyMesh>::GetPerVertexAttribute<bool>(*m, string("isFixed"));
        isSplitted = tri::Allocator<MyMesh>::GetPerVertexAttribute<bool>(*m, string("isSplitted"));
        edge_threshold = tri::Allocator<MyMesh>::GetPerMeshAttribute<double>(*m, string("edge_threshold"));
        edge_threshold() = 0.002*m->bbox.Diag();
    }

    void compute(){
        createLHS();
        createRHS();
        solveLS();
        MyMesh::PerVertexAttributeHandle<MyMesh::CoordType> sol = tri::Allocator<MyMesh>::GetPerVertexAttribute<MyMesh::CoordType>(*m, string("solution"));

        for (auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
            vi->P() = sol[vi];
        }
        update_omega();
        update_topology();
        detect_degeneracies();
        tri::Allocator<MyMesh>::CompactFaceVector(*m);
        tri::Allocator<MyMesh>::CompactVertexVector(*m);

        tri::io::ExporterOFF<MyMesh>::Save(*m, "prova.off", tri::io::Mask::IOM_FACECOLOR);

        int counter = 0;
        for(auto fit = m->face.begin(); fit != m->face.end(); fit++){
            printf("-------------FACE %lu----------------\n", tri::Index(*m, fit.base()));
            if(fit->IsD()) continue;
            for(int i=0; i<fit->VN(); i++){
                PRINTP(fit->V(i)->P());
//                f == f->FFp(e)->FFp(f->FFi(e))
                MyMesh::FacePointer temp = fit->FFp(i)->FFp(fit->FFi(i));

                if(fit.base() != temp) counter++;
            }
        }
        cout << "non manifold face: "<<counter<<endl;

    }

    void createLHS(){
        laplacian = tri::Allocator<MyMesh>::GetPerMeshAttribute<laplacian_triple>(*m, string("laplacian"));
        nrows = 3 * m->VN();
        ncols = m->VN();
        LHS.setZero();
        LHS.resize(nrows, ncols);
        RHS = Eigen::MatrixXd::Zero(nrows, 3);
        X = Eigen::MatrixXd::Zero(ncols, 3);

        vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(ncols*9);

        for (auto l = laplacian().begin(); l != laplacian().end(); l++){
            triplets.emplace_back(Eigen::Triplet<double>(tri::Index(*m, get<0>(*l)), tri::Index(*m, get<1>(*l)), get<2>(*l)));
        }

        for(auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
            triplets.emplace_back(Eigen::Triplet<double>(tri::Index(*m, *vi)+ncols, tri::Index(*m, *vi), omega_H[*vi]));
        }
        for(auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
            triplets.emplace_back(Eigen::Triplet<double>(tri::Index(*m, *vi)+2*ncols, tri::Index(*m, *vi), omega_M[*vi]));
        }
#if COLLAPESER_DEBUG
        cout << triplets.size()<<endl;
        cout << "LHS rows:" << LHS.rows()<<endl;
        cout << "LHS cols:" << LHS.cols()<<endl;
        for (auto it = triplets.begin(); it != triplets.end(); it++){
            if (it->row()<0 || it->row()>=LHS.rows() || it->col()<0 || it->col()>=LHS.cols())
                printf("col: %d\trow:%d\t value: %f\n", it->col(), it->row(), it->value());
        }

#endif
        LHS.setFromTriplets(triplets.begin(), triplets.end());
    }

    void createRHS(){
        for(auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
            MyMesh::CoordType c = (float)omega_H[vi] * vi->P();
            RHS.row(ncols + tri::Index(*m, *vi)) = Eigen::Vector3d(c.X(), c.Y(), c.Z());
        }

        for(auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
            MyMesh::CoordType c = (float) omega_M[vi] * vert_mat[vi];
            RHS.row(2*ncols + tri::Index(*m, *vi)) = Eigen::Vector3d(c.X(), c.Y(), c.Z());
        }
    }

    void solveLS(){
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(LHS.transpose() * LHS);

        X.col(0) = solver.solve(LHS.transpose() * RHS.col(0));
        X.col(1) = solver.solve(LHS.transpose() * RHS.col(1));
        X.col(2) = solver.solve(LHS.transpose() * RHS.col(2));

        if (!std::isfinite(X.norm())){
            throw EIGEN_EXCEPTIONS("Problem with linear least square solution.");
        }

        MyMesh::PerVertexAttributeHandle<MyMesh::CoordType> sol = tri::Allocator<MyMesh>::GetPerVertexAttribute<MyMesh::CoordType>(*m, string("solution"));
        for(auto vi=m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
            Eigen::Vector3d p = X.row(tri::Index(*m, *vi));
            sol[vi] = MyMesh::CoordType(p[0], p[1], p[2]);
        }
    }

    void update_omega(){
        for (auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if(vi->IsD()) continue;
            if(isFixed[vi]){
                omega_L[vi] = 0;
                omega_H[vi] = 1.0/zero_TH;
                omega_M[vi] = 0;
                continue;
            }

            omega_L[vi] = omega_L_0;
            omega_H[vi] = omega_H_0;
            omega_M[vi] = omega_M_0;

            if(isSplitted[vi]){
                omega_L[vi] = omega_L_0;
                omega_H[vi] = omega_H_0;
                omega_M[vi] = 0;
            }
        }
    }

    typedef tuple<MyMesh::VertexType*,MyMesh::VertexType*,MyMesh::VertexType*> triVert;
    void update_topology(){
        //todo
        // Pos is not suitable for mesoskeletons
        // this is not manifold, so the navigation is not good
        // a possible idea is to use only VF adjacency

        edge_collapse();
        edge_split();
        cout << "topology updated, face and vertex cleaned"<<endl;
    }

    void edge_collapse(){
        int count=0;
        //edge collap
        MyMesh::VertexType *V0, *V1;
        vector<triVert> newFaces;
        double edge_length;
//        tri::UpdateFlags<MyMesh>::FaceClearV(*m);
        tri::UpdateFlags<MyMesh>::VertexClearV(*m);

        vector<MyMesh::VertexPointer> adjVertices;
        newFaces.reserve(4);
        for(auto vit = m->vert.begin(); vit != m->vert.end(); vit++) {
            if(!vit->IsD()) {
                face::VVStarVF<MyMesh::FaceType>(vit.base(), adjVertices);
                for (auto adjv = adjVertices.begin(); adjv != adjVertices.end(); adjv++) {
                    if ((*adjv)->IsD()) continue;
                    MyMesh::FacePointer fp = vit->VFp();
                    V0 = vit.base();
                    V1 = (*adjv);

                    edge_length = Distance(V0->P(), V1->P());
                    if (edge_length < edge_threshold()) {
#if COLLAPESER_DEBUG
                        cout << "edge_length:" << edge_length << " threshold: " << edge_threshold() << endl;
#endif
                        V1->P() = (V1->P() + V0->P()) / 2;

                        auto mat0 = vert_mat[V0];
                        auto mat1 = vert_mat[V1];
                        //find the closest medial axis to V1
                        double dist1 = Distance(mat0, V1->P());
                        double dist2 = Distance(mat1, V1->P());

                        if (dist1 < dist2)
                            vert_mat[V1] = vert_mat[V0];
                        else
                            vert_mat[V1] = vert_mat[V1];

                        vcg::face::Pos<MyMesh::FaceType> p(fp, V0);
                        vector<triVert> temp;

                        MyMesh::VertexType *anchorVert = p.VFlip();
                        MyMesh::FaceType *start_face = p.F();
                        vector<MyMesh::VertexType *> to_be_connected;
                        do {
                            //todo
                            // problem with manifoldness
                            if (!p.F()->IsD()) {
                                tri::Allocator<MyMesh>::DeleteFace(*m, *p.F());
                            }
                            p.FlipE();
                            if (p.VFlip() != anchorVert && !p.VFlip()->IsD()) {
                                to_be_connected.emplace_back(p.VFlip());
                            }
                            p.FlipF();
                        } while (p.F() != start_face);

                        temp.reserve(to_be_connected.size() - 1);
                        for (int j = 0; j < to_be_connected.size() - 1; j++) {
                            temp.emplace_back(triVert(anchorVert, to_be_connected[j], to_be_connected[j + 1]));
                        }
                        printf("deleting vertex %p\n", V0);
                        tri::Allocator<MyMesh>::DeleteVertex(*m, *V0);

                        newFaces.insert(newFaces.end(), temp.begin(), temp.end());
                        count++;
                        break;// one ring faces on vit are deleted, so if there was any other edges, are collapsed together
                    }
                }
                MyMesh::FaceIterator fi = tri::Allocator<MyMesh>::AddFaces(*m, newFaces.size());
                for(auto it = newFaces.begin();fi!=m->face.end(); fi++, it++){
                    fi->V0(0) = get<0>(*it);
                    fi->V1(0) = get<1>(*it);
                    fi->V2(0) = get<2>(*it);
                }

                newFaces.clear();
                newFaces.reserve(4);
            }
        }

        tri::UpdateTopology<MyMesh>::VertexFace(*m);
//        tri::UpdateTopology<MyMesh>::FaceFace(*m);
        printf("collapsed %d edges, size: %lu\n", count, newFaces.size());
        tri::Allocator<MyMesh>::CompactFaceVector(*m);
        tri::Allocator<MyMesh>::CompactVertexVector(*m);
    }

    void edge_split(){
        tri::UpdateFlags<MyMesh>::VertexClearV(*m);
        double angle_threshold = 110*(3.14/180); //110° in radiant
        vector<triVert> newFaces;
        vector<triVert> splittingFaceVerts;
        vector<tuple<MyMesh::CoordType, double>> newCoords;

        for(auto vit = m->vert.begin(); vit != m->vert.end(); vit++){
            vit->SetV();
            face::VFIterator<MyMesh::FaceType> vfi(vit.base());
            for(;!vfi.End(); ++vfi){
                MyMesh::FacePointer faceA = vfi.F();
                if(faceA->IsD()) continue;
                int vn = vfi.F()->VN();

                MyMesh::VertexPointer from_vert = vit.base();
                MyMesh::VertexPointer to_vert = vfi.F()->V((vfi.I()+1)%vn);
                MyMesh::VertexPointer edge_oppositeVert = vfi.F()->V((vfi.I()+2)%vn);
                if(to_vert->IsV()) continue; //edge already processed
                if(Distance(from_vert->P(), to_vert->P()) < zero_TH) continue;//too short edge, should be collapsed

                double opposite_angle = Angle(from_vert->cP() - edge_oppositeVert->cP(), to_vert->cP() - edge_oppositeVert->cP());
                if(opposite_angle < angle_threshold) continue;
                MyMesh::CoordType projector = (to_vert->P()-from_vert->P()).normalized();
                MyMesh::CoordType projectee = vit->P() - from_vert->P();
                double t = projector.dot(projectee);

                newCoords.emplace_back(tuple<MyMesh::CoordType, double>(from_vert->P() + projector * t, t));
                tri::Allocator<MyMesh>::DeleteFace(*m, *faceA);
                splittingFaceVerts.emplace_back(triVert(from_vert, edge_oppositeVert, to_vert));
            }
        }

        //allocating new vertices and new faces
        tri::Allocator<MyMesh>::PointerUpdater<MyMesh::VertexPointer> vpu;
        printf("before vn: %d\t vert.size: %lu\n", m->vn, m->vert.size());
        MyMesh::VertexIterator vi = tri::Allocator<MyMesh>::AddVertices(*m, newCoords.size(), vpu);
        printf("after vn : %d\t vert.size: %lu\n", m->vn, m->vert.size());

        assert(splittingFaceVerts.size() == newCoords.size());
//        if(vpu.NeedUpdate()){
//            for(auto it = splittingFaceVerts.begin(); it != splittingFaceVerts.end(); it++){
//                vpu.Update(get<0>(*it));
//                vpu.Update(get<1>(*it));
//                vpu.Update(get<2>(*it));
//            }
//        }

        auto splittingFaceIt = splittingFaceVerts.begin();
        auto coordsIt = newCoords.begin();
        while(vi != m->vert.end()){
            MyMesh::CoordType coord = get<0>(*coordsIt);
            double t = get<1>(*coordsIt);
            vi->P() = coord;
            MyMesh::VertexPointer v0 = get<0>(*splittingFaceIt);
            MyMesh::VertexPointer v1 = get<1>(*splittingFaceIt);
            MyMesh::VertexPointer v2 = get<2>(*splittingFaceIt);
            printf("before update index: %lu %lu %lu %lu\n", tri::Index(*m, vi.base()), tri::Index(*m, v0),tri::Index(*m, v1),tri::Index(*m, v2));
            if(vpu.NeedUpdate()){
                vpu.Update(v0);
                vpu.Update(v1);
                vpu.Update(v2);
            }
            printf("after update index: %lu %lu %lu %lu\n", tri::Index(*m, vi.base()), tri::Index(*m, v0),tri::Index(*m, v1),tri::Index(*m, v2));
            newFaces.emplace_back(triVert(vi.base(), v0, v1));
            newFaces.emplace_back(triVert(vi.base(), v1, v2));

            MyMesh::CoordType mat0 = vert_mat[v0];
            MyMesh::CoordType mat1 = vert_mat[v2];
            MyMesh::CoordType mat_projector = (mat1 - mat0).normalized();
            vert_mat[vi] = mat0 + mat_projector * t;
            isSplitted[vi] = true;

            vi++;
            splittingFaceIt++;
            coordsIt++;
        }
        assert(newFaces.size() == 2*newCoords.size());
        //create all new faces at once
        MyMesh::FaceIterator  fi = tri::Allocator<MyMesh>::AddFaces(*m, newFaces.size());
        for(auto nfi = newFaces.begin(); nfi != newFaces.end() && fi!=m->face.end(); nfi++, fi++){
            fi->V(0) = get<0>(*nfi);
            fi->V(1) = get<1>(*nfi);
            fi->V(2) = get<2>(*nfi);
            int i0, i1, i2;
            i0 = tri::Index(*m, get<0>(*nfi));
            i1 = tri::Index(*m, get<1>(*nfi));
            i2 = tri::Index(*m, get<2>(*nfi));
            printf("%d %d %d\n", i0,i1,i2);
//            assert(i0 < m->vn);
//            assert(i1 < m->vn);
//            assert(i2 < m->vn);
            for(int j = 0; j< fi->VN(); j++){
                assert(fi->V(j) != NULL);
            }
            fi->C() = Color4b::Red;
        }

        printf("splitted %lu edges\n", newCoords.size());
        tri::UpdateTopology<MyMesh>::VertexFace(*m);
//        tri::UpdateTopology<MyMesh>::FaceFace(*m);
    }

    void detect_degeneracies(){
        double length = edge_threshold()/10;
        for (auto vi = m->vert.begin(); vi != m->vert.end(); vi++){
            if (isFixed[vi] || vi->IsD()) continue;

            int counter = 0;
//            vector<MyMesh::VertexPointer> adjacentVertices;
//            vcg::face::VVStarVF<MyMesh::FaceType>(vi.base(), adjacentVertices);
//
//            for(auto adjV = adjacentVertices.begin(); adjV != adjacentVertices.end(); adjV++){
//                double edge_length = Distance(vi->P(), adjV.operator*()->P());
//                if(edge_length <= length){ //here should be another constraints which check if is collapsable, but not done
//                    counter++;
//                }
//            }

            vcg::face::VFIterator<MyFace> vfi(vi.base()); //initialize the iterator tohe first face
            for(;!vfi.End();++vfi){
                if(vfi.F()->IsD()) continue;

                int vn = vfi.F()->VN();
                MyMesh::VertexPointer V1 = vfi.F()->V((vfi.I()+1)%vn);
                MyMesh::VertexPointer V2 = vfi.F()->V((vfi.I()+2)%vn);
                if( Distance(vi->P(), V1->P()) <= length ||
                    Distance(vi->P(), V2->P()) <= length){
                    counter++;
                }
            }

            if (counter >= 2) {
                isFixed[vi] = true;
            } else {
                isFixed[vi] = false;
            }
        }
    }
};


#endif //MCFSKET_COLLAPSER_H
