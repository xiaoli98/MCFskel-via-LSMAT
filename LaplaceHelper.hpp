//
// Created by boy on 27/11/22.
//

#ifndef MCF_SKET_LAPLACEHELPER_HPP
#define MCF_SKET_LAPLACEHELPER_HPP

#include <unordered_map>
#include <math.h>
#include <tuple>
#include "utils.hpp"


class LaplaceHelper{
    typedef vector<tuple<MyMesh::VertexType*, MyMesh::VertexType*, double>> laplacian_triple;
private:
    MyMesh *m;
//    MyMesh::VertContainer *vert;
//    MyMesh::FaceContainer *faces;
    laplacian_triple laplacian;
    laplacian_triple laplacian_weights;

    MyMesh::PerMeshAttributeHandle<laplacian_triple> laplacian_handle;
public:
    LaplaceHelper(MyMesh *m){
        this->m = m;
    }

    const vector<tuple<MyMesh::VertexType*, MyMesh::VertexType*, double>> &getLaplacian() const {
        return laplacian;
    }

    const vector<tuple<MyMesh::VertexType*, MyMesh::VertexType*, double>> &getLaplacianWeights() const {
        return laplacian_weights;
    }

    double compute_area(Point3d p, Point3d q, Point3d r, double angle){
        double pq = Distance(p, q);
        double pr = Distance(p, r);
        return 0.5 * pq * pr * sin(angle);
    }

    void compute_laplace() {
        int n_vert = m->vert.size();
        laplacian.clear();
        laplacian.reserve(6 * n_vert);
//        cout << "reserving: "<<6* n_vert<<endl;
//        laplacian_weights.resize(n_vert);
        tri::UpdateTopology<MyMesh>::VertexFace(*m);
//        tri::UpdateTopology<MyMesh>::FaceFace(*m);
//        MyMesh::VertexIterator vi;
//        MyMesh::FaceIterator fi;
//        int visited_vert = 0;

        ///creating mapping of vertex and the index
        //todo -> risolto, la copia di vettori e' il male
        //  non ben capito ma i puntatori nei container di vertex sono diversi dai puntatori dei vertici nelle facce
        //  quindi inizializzo la mappa con usando i puntatori dei vertici delle facce
        /// TROPPO PYTHON HA ROVINATO LA NUOVA GENERAZIONE

        MyMesh::VertexType *vert0_p, *vert1_p;
//        tri::UpdateFlags<MyMesh>::VertexClearV(*m);

#if COLLAPESER_DEBUG
        int nulls = 0;
        for(auto vit = m->vert.begin(); vit != m->vert.end(); vit++){
            if(vit->VFp() == NULL){
                printf("%p\t", vit.base());
                nulls++;
            }
        }
        cout <<"null vertexface pointers: "<<nulls<<endl;
#endif

        // todo
        // Pos is not suitable for mesoskeletons
        // this is not manifold, so the navigation is not good
        //a possible idea is to use only VF an FF adjacency
        for (auto vit = m->vert.begin(); vit != m->vert.end(); vit++){
            if (vit->IsD()) continue;
//            vcg::face::Pos<MyFace> p(start, vit.base());
            double sin_alpha, cos_alpha, alpha, cot_alpha;
            double sin_beta, cos_beta, beta,cot_beta;
            double w = 0, sum_area=0, sum = 0;
//------------------------------------------------------------------------------------------------------------------------------
            vector<tuple<MyMesh::FacePointer, int>> oneRingFaces;
            MyMesh::FacePointer leftFace, rightFace;
            face::VFIterator<MyFace> vfi(vit.base());
            face::JumpingPos<MyMesh::FaceType> jpos(vit->VFp(), vit.base());

            MyMesh::FaceType *start = vit->VFp();
            do{
//                oneRingFaces.emplace_back(tuple<MyMesh::FacePointer, int>(vfi.F(), vfi.I()));
                oneRingFaces.emplace_back(tuple<MyMesh::FacePointer, int>(jpos.f, jpos.z));
                jpos.NextFE();
            }while(jpos.f != start);


            int fn = oneRingFaces.size();
            int leftidx, leftI, rightidx, rightI;

            for(int i = 0; i < oneRingFaces.size()+1; i++){
                rightidx = i%fn;
                rightFace = get<0>(oneRingFaces[rightidx]);
                rightI = get<1>(oneRingFaces[rightidx]);

                leftidx = (i+1)%fn;
                leftFace = get<0>(oneRingFaces[leftidx]);
                leftI = get<1>(oneRingFaces[leftidx]);

                alpha = Angle(rightFace->V(rightFace->Prev(rightI+1))->cP() - vit->cP(), rightFace->V(rightFace->Next(rightI+1))->cP() - vit->cP());
                sincos(alpha, &sin_alpha, &cos_alpha);
                cot_alpha = cos_alpha / sin_alpha;
                beta = Angle(leftFace->V(leftFace->Prev(leftI-1))->cP() - vit->cP(), leftFace->V(leftFace->Next(leftI-1))->cP() - vit->cP());
                sincos(beta, &sin_beta, &cos_beta);
                cot_beta = cos_beta / sin_beta;

                w = (cot_alpha + cot_beta) * 0.5;

                vert0_p = vit.base();
                vert1_p = rightFace->V(rightidx-1);
                laplacian.emplace_back(make_tuple(vert0_p, vert1_p, w));
                sum += w;
            }

            for(auto f_i = oneRingFaces.begin(); f_i != oneRingFaces.end(); f_i++){
                MyMesh::FacePointer f = get<0>(*f_i);
                int idx = get<1>(*f_i);
                double angle = Angle(f->V(f->Prev(idx))->cP() - vit->cP(), f->V(f->Next(idx))->cP() - vit->cP());
                sum_area += compute_area(f->V(idx)->cP(), f->V(f->Prev(idx))->cP(), f->V(f->Next(idx))->cP(), angle);
            }

//------------------------------------------------------------------------------------------------------------------------------
#pragma region
//            do
//            {
//                double area = 0, angle;
//                Point3d P, Q, R; //3 vertices to be save for area computation
//                P = p.V()->P();
//                vert0_p = p.V();
//
//                //IN THE PREV EDGE
//                //go to the opposite vert on the same edge
//                p.FlipV();
//                alpha = p.AngleRad();
//                sincos(alpha, &sin_alpha, &cos_alpha);
//                cot_alpha = cos_alpha / sin_alpha;
//                Q = p.V()->P();
//                p.FlipV();//return to the original vert
//
//                //go to the NEXT EDGE and the other vertex of the edge
//                p.FlipE();
//                p.FlipV();
//                R = p.V()->P();
//                vert1_p = p.V();
//
//                p.FlipV();
//
//                angle = p.AngleRad();
//                //go to the adjacent face
//                while(p.FFlip() == NULL)
//                    p.FlipF();
//                p.FlipF();
//                p.FlipE();
//
//                //go to the opposite vert on the same edge
//                p.FlipV();
//                beta = p.AngleRad();
//                sincos(beta, &sin_beta, &cos_beta);
//                p.FlipV();
//                cot_beta = cos_beta / sin_beta;
//                w = (cot_alpha + cot_beta) * 0.5;
//
//                area = compute_area(P, Q, R, angle);
//                sum += w;
//                sum_area += area;
//                //inserting <i, j, w>
//                laplacian.emplace_back(make_tuple(vert0_p, vert1_p, w));
//
//                //NEXT EDGE becomes CURRENT EDGE
//                p.FlipE();
//            }while(p.f!=start);
#pragma endregion
            laplacian.emplace_back(make_tuple(vert0_p, vert0_p, -sum));
            laplacian_weights.emplace_back(make_tuple(vert0_p, vert0_p, 1 / (sum_area / 3)));
        }

        laplacian_handle = tri::Allocator<MyMesh>::GetPerMeshAttribute<laplacian_triple>(*m, string("laplacian"));
        laplacian_handle() = laplacian;
    }

    void print_lapacian(){
        cout << "laplacian size:"<<laplacian.size() <<endl;
        for(auto item = laplacian.begin(); item != laplacian.end(); item++){
            printf("row: %lu col:%lu value:%f\n", tri::Index(*m,get<0>(*item)), tri::Index(*m, get<1>(*item)), get<2>(*item));
        }
    }
};
#endif //MCF_SKET_LAPLACEHELPER_HPP
