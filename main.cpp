#include <iostream>

#include "meanCurvatureFlow.hpp"
#include "sphereShrinking.hpp"
#include "LSMAT.hpp"
#include "utils.hpp"
#include <vcg/space/index/kdtree/kdtree.h>

using namespace vcg;
using namespace std;

int main(int argc, char* argv[]){
    MyMesh m;

    tri::io::ImporterOFF<MyMesh>::Open(m,argv[1]);
    tri::UpdateBounding<MyMesh>::Box(m);
    tri::RequirePerVertexNormal(m);
    tri::RequireVFAdjacency(m);
    tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);

    cout << "mesh name: " <<argv[1] <<endl;
    cout << "FN: "<< m.FN() << " VN: "<<m.VN()<<endl;

#if LS_MAT
    LSMAT lsmat(&m);
    for(int j = 0; j < 10; j++){
        cout << "*****************LSMAT iteration "<< j <<"******************"<<endl;
        string filename = to_string(*basename(argv[1])) + "_LSMAT_"+ to_string(j)+".off";
        vector<Sphere3d> spheres = lsmat.compute(5);

//        for(int i = 0; i<m.VN(); i++){
//            m.vert[i].P() = spheres[i].Center();
//        }
//        tri::io::ExporterOFF<MyMesh>::Save(m, filename.c_str(), tri::io::Mask::IOM_FACECOLOR);
    }
#endif

#if SPHERE_SHRINK
    VertexConstDataWrapper<MyMesh> ww(m);
    KdTree<float> tree(ww);
    SphereShrinking ss = SphereShrinking(&m, tree);
    ss.compute_ma_point();
    vector<Sphere3d> medial = ss.getMedialSpheres();
#if DEBUG
    for(auto m = medial.begin(); m < medial.end(); m++){
        cout<<m->Radius()<<"\t\t";
        PRINTP(m->Center())
    }
#endif
    int i= 0;
    for(auto vert = m.vert.begin(); vert != m.vert.end(); vert++){
        vert->P() = medial[i].Center();
        i++;
    }
    string outFileName = strcat(basename(argv[1]), "_SS.off") ;
    tri::io::ExporterOFF<MyMesh>::Save(m, outFileName.c_str(), tri::io::Mask::IOM_FACECOLOR);
#endif

#if MCF
    MeanCurvatureFlow mcf(&m);
    mcf.compute_skel();
#endif
}

