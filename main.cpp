#include <iostream>

#include "RboxPoints.h"
#include "QhullError.h"
#include "QhullQh.h"
#include "QhullFacet.h"
#include "QhullFacetList.h"
#include "QhullFacetSet.h"
#include "QhullRidge.h"
#include "QhullLinkedList.h"
#include "QhullVertex.h"
#include "Qhull.h"

#include "sphereShrinking.hpp"
#include "LSMAT.hpp"
#include "utils.hpp"

using namespace vcg;
using namespace std;
using namespace orgQhull;

class MAT{
public:
//    MyMesh m;
//    ScalarVertexProperty vangle;
//    ScalarVertexProperty vradii;
//    StarlabDrawArea* drawArea;

    std::vector< Point3f > loci;
    std::vector< std::vector<uint> > cells;

    std::vector<int> poleof;				// The index of the voronoi pole the vertex refers to
    std::vector< std::vector<int> > scorr; 	// The index the surface sample to which a pole corresponds to

    std::vector<double> alpha;
    std::vector<double> radii;

    coordT *points;
    coordT *normals;

    int nvertices{}, nvornoi{};

    MAT(){}

    coordT *getPoints(MyMesh &m, int dimension){
        int cnt=0;
        coordT *coords;
        coords = this->points = this->normals = (coordT*)malloc((m.VN())*(dimension)*sizeof(coordT));

        MyMesh::VertexIterator vi;
        for(vi = m.vert.begin(); vi != m.vert.end(); ++vi){
//        cout <<"Point coords: ";
            for(int ii = 0;ii < dimension; ++ii){
                *coords = (*vi).P()[ii];
                *normals = (*vi).N()[ii];
                coords++;
                this->normals++;
            }
//        cout <<endl;
            ++cnt;
//        cout << "is in : "<<m.bbox.IsIn((*vi).P())<<endl;
        }
        assert(cnt==m.VN());
        return points;
    }

    void searchPoles(qhT *qh, int nvertices, int nvoronoi){
        std::vector<int> poleof(nvertices, 0);
        std::vector< std::vector<int> > scorr(nvoronoi, std::vector<int>(4, 0));
        std::vector<int> counter(nvoronoi, 0);

        MyVertex p;
        vertexT *vertex;
        int dim = 3;
        uint sidx = 0;
        FORALLvertices {
            double max_dist = 0;
            double max_dist_i = 0;

//            vertex->;

            Point3f voro_vertex;
            for(int j = 0; j < (int)cells[sidx].size(); j++){
                int vidx = cells[sidx][j];
                voro_vertex = loci[vidx];

                int freesub = counter[vidx];
                if( freesub < 4 ){
                    counter[vidx]++;
                    scorr[vidx][freesub] = sidx;
                }

                double dist = qh_pointdist(vertex->point, (double*)voro_vertex.V(), dim);
                if(dist > max_dist){
                    max_dist = dist;
                    max_dist_i = vidx;
                }
            }
            poleof[sidx] = max_dist_i;
            sidx++;
        }
    }

    void getMedialSpokeAngleAndRadii(int nvoronoi, coordT *points)
    {
        std::vector<double> alpha;
        std::vector<double> radii;
        //--- Temp data
        Point3f surf_vertex;
        Point3f voro_vertex;
        Point3f s; // temp spoke data
        double curralpha;

        alpha = std::vector<double> ( nvoronoi );
        radii = std::vector<double> ( nvoronoi );

        //--- For every voronoi vertex
        for(int vidx = 0; vidx < nvoronoi; vidx++)
        {
            /// Do not use invalid poles
            //if( ispole.get(vidx) == 0.0 ) continue;

            // Retrieve vertex coordinates
            voro_vertex = loci[vidx];

            // Create the medial spokes
            // General positions => only 4 vertices / loci
            std::vector<Point3f> spokes;
            for(int i_sidx = 0; i_sidx < 4; i_sidx++)
            {
                // Retrieve surface coordinate
                int sidx = scorr[vidx][i_sidx];
//                surf_vertex = points[Vertex(sidx)];

                // Create spoke
                s = surf_vertex.operator-(voro_vertex);

                // Spoke length (shouldn't this be same as we are voronoi loci?)
                radii[vidx] = s.Norm();

                // Normalize spoke
                s.normalized();
                spokes.push_back(s);
            }

            // Measure largest spoke aperture
            // and store it in output
            double alpha_max = 0;
            for(int i=0; i<4; i++){
                Point3f& s1 = spokes[i];
                for(int j=0; j < 4; j++){
                    Point3f& s2 = spokes[j];
                    curralpha = Angle(s1,s2);
                    if( curralpha > alpha_max ){
                        alpha_max = curralpha;
                    }
                }
            }
            alpha[vidx] = alpha_max;
        }
    }
};

void searchFirstPole(qhT *qh, int dim, MyMesh &m, double threshold){
    qh_setvoronoi_all(qh);
    vector<double*> voronoi_vertices;
    vector<double*> normals;
    setT* poles_set= qh_settemp(qh, qh->num_facets);

    vertexT *vertex;
    FORALLvertices {
        double *voronoi_vertex;
        double *first_pole;
        voronoi_vertices.clear();
        normals.clear();
        double max_dist=0; //distance from first_pole to vertex

        if (qh->hull_dim == 3)
            qh_order_vertexneighbors(qh, vertex);

        bool is_on_convexhull =false;

        //Considering the neighboring facets of the vertex in order to
        //compute the Voronoi region for that vertex
        facetT *neighbor, **neighborp;

        //Finding first_pole.
        FOREACHneighbor_(vertex) {
            if (neighbor->upperdelaunay)
                is_on_convexhull =true;
            else {
                voronoi_vertex = neighbor->center;
                voronoi_vertices.push_back(neighbor->center);

                if (neighbor->toporient)
                    normals.push_back(neighbor->normal);
                double* vertex1= vertex->point;
                double* vertex2= voronoi_vertex;
                double dist = qh_pointdist(vertex1, vertex2,dim);
                if(dist>max_dist){
                    max_dist=dist;
                    first_pole=voronoi_vertex;
                }
            }
        }

        if(is_on_convexhull){ //This is a Voronoy vertex at infinity
            //Compute normals average
            double* avg_normal= new double[3];
            for(int i=0;i<3;i++){
                avg_normal[i]=0;
                for(size_t j=0; j< normals.size();j++)
                    avg_normal[i]+= (normals[j])[i];
            }
            first_pole=avg_normal;

//#if(TestMode)
//            firstTest << first_pole[0]<< " " << first_pole[1]<<" " << first_pole[2]<<"\n ";
//#endif

        }
        else{
            assert(first_pole!=NULL);
            if(first_pole!=NULL){
                //Test if the Voronoi vertex is too far from the origin.
                //It can cause problems to the qhull library.
                bool discard=false;
                for(int i =0;i<3;i++)
                {
                    double* bbCenter = new double[dim];
                    double* pole = first_pole;
                    for(int i=0;i<dim;i++)
                        bbCenter[i] = m.bbox.Center()[i];
                    if(qh_pointdist(bbCenter,pole,dim)>(threshold*m.bbox.Diag()))
                        discard=true;
                }
                if(!discard)
                    qh_setunique(qh, &poles_set, first_pole);
            }
        }

        //vector vertex-first_pole
        double* sp1= new double[3];
        for(int i=0; i<3;i++)
            sp1[i]= first_pole[i] - vertex->point[i];
    }
}

//int main( int argc, char **argv )
//{
//    if(argc<2)
//    {
//        printf("Usage trimesh_base <meshfilename.off>\n");
//        return -1;
//    }
//
//    MyMesh m;
//
//    MyMesh m2;
//    int loadmask = 0;
////    tri::io::ImporterOBJ<MyMesh>::Open(m,argv[1], loadmask);
//    tri::io::ImporterOFF<MyMesh>::Open(m,argv[1]);
//    Point3f p;
//    printf("Input mesh '%s'  vn:%i fn:%i\n",argv[1],m.VN(),m.FN());
//
//    tri::Allocator<MyMesh>::CompactVertexVector(m);
//    tri::UpdateBounding<MyMesh>::Box(m);
//
//    char flags[]="v Qbb";
//    char comment[]="";
//    int dimension = 3;
//    MAT mat;
//    coordT *points = mat.getPoints(m, dimension);
//
//    Qhull myqhull(comment,dimension,m.VN(), points, flags);
//
//    std::vector< std::vector<uint>> cells;
//    std::vector<Point3f> loci;
//
//    int i = 0;
//    cout << "facetList size: "<<myqhull.facetList().size()<<endl;
//    for(QhullFacet f: myqhull.facetList())
//    {
////        cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
//        //SE E' DELAUNAY ALLORA E' UN VERTICE DEL VORONOI
//        if(f.isUpperDelaunay()) {
////            cout <<"isDelaunay"<<endl;
//            continue;
//        }
//        QhullPoint voronoiPoints = f.voronoiVertex();
////        cout << "VORONOI P "<<*voronoiPoints.coordinates()<<endl;
////        cout<<"added point to new mesh"<<endl;
//        Point3f p(voronoiPoints[0],voronoiPoints[1],voronoiPoints[2]);
//        tri::Allocator<MyMesh>::AddVertex(m2, p, Color4b::Red);
//
//        if(!m.bbox.IsIn(p)) {
////            PRINTP(p);
////            cout<<"not in the bounding box"<<endl;
//            continue;
//        }
////        cout <<"here"<<endl;
//        loci.emplace_back(p);
//
//        cells.resize(m.VN());
//        for(QhullVertex v: f.vertices())
//            cells[v.point().id()].push_back(i);
//        i++;
//    }
//
//    cout << "loci size: "<<loci.size()<<endl;
////    for(Point3f temp: loci){
////        PRINTP(temp);
////    }
//
//    int nvertices = m.VN();
//    int nvoronoi = loci.size();
//
//    auto angle = [](Point3f a, Point3f b)->double {
//        return acos(a.X()*b.X() + a.Y()*b.Y() + a.Z()*b.Z());
//    };
//
//    //TODO compute medial axis
//        //get medial angle and radii
//        //set surfaces to medial
////    tri::io::ExporterOBJ<MyMesh>::Save(m,"cactusTest.obj",tri::io::Mask::IOM_FACECOLOR);
//    cout << "mesh2 vertices: "<< m2.VN()<<endl;
//    tri::io::ExporterOFF<MyMesh>::Save(m2,"torusTest.off",tri::io::Mask::IOM_FACECOLOR);
//
//    return 0;
//}


int main(int argc, char* argv[]){
    MyMesh m;

    tri::io::ImporterOFF<MyMesh>::Open(m,argv[1]);
    tri::UpdateBounding<MyMesh>::Box(m);
    tri::RequirePerVertexNormal(m);
    tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);

//    LSMAT lsmat(&m);
//    for(int j = 0; j < 10; j++){
//        cout << "*****************LSMAT iteration "<< j <<"******************"<<endl;
//        string filename = to_string(*basename(argv[1])) + "_LSMAT_"+ to_string(j)+".off";
//        vector<Sphere3d> spheres = lsmat.compute(1);
//
//        for(int i = 0; i<m.VN(); i++){
//            m.vert[i].P() = spheres[i].Center();
//        }



//        tri::io::ExporterOFF<MyMesh>::Save(m, filename.c_str(), tri::io::Mask::IOM_FACECOLOR);
//    }
    //TODO
    //  esportare le sfere con addSphere()
    //  mappare le sfere con i colori

    SphereShrinking ss = SphereShrinking(&m);
//    for(auto point = m.vert.begin(); point < m.vert.end(); point++){
//        point->P() = ss.compute_ma_point(point->P(), point->N());
//    }

    ss.compute_ma_point();
    vector<Sphere3d> medial = ss.getMedialSpheres();

    for(auto m = medial.begin(); m < medial.end(); m++){
        cout<<m->Radius()<<"\t\t";
        PRINTP(m->Center())
    }

    int i = 0;
    for(auto vert = m.vert.begin(); vert != m.vert.end(); vert++){
        if(medial[i].Radius() > 0)
            vert->P() = medial[i].Center();
        i++;
    }
    string outFileName = strcat(basename(argv[1]), "_SS.off") ;
    tri::io::ExporterOFF<MyMesh>::Save(m, outFileName.c_str(), tri::io::Mask::IOM_FACECOLOR);

}

