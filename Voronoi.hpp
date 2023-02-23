#ifndef MCFSKET_VORONOI_H
#define MCFSKET_VORONOI_H

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
#include "utils.hpp"

using namespace orgQhull;
using namespace vcg;
using namespace std;
class Voronoi_MAT{
public:
    MyMesh *m;
//    ScalarVertexProperty vangle;
//    ScalarVertexProperty vradii;
//    StarlabDrawArea* drawArea;

    std::vector< Point3f > loci;
    std::vector< std::vector<uint> > cells;

    std::vector<int> poleof;				// The index of the voronoi pole the vertex refers to
    std::vector< std::vector<int> > scorr; 	// The index the surface sample to which a pole corresponds to

    std::vector<double> alpha;
    std::vector<double> radii;

    MyMesh::PerVertexAttributeHandle<float> vangle;
    MyMesh::PerVertexAttributeHandle<float> vradii;

    coordT *points;
    coordT *normals;

    int nvertices{}, nvornoi{};

    Voronoi_MAT(MyMesh *m){
        this->m = m;
        vangle = tri::Allocator<MyMesh>::GetPerVertexAttribute<float>(*m, string("vangle"));
        vradii = tri::Allocator<MyMesh>::GetPerVertexAttribute<float>(*m, string("vradii"));
    }

    void computeVoronoiDiagram(){
        points = getPoints();
        Qhull qhull("", 3, m->VN(), points, "v Qbb");
        //QhullFacetList facets= qhull.facetList();
        //std::cout << facets;

        // Extract vertices & cell structure
        cells.resize( m->VN() );
        int i = 0;

        for(QhullFacet f: qhull.facetList())
        {
            if(f.isUpperDelaunay()) continue;

            QhullPoint voronoiPoints = f.voronoiVertex();
            Point3f p(voronoiPoints[0],voronoiPoints[1],voronoiPoints[2]);

            if(!m->bbox.IsIn(p)) continue;

            loci.emplace_back(p);

            for(QhullVertex v: f.vertices())
                cells[v.point().id()].emplace_back(i);

            i++;
        }

        // Used later
        nvertices = this->m->VN();
        nvornoi = loci.size();
    }

    coordT *getPoints(int dimension=3){
        int cnt=0;
        coordT *coords;
        coords = this->points = this->normals = (coordT*)malloc((m->VN())*(dimension)*sizeof(coordT));

        MyMesh::VertexIterator vi;
        for(vi = m->vert.begin(); vi != m->vert.end(); ++vi){
            for(int ii = 0;ii < dimension; ++ii){
                *coords = (*vi).P()[ii];
                *normals = (*vi).N()[ii];
                coords++;
                this->normals++;
            }
            ++cnt;
        }
        assert(cnt==m->VN());
        return points;
    }

    struct Spoke{
        double x;
        double y;
        double z;

        Spoke(){ x=y=z=0; }
        Spoke( double x, double y, double z ){
            this->x = x;
            this->y = y;
            this->z = z;
        }

        double angle( const Spoke& s ){
            return acos( (this->x)*(s.x)+(this->y)*(s.y)+(this->z)*(s.z) );
        }
    };

    void searchVoronoiPoles()
    {
        Point3f surf_vertex;
        Point3f voro_vertex;
        Point3f surf_normal;

        poleof = std::vector<int>(nvertices, 0);
        scorr  = std::vector< std::vector<int> > (nvornoi, std::vector<int>(4, 0));

        //--- Keeps track of how many
        //    surface points a voronoi loci has been
        //    associated with. Assuming general positions
        //    this number should be always 4.
        std::vector<int> counter(nvornoi, 0);

        uint sidx;
        //--- For every voronoi cell
        for(auto vit = m->vert.begin(); vit != m->vert.end(); vit++)
        {
            // Retrieve surface vertex
            sidx = tri::Index(*m, vit.base());
            surf_vertex = vit->cP();
            surf_normal = vit->N();

            // Index and distance to furthest voronoi loci
            double max_neg_t = DBL_MAX;
            double max_neg_i = 0;

            // For each element of its voronoi cell
            for(int j = 0; j < (int)cells[sidx].size(); j++)
            {
                int vidx = cells[sidx][j];

                voro_vertex = loci[vidx];

                // Mark the fact that (voronoi) vidx corresponds to (surface) sidx
                // (in the next free location and update this location)
                // the freesub-th correspondent of the voronoi vertex is vidx
                int freesub = counter[vidx];
                if( freesub < 4 ){
                    counter[vidx]++;
                    scorr[vidx][freesub] = sidx;
                }

                // Project the loci on the vertex normal & Retain furthest
                double t = Point3d(voro_vertex - surf_vertex) * surf_normal;
                if(t < 0 && t < max_neg_t){
                    //drawArea->drawPoint(voro_vertex);
                    max_neg_t = t;
                    max_neg_i = vidx;
                }
            }

            // Save pole to which surface corresponds
            // Store index (possibly nan!! into buffer)
            poleof[sidx] = max_neg_i;
        }
    }

//    void searchPoles(qhT *qh, int nvoronoi){
//        std::vector<int> poleof(nvertices, 0);
//        std::vector< std::vector<int> > scorr(nvoronoi, std::vector<int>(4, 0));
//        std::vector<int> counter(nvoronoi, 0);
//
//        MyVertex p;
//        vertexT *vertex;
//        int dim = 3;
//        uint sidx = 0;
//        FORALLvertices {
//                double max_dist = 0;
//                double max_dist_i = 0;
//
//                Point3f voro_vertex;
//                for(int j = 0; j < (int)cells[sidx].size(); j++){
//                    int vidx = cells[sidx][j];
//                    voro_vertex = loci[vidx];
//
//                    int freesub = counter[vidx];
//                    if( freesub < 4 ){
//                        counter[vidx]++;
//                        scorr[vidx][freesub] = sidx;
//                    }
//
//                    double dist = qh_pointdist(vertex->point, (double*)voro_vertex.V(), dim);
//                    if(dist > max_dist){
//                        max_dist = dist;
//                        max_dist_i = vidx;
//                    }
//                }
//                poleof[sidx] = max_dist_i;
//                sidx++;
//        }
//    }

    void getMedialSpokeAngleAndRadii()
    {
        //--- Temp data
        Point3f surf_vertex;
        Point3f voro_vertex;
        Point3f s; // temp spoke data
        double curralpha;

        alpha = std::vector<double> ( nvornoi );
        radii = std::vector<double> ( nvornoi );

        //--- For every voronoi vertex
        for(int vidx = 0; vidx < nvornoi; vidx++)
        {
            /// Do not use invalid poles
            //if( ispole.get(vidx) == 0.0 ) continue;

            // Retrieve vertex coordinates
            voro_vertex = loci[vidx];

            // Create the medial spokes
            // General positions => only 4 vertices / loci
            std::vector<Spoke> spokes;
            for(int i_sidx = 0; i_sidx < 4; i_sidx++)
            {
                // Retrieve surface coordinate
                int sidx = scorr[vidx][i_sidx];
                surf_vertex = m->vert[sidx].cP();

                // Create spoke
                s = surf_vertex - voro_vertex;

                // Spoke length (shouldn't this be same as we are voronoi loci?)
                radii[vidx] = s.Norm();

                // Normalize spoke
                s.normalize();
                spokes.push_back( Spoke(s[0],s[1],s[2]) );
            }


            // Measure largest spoke aperture
            // and store it in output
            double alpha_max = 0;
            for(int i=0; i<4; i++){
                Spoke& s1 = spokes[i];
                for(int j=0; j < 4; j++){
                    Spoke& s2 = spokes[j];
                    curralpha = s1.angle(s2);
                    if( curralpha > alpha_max ){
                        alpha_max = curralpha;
                    }
                }
            }
            alpha[vidx] = alpha_max;
        }
    }
    void setMedialAttribute(){
        MyMesh::PerVertexAttributeHandle<MyMesh::CoordType> vert_mat = tri::Allocator<MyMesh>::GetPerVertexAttribute<MyMesh::CoordType>(*m, string("vert_mat"));
        for(auto vit = m->vert.begin(); vit != m->vert.end(); vit++){
            vert_mat[vit.base()] = loci[poleof[tri::Index(*m, vit.base())]];
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
#endif //MCFSKET_VORONOI_H

