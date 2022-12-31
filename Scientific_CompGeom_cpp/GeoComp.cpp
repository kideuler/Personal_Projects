#include "GeoComp.hpp"
#include <unistd.h>
using namespace std;

/**
 * @brief Converts half-facet id to element id
 * 
 * @param hfid Half facet id
 * @return element id
 */
int hfid2eid(int hfid){
    return (hfid >> 8);
}
/**
 * @brief Convert half-facet id to local edge id
 * 
 * @param hfid Half facet id
 * @return local edge id
 */
int hfid2lid(int hfid){
    return (hfid&255) + 1;
}
/**
 * @brief Converts element id and local edge id to half-facet id
 * 
 * @param eid element id
 * @param lid local edge id
 * @return half-facet id
 */
int elids2hfid(int eid, int lid){
    return ((eid << 8) + lid - 1);
}

/// subfunctions for delaunay triangulation
double eval_alpha(const vector<vector<double>> xs,double r_ref);
static vector<double> circumcenter(const vector<vector<double>> xs);
static bool inside_tri(const Mat &xs, const vector<double> &ps);
static bool inside_circumtri(const Mat xs, const vector<double> ps);
static double area_tri(const Mat &xs);
void find_bad_tri_recursive(Triangulation* DT, int tri, int *nbad, int vid);
void find_bad_tri_recursive(Triangulation* DT, int tri, int *nbad, int vid, bool* exitf);
void find_bad_tri(Triangulation* DT, int tri, int *nbad, int vid);
vector<int> facet_reorder(Triangulation *DT, int *nsegs);


/// array functions for point data
static vector<double> min_array(const vector<vector<double>> &xs);
static vector<double> max_array(const vector<vector<double>> &xs);
static vector<double> find_center(const vector<vector<double>> &xs);
void reorder(vector<vector<double>> &xs);

/**
 * @brief Refine a delaunay triangulation using Chews second algorithm
 * 
 * @param DT Triangultion data structure
 * @param r_ref radius of circumcircle (double)
 */
void GeoComp_refine(Triangulation* DT, double r_ref){
    auto f = [r_ref](vector<double> xs) {return r_ref; };
    GeoComp_refine(DT, f);
}

/**
 * @brief Refine a delaunay triangulation using Chews second algorithm
 * 
 * @param DT Triangultion data structure
 * @param r_ref radius of circumcircle (function of position)
 */
void GeoComp_refine(Triangulation* DT, function<double(vector<double>)> r_ref){
    // random number generator
    default_random_engine re;
    uniform_real_distribution<double> unif(-1, 1);

    int n,i;
    bool exitl;
    double alpha;
    int nv = (*DT).coords.size();
    Mat ps = Zeros(3,2);
    int tri;
    (*DT).coords.resize((*DT).elems.size());
    vector<double> C;

    // main loop 
    n=0;
    while (n<(*DT).nelems){
        if (!(*DT).delete_elem[n]){
            ps[0] = (*DT).coords[(*DT).elems[n][0]];
            ps[1] = (*DT).coords[(*DT).elems[n][1]];
            ps[2] = (*DT).coords[(*DT).elems[n][2]];
            C = circumcenter(ps);
            // wiggle circumcenter for stability
            C[0] += 1e-5*unif(re);
            C[1] += 1e-5*unif(re);
            // sanity check
            if (abs(C[0] > 1e6 || abs(C[1]) > 1e6)){
                cout << "abnormally large point added from points " << C[0] << "," << C[1] << endl;
                printMat(ps); 
            }
            alpha = eval_alpha(ps,r_ref(C));
            if (alpha > 1.2){
                // add circumcircle to triangulation
                (*DT).coords[nv] = C;
                Bowyer_watson2d(DT,nv,n,true);
                nv++;

                if ((double) DT->nelems >= (0.8)*((double) DT->elems.size())){
                    cout << "approaching size bound, freeing up space" << endl;
                    delete_tris(DT,&n);
                }
            }
        }
        n++;
    }

    (*DT).coords.resize(nv);
    delete_tris(DT);
    cout << "created " << (*DT).nelems << " triangles from refining the mesh" << endl;
}

/**
 * @brief Create constrained Delaunay triangulation from PSLG
 * 
 * @param segs Constrained segments (nsegs -by- 2)
 * @param xs Point data (nv -by- 2)
 * @return Constrained delaunay triangulation
 */
Triangulation GeoComp_Delaunay_Triangulation(const vector<vector<int>> &segs, vector<vector<double>> &xs){

    // segs define boundary segments for the mesh
    Triangulation DT = GeoComp_Delaunay_Triangulation(xs);
    int nv = xs.size();
    vector<vector<int>> onering(nv);
    vector<int> numonering(nv);

    int vid,i,j,hfid;
    for (i = 0; i<DT.nelems; i++){
        for (j = 0; j<3; j++){
            vid = DT.elems[i][j];
            onering[vid].resize(numonering[vid]+1);
            onering[vid][numonering[vid]] = elids2hfid(i+1,j+1);
            numonering[vid]++;
        }
    }

    int ns = segs.size();
    int oppeid,eid,lnid;
    bool exitf;
    for (i = 0; i<ns; i++){
        vid = segs[i][0];
        j=0;
        exitf = false;
        while(j<numonering[vid] && !exitf){
            eid = hfid2eid(onering[vid][j])-1;
            lnid = hfid2lid(onering[vid][j])-1;
            if (DT.elems[eid][(lnid+1)%3] == segs[i][1]){
                hfid = DT.sibhfs[eid][(lnid)%3];
                if (hfid != 0){
                    oppeid = hfid2eid(hfid)-1;
                    DT.delete_elem[oppeid] = true;
                }
                exitf = true;
            }
            j++;
        }
        if (!exitf){
            cout << "no constrained edge found for edge: " << segs[i][0] << "," << segs[i][1] << endl;
        }
    }

    delete_tris(&DT);

    return DT;
}

/**
 * @brief create Delaunay triangulation from point set
 * 
 * @param xs Point data (nv -by- 2)
 * @return Delaunay Triangulation 
 */
Triangulation GeoComp_Delaunay_Triangulation(vector<vector<double>> &xs){
    // size checking
    default_random_engine re;
    uniform_real_distribution<double> unif(-1, 1);
    int nv = xs.size();
    int n,i,j;

    Triangulation DT;
    int ub = 2*nv*nv;
    DT.coords = Zeros(nv+3,2);
    DT.elems = Zerosi(ub,3);
    DT.sibhfs = Zerosi(ub,3);
    DT.facets = Zerosi(ub,2);
    DT.bwork.resize(ub);
    DT.delete_elem.resize(ub);
    DT.vedge.resize(ub);
    DT.badtris.resize(ub);
    DT.on_boundary.resize(ub);
    vector<double> a = min_array(xs);
    vector<double> b = max_array(xs);

    // reorder points for optimal triangle search O(N*log(N))
    //reorder(xs);

    for (n=0; n<nv; n++){
        DT.coords[n][0] = xs[n][0] + (b[0]-a[0])*1e-4*unif(re);
        DT.coords[n][1] = xs[n][1] + (b[1]-a[1])*1e-4*unif(re);
    }
    double dx = (b[0]-a[0])/10;
    double dy = (b[1]-a[1])/10;
    a[0] = a[0] - dx;
    a[1] = a[1] - dy;
    b[0] = b[0] + dx;
    b[1] = b[1] + dy;
    DT.coords[nv][0] = a[0];
    DT.coords[nv][1] = a[1];
    DT.coords[nv+1][0] = a[0] + 2*(b[0]-a[0]);
    DT.coords[nv+1][1] = a[1];
    DT.coords[nv+2][0] = a[0];
    DT.coords[nv+2][1] = a[1] + 2*(b[1]-a[1]);

    DT.elems[0] = {nv,nv+1,nv+2};
    DT.sibhfs[0] = {0,0,0};

    // loop inserting each point into triangulation
    DT.nelems = 1;
    int tri = -1;
    bool exitl;
    vector<vector<double>> ps = {{0,0},{0,0},{0,0}};
    for (int n = 0; n < nv; n++){
        i=DT.nelems-1;
        exitl = false;
        while ((i >= 0 ) && !exitl){
            if (!DT.delete_elem[i]){
                ps[0] = DT.coords[DT.elems[i][0]];
                ps[1] = DT.coords[DT.elems[i][1]];
                ps[2] = DT.coords[DT.elems[i][2]];
                if (inside_tri(ps,DT.coords[n])){
                    tri = i;
                    exitl = true;
                }
                i--;
            } else {
                i--;
            }
        }
        // inserting node into the triangulation using Bowyer-Watson algorithm
        Bowyer_watson2d(&DT,n,tri,false);
        if ((double) DT.nelems >= (0.8)*((double) ub)){
            cout << "approaching size bound, freeing up space" << endl;
            delete_tris(&DT);
        }
    }

    for (int n = 0; n<DT.nelems; n++){
        for (int i = 0; i<3; i++){
            if (DT.elems[n][i] >= nv){
                DT.delete_elem[n] = true;
            }
        }
    }

    delete_tris(&DT);
    DT.coords.resize(nv);
    for (n = 0; n<nv; n++){
        DT.coords[n] = xs[n];
    }
    cout << "created " << DT.nelems << " triangles from initial points" << endl;
    return DT;
}

/// mod function for Bowyer-Watson
static int modi(int a, int b){
    int z = a - b*(a/b);
    if (z<0){
        z=z+b;
    }
    return z;
}

/**
 * @brief Add point to triangulation using Bowyer-Watson Algorithm
 * 
 * @param DT Triangulation (Pass by reference)
 * @param vid Point in triangulation to be added
 * @param tri_s Starting triangle
 * @param refine whether the mesh is being refined flag
 */
void Bowyer_watson2d(Triangulation* DT, int vid, int tri_s,bool refine){
    int nbad = 0;
    int i,j;
    vector<vector<double>> ps = {{0,0},{0,0},{0,0}};
    if (refine) {
        bool exitf = false;
        int eid,lid;
        find_bad_tri_recursive(DT, tri_s, &nbad, vid, &exitf);
        if (exitf) {
            eid = hfid2eid(nbad)-1;
            lid = hfid2lid(nbad)-1;
            DT->coords[vid] = (DT->coords[DT->elems[eid][(lid+1)%3]] + DT->coords[DT->elems[eid][lid]])/2;
            nbad = 0;
            find_bad_tri_recursive(DT, eid, &nbad, vid);
        }
    } else {
        find_bad_tri_recursive(DT, tri_s, &nbad, vid);
    }

    
    for (i = 0; i<nbad; i++){
        DT->bwork[i] = true;
    }

    // checking for duplicate elements
    for (i=0; i<nbad; i++){
        if ((*DT).bwork[i]){
            for (j=i+1; j<nbad; j++){
                if (DT->elems[DT->badtris[i]][0] == DT->elems[DT->badtris[j]][0] && DT->elems[DT->badtris[i]][1] == DT->elems[DT->badtris[j]][1] && DT->elems[DT->badtris[i]][2] == DT->elems[DT->badtris[j]][2]){
                    cout << "found duplicate elements, sibhfs is: " << endl;
                    cout <<  DT->sibhfs[DT->badtris[i]][0] << " " << DT->sibhfs[DT->badtris[i]][1] << " " << DT->sibhfs[DT->badtris[i]][2] << endl;
                    cout <<  DT->sibhfs[DT->badtris[j]][0] << " " << DT->sibhfs[DT->badtris[j]][1] << " " << DT->sibhfs[DT->badtris[j]][2] << endl;
                    DT->bwork[j] = false;
                }   
            }
        }
    }

    // adding facets to buffer
    int nsegs = 0;
    for (i=0; i<nbad; i++){
        for (j=0; j<3; j++){
            (*DT).facets[nsegs][0] = (*DT).elems[(*DT).badtris[i]][j];
            (*DT).facets[nsegs][1] = (*DT).elems[(*DT).badtris[i]][(j+1)%3];
            (*DT).vedge[nsegs] = (*DT).sibhfs[(*DT).badtris[i]][j];
            nsegs++;
        }
    }
    // find unique facets and reordering
    vector<int> order = facet_reorder(DT, &nsegs);

    for (i = 0; i<nsegs; i++){
        for (j=i+1; j<nsegs; j++){
            if (DT->facets[i][0] == DT->facets[i][1] && DT->facets[i][1] == DT->facets[i][1]){
                cout << "duplicate segments found: " << DT->facets[i][0] << "," << DT->facets[i][1] << endl;
            }
        }
    }

    int hfid;
    int ntri_b = (*DT).nelems;
    // Adding new triangles from Polygonal hole
    for (i=0; i<nsegs; i++){
        ps[0] = DT->coords[(*DT).facets[order[i]][0]];
        ps[1] = DT->coords[(*DT).facets[order[i]][1]];
        ps[2] = DT->coords[vid];
        if (abs(area_tri(ps)) > 1e-5) {
        (*DT).elems[(*DT).nelems] = {(*DT).facets[order[i]][0],(*DT).facets[order[i]][1],vid};
        (*DT).sibhfs[(*DT).nelems][0] = (*DT).vedge[order[i]];
        if ((*DT).vedge[order[i]] != 0){
            hfid = (*DT).vedge[order[i]];
            (*DT).sibhfs[hfid2eid(hfid)-1][hfid2lid(hfid)-1] = elids2hfid((*DT).nelems+1,1);
        } else {
            (*DT).on_boundary[(*DT).nelems] = true;
        }
        int eid1 = ntri_b + modi(i+1,nsegs) + 1;
        int eid2 = ntri_b + modi(i-1,nsegs) + 1;
        (*DT).sibhfs[(*DT).nelems][1] = elids2hfid(ntri_b + modi(i+1,nsegs) + 1, 3);
        (*DT).sibhfs[(*DT).nelems][2] = elids2hfid(ntri_b + modi(i-1,nsegs) + 1, 2);
        (*DT).nelems++;
        }
    }
}

/**
 * @brief Find triangles which have point in their circumcircle (optimal recursive algorithm)
 * 
 * @param DT Trangulation passed by reference
 * @param tri triangle to be focused on
 * @param nbad number of bad triangle passed by reference
 * @param vid point in question
 */
void find_bad_tri_recursive(Triangulation* DT, int tri, int *nbad, int vid){
    vector<vector<double>> xs = {{0,0},{0,0},{0,0}};
    vector<double> ps = (*DT).coords[vid];
    xs[0] = (*DT).coords[(*DT).elems[tri][0]];
    xs[1] = (*DT).coords[(*DT).elems[tri][1]];
    xs[2] = (*DT).coords[(*DT).elems[tri][2]];
    int eid;
    if (inside_circumtri(xs,ps) && !DT->delete_elem[tri]){
        (*DT).delete_elem[tri] = true;
        (*DT).badtris[*nbad] = tri;
        (*nbad)++;
        eid = hfid2eid((*DT).sibhfs[tri][0]);
        if (eid!=0 && !(*DT).delete_elem[eid-1]){
            find_bad_tri_recursive(DT,eid-1,nbad,vid);
        }
        eid = hfid2eid((*DT).sibhfs[tri][1]);
        if (eid!=0 && !(*DT).delete_elem[eid-1]){
            find_bad_tri_recursive(DT,eid-1,nbad,vid);
        }
        eid = hfid2eid((*DT).sibhfs[tri][2]);
        if (eid!=0 && !(*DT).delete_elem[eid-1]){
            find_bad_tri_recursive(DT,eid-1,nbad,vid);
        }
    } else {
        DT->delete_elem[tri] = false;
        return;
    }
}

/**
 * @brief Find triangles which have point in their circumcircle (optimal recursive algorithm) and exit if the point is outside domain
 * 
 * @param DT Trangulation passed by reference
 * @param tri triangle to be focused on
 * @param nbad number of bad triangle passed by reference
 * @param vid point in question
 * @param exitf Flag if the point lies outside the domain
 */
void find_bad_tri_recursive(Triangulation* DT, int tri, int *nbad, int vid, bool* exitf){
    if (!*exitf) {
        vector<vector<double>> xs = {{0,0},{0,0},{0,0}};
        double d;
        int j;
        vector<double> u(2);
        vector<double> v(2);
        vector<double> ps = (*DT).coords[vid];
        xs[0] = (*DT).coords[(*DT).elems[tri][0]];
        xs[1] = (*DT).coords[(*DT).elems[tri][1]];
        xs[2] = (*DT).coords[(*DT).elems[tri][2]];
        int eid,lid;
        if (inside_circumtri(xs,ps) && !DT->delete_elem[tri]){
            (*DT).delete_elem[tri] = true;
            (*DT).badtris[*nbad] = tri;
            (*nbad)++;

            for(j = 0; j<3; j++){
            eid = hfid2eid((*DT).sibhfs[tri][j]);
            if (eid!=0 && !(*DT).delete_elem[eid-1]){
                find_bad_tri_recursive(DT,eid-1,nbad,vid,exitf);
            } else if (!(*DT).delete_elem[eid-1]){
                u =  DT->coords[DT->elems[tri][(j+1)%3]] - DT->coords[DT->elems[tri][j]];
                v =  DT->coords[vid] - DT->coords[DT->elems[tri][j]];
                *exitf = u[0]*v[1] - u[1]*v[0] < 0;
                if (*exitf) {
                    *nbad = elids2hfid(tri+1,j+1);
                    (*DT).delete_elem[tri] = false;
                    return;
                } 
            }
            }
        } else {
            DT->delete_elem[tri] = false;
            return;
        } 
    } else {
        return;
    }
}

/**
 * @brief Find triangles which have point in their circumcircle (unoptimal slow algorithm)
 * 
 * @param DT Trangulation passed by reference
 * @param tri triangle to be focused on
 * @param nbad number of bad triangle passed by reference
 * @param vid point in question
 */
void find_bad_tri(Triangulation* DT, int tri, int *nbad, int vid){
    vector<vector<double>> xs = {{0,0},{0,0},{0,0}};
    vector<double> ps = (*DT).coords[vid];
    for (int i = 0; i<DT->nelems; i++){
        if (!DT->delete_elem[i]){
            xs[0] = (*DT).coords[(*DT).elems[i][0]];
            xs[1] = (*DT).coords[(*DT).elems[i][1]];
            xs[2] = (*DT).coords[(*DT).elems[i][2]];
            if (inside_circumtri(xs,ps)){
                (*DT).delete_elem[i] = true;
                (*DT).badtris[*nbad] = i;
                (*nbad)++;
            }
        }
    }
}

/**
 * @brief Delete triangles and reorganize data in Triangulation DT
 * 
 * @param DT Triangulation DT passed by reference
 */
void delete_tris(Triangulation* DT){
    int i,j;
    int nelems = 0;
    int sz = (*DT).sibhfs[0].size();
    vector<int> idx((*DT).nelems);
    vector<int> idx_rev((*DT).nelems);
    fill(idx_rev.begin(), idx_rev.end(),-1);

    // delete triangles to be deleted
    for (i = 0; i<(*DT).nelems; i++){
        if (!(*DT).delete_elem[i]){
            (*DT).elems[nelems] = (*DT).elems[i];
            (*DT).sibhfs[nelems] = (*DT).sibhfs[i];
            (*DT).delete_elem[nelems] = false;
            idx[nelems] = i;
            nelems++;
        } else {
            DT->delete_elem[i] = false;
        }
    }
    (*DT).nelems = nelems;

    for (i = 0; i<(*DT).nelems; i++){
        idx_rev[idx[i]] = i;
    }

    int nside;
    int hfid, eid, lid;
    for (i = 0; i<nelems; i++){
        nside = 0;
        for (j = 0; j < sz; j++){
            hfid = (*DT).sibhfs[i][j];
            if (!hfid == 0){
                eid = hfid2eid(hfid);
                lid = hfid2lid(hfid);
                if (!DT->delete_elem[eid-1]){
                    if (idx_rev[eid-1] == -1){
                        DT->sibhfs[i][j] = 0;
                    } else {
                        (*DT).sibhfs[i][j] = elids2hfid(idx_rev[eid-1]+1,lid);
                    }
                } else {
                    DT->sibhfs[i][j] = 0;
                }
                nside++;
            } else {
                DT->sibhfs[i][j] = 0;
            }
        }
        if (nside == sz){
            (*DT).on_boundary[i] = false;
        } else {
            (*DT).on_boundary[i] = true;
        }
        (*DT).delete_elem[i] = false;
    }
}
/**
 * @brief Delete triangles and reorganize data in Triangulation DT and keep track of specific triangle tri
 * 
 * @param DT Triangulation DT passed by reference
 * @param tri triangle pased by reference
 */
void delete_tris(Triangulation* DT, int* tri){
    int i,j;
    int nelems = 0;
    int sz = (*DT).sibhfs[0].size();
    vector<int> idx((*DT).nelems);
    vector<int> idx_rev((*DT).nelems);

    // delete triangles to be deleted
    for (i = 0; i<(*DT).nelems; i++){
        if (!(*DT).delete_elem[i]){
            (*DT).elems[nelems] = (*DT).elems[i];
            (*DT).sibhfs[nelems] = (*DT).sibhfs[i];
            (*DT).delete_elem[nelems] = false;
            idx[nelems] = i;
            nelems++;
        } else {
            DT->delete_elem[i] = false;
        }
    }
    (*DT).nelems = nelems;

    for (i = 0; i<(*DT).nelems; i++){
        idx_rev[idx[i]] = i;
    }

    int nside;
    int hfid, eid, lid;
    for (i = 0; i<nelems; i++){
        nside = 0;
        for (j = 0; j < sz; j++){
            hfid = (*DT).sibhfs[i][j];
            if (!hfid == 0){
                eid = hfid2eid(hfid);
                lid = hfid2lid(hfid);
                if (!DT->delete_elem[eid-1]){
                    if (idx_rev[eid-1] == -1){
                        DT->sibhfs[i][j] = 0;
                    } else {
                        (*DT).sibhfs[i][j] = elids2hfid(idx_rev[eid-1]+1,lid);
                    }
                } else {
                    DT->sibhfs[i][j] = 0;
                }
                nside++;
            } else {
                DT->sibhfs[i][j] = 0;
            }
        }
        if (nside == sz){
            (*DT).on_boundary[i] = false;
        } else {
            (*DT).on_boundary[i] = true;
        }
        (*DT).delete_elem[i] = false;
    }
    *tri = idx_rev[*tri];
}

// sub-functions necessary for delaunay
/// evalute alpha for delaunay refinement
double eval_alpha(const vector<vector<double>> xs,double r_ref){
    double a = norm(xs[1]-xs[0]);
    double b = norm(xs[2]-xs[1]);
    double c = norm(xs[2]-xs[0]);
    double s = 0.5*(a+b+c);
    double A = sqrt(s*(s-a)*(s-b)*(s-c));
    double r = a*b*c/(4*A);
    return r/r_ref;
}
/// find circumcenter for a triangle
static vector<double> circumcenter(const vector<vector<double>> xs){
    double ax = xs[0][0];
    double ay = xs[0][1];
    double bx = xs[1][0];
    double by = xs[1][1];
    double cx = xs[2][0];
    double cy = xs[2][1];
    double D = 2*(ax*(by-cy) + bx*(cy-ay) + cx*(ay-by));
    double ux = (ax*ax + ay*ay)*(by-cy) + \
        (bx*bx + by*by)*(cy-ay) + \
        (cx*cx + cy*cy)*(ay-by);
    double uy = (ax*ax + ay*ay)*(cx-bx) + \
        (bx*bx + by*by)*(ax-cx) + \
        (cx*cx + cy*cy)*(bx-ax);
    return {ux/D, uy/D};
}
/// find whether point is inside the circumcircle of a triangle
static bool inside_circumtri(const Mat xs, const vector<double> ps){
    vector<double> C = circumcenter(xs);
    double R = pow(xs[0][0]-C[0],2) + pow(xs[0][1] - C[1],2);
    bool D = pow(ps[0]-C[0],2) + pow(ps[1] - C[1],2) < R;
    return (D);
}
/// find whether point is inside a triangle
static bool inside_tri(const Mat &xs, const vector<double> &ps){
    double val1 = (ps[0]-xs[1][0])*(xs[0][1]-xs[1][1]) - (xs[0][0]-xs[1][0])*(ps[1]-xs[1][1]);
    double val2 = (ps[0]-xs[2][0])*(xs[1][1]-xs[2][1]) - (xs[1][0]-xs[2][0])*(ps[1]-xs[2][1]);
    double val3 = (ps[0]-xs[0][0])*(xs[2][1]-xs[0][1]) - (xs[2][0]-xs[0][0])*(ps[1]-xs[0][1]);
    bool has_neg = (val1 <= 0) || (val2<=0) || (val3<=0);
    bool has_pos = (val1 >= 0) || (val2>=0) || (val3>=0);
    return !(has_neg && has_pos);
}
/// Find unique facets and reorder them for bowyer watson algorithm
vector<int> facet_reorder(Triangulation* DT, int *nsegs){
    int nsegs2 = 0;
    int i,j;
    vector<vector<int>> segs2 = Zerosi(*nsegs,2);
    vector<int> vedge2(*nsegs);
    for (i=0;i<(*nsegs);i++){
        (*DT).bwork[i] = true;
    }
    for (i=0; i<(*nsegs); i++){
        if ((*DT).bwork[i]){
            for (j=0; j<(*nsegs); j++){
                if ((((*DT).facets[i][0] == (*DT).facets[j][0] && (*DT).facets[i][1] == (*DT).facets[j][1]) || \
                ((*DT).facets[i][0] == (*DT).facets[j][1] && (*DT).facets[i][1] == (*DT).facets[j][0])) && i!=j){
                    (*DT).bwork[i] = false;
                    (*DT).bwork[j] = false;
                }
            }
        }

        if ((*DT).bwork[i]){
            segs2[nsegs2][0] = (*DT).facets[i][0];
            segs2[nsegs2][1] = (*DT).facets[i][1];
            vedge2[nsegs2] = (*DT).vedge[i];
            nsegs2++;
        }
    }

    if (nsegs2 == 0){
        for (i=0; i<(*nsegs)/3; i++){
            cout << DT->badtris[i] << " ";
        }
        cout << endl;
        for (i=0; i<(*nsegs); i++){
            cout << (*DT).facets[i][0] << " " << (*DT).facets[i][1] << endl;
        }
    }

    *nsegs = nsegs2;
    for (i=0; i<(*nsegs); i++){
        (*DT).facets[i][0] = segs2[i][0];
        (*DT).facets[i][1] = segs2[i][1];
        (*DT).vedge[i] = vedge2[i];
        (*DT).bwork[i] = false;
    }
    vector<int> order(*nsegs);

    order[0] = 0;
    int vid;
    bool exitf;
    for (i=1; i<(*nsegs); i++){
        vid = (*DT).facets[order[i-1]][1];
        exitf = false;
        j=0;
        while ((j<(*nsegs)) && !exitf){
            if (vid == (*DT).facets[j][0]){
                order[i] = j;
                exitf = true;
            } else {
                j++;
            }
        }
    }
    return order;
}
/// find area of a triangle
static double area_tri(const Mat &xs){
    vector<double> u(2);
    vector<double> v(2);
    u =  xs[1]-xs[0];
    v =  xs[2]-xs[0];
    return abs(u[0]*v[1] - u[1]*v[0])/2;
}



// tiny private functions for array operations
/// find smallest x and y values in point set
static vector<double> min_array(const vector<vector<double>> &xs){
    int ndims = xs[0].size();
    vector<double> min_(ndims);

    for (int j = 0; j<ndims; j++){
        min_[j] = xs[0][j];
    }

    for (int i = 1; i<xs.size(); i++){
        for (int j = 0; j<ndims; j++){
            if (xs[i][j] < min_[j]){
                min_[j] = xs[i][j];
            }
        }
    }

    return min_;
}
/// find largest x and y values in point set
static vector<double> max_array(const vector<vector<double>> &xs){
    int ndims = xs[0].size();
    vector<double> min_(ndims);

    for (int j = 0; j<ndims; j++){
        min_[j] = xs[0][j];
    }

    for (int i = 1; i<xs.size(); i++){
        for (int j = 0; j<ndims; j++){
            if (xs[i][j] > min_[j]){
                min_[j] = xs[i][j];
            }
        }
    }

    return min_;
}
/// find center in a point set
static vector<double> find_center(const vector<vector<double>> &xs){
    int nv = xs.size();
    int ndims = 2;
    vector<double> center(2);
    for (int i=0;i<nv; i++){
        for (int j=0;j<ndims;j++){
            center[j] += xs[i][j];
        }
    }
    for (int j = 0; j<ndims; j++){
        center[j] = center[j]/((double) nv);
    }
    return center;
}

/// reorder pointset by angle
int myrandom (int i) { return std::rand()%i;}
void reorder(vector<vector<double>> &xs){
    srand ( unsigned ( time(0) ) );
    int nv = xs.size();
    int ndims = xs[0].size();
    assert(ndims == 2);
    vector<double> quantity(nv);
    random_shuffle(xs.begin(), xs.end(),myrandom);
    return;

    vector<double> center = find_center(xs);

    // measuring some quantity for each point
    vector<double> e = {1,0};
    vector<double> u(ndims);
    for (int n = 0; n<nv; n++){
        u[0] = xs[n][0] - center[0];
        u[1] = xs[n][1] - center[1];
        quantity[n] = acos(inner(u,e)/norm(u));
    }

    // basic sorting algo for order points based on some quantity
    int min_idx; 
    int j,n;
    for (n = 0; n<nv-1; n++){
        min_idx = n;
        for (j=n+1; j<nv; j++){
            if (quantity[j] < quantity[min_idx]){
                min_idx = j;
            }

            if (min_idx != n){
                swap(quantity[min_idx],quantity[n]);
                swap(xs[min_idx],xs[n]);
            }
        }
    }
}



// Sanity functions
/**
 * @brief Check whether triangulation has valid half-faces
 * 
 * @param DT Triangulation passed by reference
 * @return true 
 * @return false 
 */
bool check_sibhfs(Triangulation* DT){
    const vector<vector<int>> edges = {{0,1},{1,2},{2,0}};
    int nelems = DT->nelems;
    int hfid,eid,lid;
    cout << "nelems: " << nelems << endl;
    for (int i = 0; i<nelems; i++){
        for (int j = 0; j<3; j++){
            hfid = DT->sibhfs[i][j];
            if (hfid != 0){
            eid = hfid2eid(hfid);
            lid = hfid2lid(hfid);
            if (eid > nelems || lid > 3 || eid < 0 || DT->delete_elem[eid-1]){
                cout << "sibhfs is wrong at elem: " << i << " face: " << j << " oppeid: " << eid-1 << " opplid: " << lid-1 << endl;
            } 
            if (DT->elems[i][edges[j][0]] != DT->elems[eid-1][edges[lid-1][1]] || DT->elems[i][edges[j][1]] != DT->elems[eid-1][edges[lid-1][0]]){
                cout << "sides dont match is wrong at elem: " << i << " face: " << j << " oppeid: " << eid-1 << " opplid: " << lid-1 << " sibhfs: " << hfid << endl;
            }   
            }
        }
    }
} 
/**
 * @brief Find whether triangulation has valid Jacobian determinants for all triangles
 * 
 * @param DT Triangulation passed by reference
 * @return true 
 * @return false 
 */
bool check_jacobians(Triangulation* DT){
    const Mat dphi = {{-1,-1},{1,0},{0,1}};
    Mat ps = {{0,0},{0,0},{0,0}};
    Mat J;
    double detJ;
    bool check = true;
    for (int i = 0; i<DT->nelems; i++){
        ps[0] = DT->coords[DT->elems[i][0]];
        ps[1] = DT->coords[DT->elems[i][1]];
        ps[2] = DT->coords[DT->elems[i][2]];
        J = Transpose(ps)*dphi;
        detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
        if (detJ < 0){
            cout << "negative jacobian at eid: " << i << endl;
            check = false;
        }
    }
    return false;
}

// write mesh to files
void WrtieVtk_tri(const mesh &msh){
    FILE *fid;
    fid = fopen("test.vtk","w");
    fprintf(fid,"# vtk DataFile Version 3.0\n");
    fprintf(fid,"This file was written using writevtk_unstr.m\n");
    fprintf(fid,"ASCII\n");

    int nv = msh.coords.size();
    int ndims = 2;
    int nelems = msh.elems.size();

    // header for points
    fprintf(fid, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fid, "POINTS %i double",nv);

    // write out vertices
    if (ndims == 2){
        for (int i=0; i<nv;i++){
            fprintf(fid,"\n%g %g %g",msh.coords[i][0],msh.coords[i][1],0.0);
        }
    } else {
        for (int i=0; i<nv;i++){
            fprintf(fid,"\n%g %g %g",msh.coords[i][0],msh.coords[i][1],msh.coords[i][2]);
        }
    }

    // write out connectivity header
    fprintf(fid,"\n\nCELLS %i %i", nelems, 4*nelems);
    for (int i = 0; i<nelems; i++){
        fprintf(fid,"\n%d %d %d %d",3,msh.elems[i][0],msh.elems[i][1],msh.elems[i][2]);
    }

    // write out cell types
    fprintf(fid, "\n\nCELL_TYPES %i", nelems);
    for (int i = 0; i<nelems; i++){
        fprintf(fid,"\n%i",5);
    }

    fclose(fid);
}
void WrtieVtk_tri(const Triangulation &msh){
    FILE *fid;
    fid = fopen("test.vtk","w");
    fprintf(fid,"# vtk DataFile Version 3.0\n");
    fprintf(fid,"This file was written using writevtk_unstr.m\n");
    fprintf(fid,"ASCII\n");

    int nv = msh.coords.size();
    int ndims = msh.coords[0].size();
    int nelems = msh.nelems;
    int nelems_2=0;
    for (int i=0; i<nelems; i++){
        if (!msh.delete_elem[i]){
            nelems_2++;
        }
    }

    // header for points
    fprintf(fid, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fid, "POINTS %i double",nv);

    // write out vertices
    if (ndims == 2){
        for (int i=0; i<nv;i++){
            fprintf(fid,"\n%g %g %g",msh.coords[i][0],msh.coords[i][1],0.0);
        }
    } else {
        for (int i=0; i<nv;i++){
            fprintf(fid,"\n%g %g %g",msh.coords[i][0],msh.coords[i][1],msh.coords[i][2]);
        }
    }

    // write out connectivity header
    fprintf(fid,"\n\nCELLS %i %i", nelems_2, 4*nelems_2);
    for (int i = 0; i<nelems; i++){
        if (!msh.delete_elem[i]){
            fprintf(fid,"\n%d %d %d %d",3,msh.elems[i][0],msh.elems[i][1],msh.elems[i][2]);
        }
    }


    // write out cell types
    fprintf(fid, "\n\nCELL_TYPES %i", nelems_2);
    for (int i = 0; i<nelems; i++){
        if (!msh.delete_elem[i]){
            fprintf(fid,"\n%i",5);
        }
    }

    fclose(fid);
}

void WrtieVtk_tri(const Triangulation &msh, const vector<double> &data){
    FILE *fid;
    fid = fopen("test.vtk","w");
    fprintf(fid,"# vtk DataFile Version 3.0\n");
    fprintf(fid,"This file was written using writevtk_unstr.m\n");
    fprintf(fid,"ASCII\n");

    int nv = msh.coords.size();
    assert(nv == data.size());
    int ndims = msh.coords[0].size();
    int nelems = msh.elems.size();

    // header for points
    fprintf(fid, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fid, "POINTS %i double",nv);

    // write out vertices
    if (ndims == 2){
        for (int i=0; i<nv;i++){
            fprintf(fid,"\n%g %g %g",msh.coords[i][0],msh.coords[i][1],0.0);
        }
    } else {
        for (int i=0; i<nv;i++){
            fprintf(fid,"\n%g %g %g",msh.coords[i][0],msh.coords[i][1],msh.coords[i][2]);
        }
    }

    // write out connectivity header
    fprintf(fid,"\n\nCELLS %i %i", nelems, 4*nelems);
    for (int i = 0; i<nelems; i++){
        fprintf(fid,"\n%d %d %d %d",3,msh.elems[i][0],msh.elems[i][1],msh.elems[i][2]);
    }

    // write out cell types
    fprintf(fid, "\n\nCELL_TYPES %i", nelems);
    for (int i = 0; i<nelems; i++){
        fprintf(fid,"\n%i",5);
    }

    fprintf(fid, "\n\nPOINT_DATA %i",nv);
    fprintf(fid, "\nSCALARS meshdata double\n");
    fprintf(fid, "LOOKUP_TABLE default\n");
    for (int i = 0; i<nv; i++){
        fprintf(fid,"%g\n",data[i]);
    }

    fclose(fid);
}
void WriteObj_mesh(const mesh &msh){
    FILE *fid;
    fid = fopen("test.obj","w");

    int nv = msh.coords.size();
    int ndims = msh.coords[0].size();
    int nelems = msh.elems.size();

    if (ndims == 2){
        for (int i=0; i<nv;i++){
            fprintf(fid,"\nv %f %f %f",msh.coords[i][0],msh.coords[i][1],0.0);
        }
    } else {
        for (int i=0; i<nv;i++){
            fprintf(fid,"\n%g %g %g",msh.coords[i][0],msh.coords[i][1],msh.coords[i][2]);
        }
    }

    for (int i = 0; i<nelems; i++){
        fprintf(fid,"\nf %d %d %d",msh.elems[i][0]+1,msh.elems[i][1]+1,msh.elems[i][2]+1);
    }

    fclose(fid);
}