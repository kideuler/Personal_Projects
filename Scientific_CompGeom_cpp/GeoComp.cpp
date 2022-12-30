#include "GeoComp.hpp"
#include <unistd.h>
using namespace std;

int hfid2eid(int hfid){
    return (hfid >> 8);
}
int hfid2lid(int hfid){
    return (hfid&255) + 1;
}
int elids2hfid(int eid, int lid){
    return ((eid << 8) + lid - 1);
}

// Functions for element linked list
void push_elem(Elems** head, vector<int> elem){
    Elems* new_Elem = new Elems;
    int sz = elem.size();
    new_Elem->elem.resize(sz);
    for (int i = 0; i<sz; i++){
        new_Elem->elem[i] = elem[i];
    }
    new_Elem->next = (*head);
    new_Elem->prev = NULL;
    if ((*head) != NULL){
        (*head)->prev = new_Elem;
    }
    (*head) = new_Elem;
}
Elems* pop_elem(Elems* head){
    if (head == NULL){
        return NULL;
    }
    Elems* temp = head;
    head = head->next;
    if (head != NULL){
        head->prev = NULL;
    }
    delete temp;
    return head;
}
void insert_elem(Elems* prev_elem, vector<int> elem){
    if (prev_elem == NULL){
        push_elem(NULL, elem);
    } else {
        Elems* new_Elem = new Elems;
        int sz = elem.size();
        new_Elem->elem.resize(sz);
        for (int i = 0; i<sz; i++){
            new_Elem->elem[i] = elem[i];
        }
        new_Elem->next = prev_elem->next;
        new_Elem->prev = prev_elem;
        prev_elem->next = new_Elem;
        new_Elem->next->prev = new_Elem;
    }
}
void delete_elem(Elems** head, Elems** elem){
    if ((*head) == NULL || *elem == NULL){
        return;
    }

    if (*head == *elem){
        *head = (*elem)->next;
    }

    if ((*elem)->next != NULL){
        (*elem)->next->prev = (*elem)->prev;
    }

    if ((*elem)->prev != NULL){
        (*elem)->prev->next = (*elem)->next;
    }

    Elems* temp = *elem;
    *elem = (*elem)->next;
    delete temp;
}

void print_LLElems(Elems* head){
    Elems* curr = head;
    while (curr != NULL){
        for (int i = 0; i<curr->elem.size(); i++){
            cout << curr->elem[i] << " ";
        }
        cout << endl;
        curr = curr->next;
    }
}

// subfunctions for delaunay triangulation
double eval_alpha(const vector<vector<double>> xs,double r_ref);
static vector<double> circumcenter(const vector<vector<double>> xs);
static bool inside_tri(const Mat &xs, const vector<double> &ps);
static bool inside_circumtri(const Mat xs, const vector<double> ps);
static double area_tri(const Mat &xs);
void find_bad_tri_recursive(Triangulation* DT, int tri, int *nbad, int vid);
void find_bad_tri_recursive(Triangulation* DT, int tri, int *nbad, int vid, bool* exitf);
vector<int> facet_reorder(Triangulation *DT, int *nsegs);


// array functions for xs
static vector<double> min_array(const vector<vector<double>> &xs);
static vector<double> max_array(const vector<vector<double>> &xs);
static vector<double> find_center(const vector<vector<double>> &xs);
void reorder(vector<vector<double>> &xs);

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
void GeoComp_refine(Triangulation* DT, double r_ref){
    auto f = [r_ref](vector<double> xs) {return r_ref; };
    GeoComp_refine(DT, f);
}
void GeoComp_refine(Triangulation* DT, function<double(vector<double>)> r_ref){
    int n,i;
    bool exitl;
    double alpha;
    int nv = (*DT).coords.size();
    Mat ps = Zeros(3,2);
    int tri;
    (*DT).coords.resize((*DT).elems.size());
    vector<double> C;
    n=0;
    while (n<(*DT).nelems){
        if (!(*DT).delete_elem[n]){
            ps[0] = (*DT).coords[(*DT).elems[n][0]];
            ps[1] = (*DT).coords[(*DT).elems[n][1]];
            ps[2] = (*DT).coords[(*DT).elems[n][2]];
            C = circumcenter(ps);
            if (abs(C[0] > 1e6 || abs(C[1]) > 1e6)){
                cout << "abnormally large point added from points " << C[0] << "," << C[1] << endl;
                printMat(ps); 
            }
            alpha = eval_alpha(ps,r_ref(C));
            if (alpha > 1.1){
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

Triangulation GeoComp_Delaunay_Triangulation(vector<vector<double>> &xs){
    // size checking
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
    reorder(xs);

    for (n=0; n<nv; n++){
        DT.coords[n] = xs[n];
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
    cout << "created " << DT.nelems << " triangles from initial points" << endl;
    return DT;
}

static int modi(int a, int b){
    int z = a - b*(a/b);
    if (z<0){
        z=z+b;
    }
    return z;
}
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
    
    if (nbad == 0){
        cout << "error occured no triangles found recursively, starting tri: " << tri_s << endl;
         
    }

    int nsegs = 0;
    for (i=0; i<nbad; i++){
        for (j=0; j<3; j++){
            (*DT).facets[nsegs][0] = (*DT).elems[(*DT).badtris[i]][j];
            (*DT).facets[nsegs][1] = (*DT).elems[(*DT).badtris[i]][(j+1)%3];
            (*DT).vedge[nsegs] = (*DT).sibhfs[(*DT).badtris[i]][j];
            nsegs++;
        }
    }
    vector<int> order = facet_reorder(DT, &nsegs);

    int hfid;
    int ntri_b = (*DT).nelems;

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

void find_bad_tri_recursive(Triangulation* DT, int tri, int *nbad, int vid){
    vector<vector<double>> xs = {{0,0},{0,0},{0,0}};
    vector<double> ps = (*DT).coords[vid];
    xs[0] = (*DT).coords[(*DT).elems[tri][0]];
    xs[1] = (*DT).coords[(*DT).elems[tri][1]];
    xs[2] = (*DT).coords[(*DT).elems[tri][2]];
    int eid;
    if (inside_circumtri(xs,ps)){
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
        if (inside_circumtri(xs,ps)){
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
/*
void find_bad_tri(Triangulation* DT, int tri, int *nbad, int vid, bool* exit){
    for (int i = 0; i<)
}
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
double eval_alpha(const vector<vector<double>> xs,double r_ref){
    double a = norm(xs[1]-xs[0]);
    double b = norm(xs[2]-xs[1]);
    double c = norm(xs[2]-xs[0]);
    double s = 0.5*(a+b+c);
    double A = sqrt(s*(s-a)*(s-b)*(s-c));
    double r = a*b*c/(4*A);
    return r/r_ref;
}
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
static bool inside_circumtri(const Mat xs, const vector<double> ps){
    vector<double> C = circumcenter(xs);
    double R = norm(xs[0]-C);
    return (norm(ps-C) <= R);
}
static bool inside_tri(const Mat &xs, const vector<double> &ps){
    double val1 = (ps[0]-xs[1][0])*(xs[0][1]-xs[1][1]) - (xs[0][0]-xs[1][0])*(ps[1]-xs[1][1]);
    double val2 = (ps[0]-xs[2][0])*(xs[1][1]-xs[2][1]) - (xs[1][0]-xs[2][0])*(ps[1]-xs[2][1]);
    double val3 = (ps[0]-xs[0][0])*(xs[2][1]-xs[0][1]) - (xs[2][0]-xs[0][0])*(ps[1]-xs[0][1]);
    bool has_neg = (val1 < 0) || (val2<0) || (val3<0);
    bool has_pos = (val1 > 0) || (val2>0) || (val3>0);
    return !(has_neg && has_pos);
}
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
            for (j=i+1; j<(*nsegs); j++){
                if (((*DT).facets[i][0] == (*DT).facets[j][0] && (*DT).facets[i][1] == (*DT).facets[j][1]) || \
                ((*DT).facets[i][0] == (*DT).facets[j][1] && (*DT).facets[i][1] == (*DT).facets[j][0])){
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
        for (i=0; i<(*nsegs); i++){
            cout << (*DT).facets[i][0] << " " << (*DT).facets[i][1] << endl;
            vector<int> order(nsegs2);
            return order;
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
static double area_tri(const Mat &xs){
    vector<double> u(2);
    vector<double> v(2);
    u =  xs[1]-xs[0];
    v =  xs[2]-xs[0];
    return abs(u[0]*v[1] - u[1]*v[0])/2;
}



// tiny private functions for array operations
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
void reorder(vector<vector<double>> &xs){
    int nv = xs.size();
    int ndims = xs[0].size();
    assert(ndims == 2);
    vector<double> quantity(nv);

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