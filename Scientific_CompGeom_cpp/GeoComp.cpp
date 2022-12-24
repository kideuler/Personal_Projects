#include "GeoComp.hpp"
using namespace std;

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

static vector<double> min_array(const vector<vector<double>> &xs);
static vector<double> max_array(const vector<vector<double>> &xs);

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

mesh GeoComp_DT2mesh(Triangulation &DT){
    mesh msh;
    int i,j,sz,k,nv,ndims;
    nv = DT.coords.size();
    ndims = DT.coords[0].size();
    msh.coords = Zeros(nv,ndims);
    for (i=0; i<nv; i++){
        for (j=0; j<ndims; j++){
            msh.coords[i][j] = DT.coords[i][j];
        }
    }
    msh.elems.resize(DT.nelems);
    Elems* curr = DT.elem_head;
    k=0;
    while (curr != NULL){
        sz = curr->elem.size();
        msh.elems[k].resize(sz);
        for (i = 0; i<sz;i++){
            msh.elems[k][i] = curr->elem[i];
        }
        k++;
        curr = curr->next;
    }
    return msh;
}

Triangulation GeoComp_Delaunay_Triangulation(vector<vector<double>> &xs){
    // size checking
    int nv = xs.size();
    int n;

    Triangulation DT;
    DT.coords = Zeros(nv+3,2);
    DT.elem_head = NULL;
    DT.sibhfs_head = NULL;
    DT.facets = Zerosi(nv*nv,2);
    DT.bwork.resize(nv);
    vector<double> a = min_array(xs);
    vector<double> b = max_array(xs);

    // reorder points for optimal triangle search O(N*log(N))
    reorder(xs);

    for (n=0; n<nv; n++){
        DT.coords[n][0] = xs[n][0];
        DT.coords[n][1] = xs[n][1];
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

    push_elem(&DT.elem_head,{nv,nv+1,nv+2});
    push_elem(&DT.sibhfs_head,{0,0,0});

    // loop inserting each point into triangulation
    DT.nelems = 1;
    printMat(DT.coords);
    for (int n = 0; n<nv; n++){
        // inserting node into the triangulation using Bowyer-Watson algorithm
        Bowyer_watson_insert_point2d(DT,n);
        print_LLElems(DT.elem_head);
        cout << endl;
    }

    return DT;
}

bool inside_circumtri(const Triangulation &DT, const Elems* tri, int vid){
    vector<vector<double>> J = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
    vector<double> u(2);
    cout << tri->elem[0] << " " << tri->elem[1] << " " << tri->elem[2] << " incirc" << endl;
    for (int i = 0; i<3; i++){
        for (int j = 0; j<2; j++){
            J[i][j] = DT.coords[tri->elem[i]][j] - DT.coords[vid][j];
            u[j] = DT.coords[tri->elem[i]][j] - DT.coords[vid][j]; 
        }
        J[i][2] = inner(u,u);
    }

    double D = J[0][2]*(J[1][0]*J[2][1] - J[2][0]*J[1][1]) + \
        J[1][2]*(J[2][0]*J[0][1] - J[0][0]*J[2][1]) + \
        J[2][2]*(J[0][0]*J[1][1] - J[1][0]*J[0][1]);
    return (D > 0);
}

vector<int> facet_reorder(Triangulation &DT, int *nsegs){
    int nsegs2 = 0;
    int i,j;
    vector<vector<int>> segs2 = Zerosi(*nsegs,2);
    for (i=0;i<(*nsegs);i++){
        DT.bwork[i] = true;
    }
    for (i=0; i<(*nsegs);i++){
        if (DT.bwork[i]){
            for (j=i+1;j<(*nsegs);j++){
                if (((DT.facets[i][0] == DT.facets[j][0]) && \
                (DT.facets[i][1] == DT.facets[j][1])) || \
                ((DT.facets[i][0] == DT.facets[j][1]) && \
                (DT.facets[i][1] == DT.facets[j][0]))){
                    DT.bwork[i] = false;
                    DT.bwork[j] = false;
                }
            }
        }

        if (DT.bwork[i]){
            segs2[nsegs2][0] = DT.facets[i][0];
            segs2[nsegs2][1] = DT.facets[i][1];
            nsegs2++;
        }
    }

    *nsegs = nsegs2;
    for (i=0; i<(*nsegs); i++){
        DT.facets[i][0] = segs2[i][0];
        DT.facets[i][1] = segs2[i][1];
        DT.bwork[i] = false;
    }

    vector<int> order(*nsegs);
    order[0] = 0;
    int vid;
    bool exitf;
    for (i=1;i<(*nsegs);i++){
        vid = DT.facets[order[i-1]][1];
        exitf = false;
        j=0;
        while ((j<(*nsegs)) && !exitf){
            if (vid == DT.facets[j][0]){
                order[i] = j;
                exitf = true;
            } else {
                j++;
            }
        }
    }

    return order;
}


void Bowyer_watson_insert_point2d(Triangulation &DT, int vid){
    static vector<vector<int>> edges = {{0,1},{1,2},{2,0}};
    int i,j;
    // find all triangles whos point is inside circumcircle 
    Elems* curr = DT.elem_head;
    int nsegs = 0;
    while (curr != NULL){
        if (inside_circumtri(DT, curr, vid)){
            cout << "current elem: " << curr->elem[0] << " " << curr->elem[1] << " " << curr->elem[2] << endl;
            for (i=0; i<3; i++){
                DT.facets[nsegs][0] = curr->elem[edges[i][0]];
                DT.facets[nsegs][1] = curr->elem[edges[i][1]];
                nsegs++;
            }
            cout << curr << endl;
            delete_elem(&DT.elem_head, &curr);
            DT.nelems--;
            cout << curr << endl;
        } else {
            curr = curr->next;
        }
    }
    // reorder facets
    cout << "nsegs: " << nsegs << endl;
    vector<int> order = facet_reorder(DT,&nsegs);

    for (i=0; i<nsegs; i++){
        cout << DT.facets[order[i]][0] << " " << DT.facets[order[i]][1] << " nsegs: " << nsegs << endl;
        push_elem(&DT.elem_head,{DT.facets[order[i]][0], DT.facets[order[i]][1], vid});
        DT.nelems++;
    }
}


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


// tiny private functions
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
/*
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
*/