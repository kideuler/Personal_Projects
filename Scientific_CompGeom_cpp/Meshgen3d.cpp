#include "GeoComp.hpp"
#include <unistd.h>
using namespace std;


bool find_enclosing_tet(Mesh* DT, int* tet, int vid);


Mesh GeoComp_Delaunay_Mesh3d(vector<vector<double>> &xs){
    // size checking
    default_random_engine re;
    uniform_real_distribution<double> unif(-1, 1);
    int nv = xs.size();
    int n,i,j;

    Mesh DT;
    int ub = 2*nv*nv;
    DT.coords = Zeros(nv+3,3);
    DT.elems = Zerosi(ub,4);
    DT.sibhfs = Zerosi(ub,4);
    DT.facets = Zerosi(ub,3);
    DT.delete_elem.resize(ub);
    DT.on_boundary.resize(ub);
    vector<double> a = min_array(xs);
    vector<double> b = max_array(xs);
    double dx = (b[0]-a[0])/10;
    double dy = (b[1]-a[1])/10;
    double dz = (b[2]-a[2])/10;
    double d = max(max(dx,dy),dz);

    for (n=0; n<nv; n++){
        DT.coords[n][0] = (xs[n][0]-a[0])/d + 1e-4*unif(re);
        DT.coords[n][1] = (xs[n][1]-a[1])/d + 1e-4*unif(re);
        DT.coords[n][2] = (xs[n][2]-a[2])/d + 1e-4*unif(re);
    }



    DT.coords[nv][0] = -100;
    DT.coords[nv][1] = -100;
    DT.coords[nv][2] = -100;
    DT.coords[nv+1][0] = 100;
    DT.coords[nv+1][1] = -100;
    DT.coords[nv+1][2] = -100;
    DT.coords[nv+2][0] = -100;
    DT.coords[nv+2][1] = 100;
    DT.coords[nv+2][2] = -100;
    DT.coords[nv+3][0] = 0;
    DT.coords[nv+3][1] = 0;
    DT.coords[nv+3][2] = 100;

    DT.elems[0] = {nv,nv+1,nv+2,nv+3};
    DT.sibhfs[0] = {0,0,0,0};

    // loop inserting each point into Mesh
    DT.nelems = 1;
    int tet = -1;
    bool exitl,inside;
    int vid;
    for (int n = 0; n < nv; n++){
        vid = n;
        tet = DT.nelems-1;
        inside = find_enclosing_tet(&DT, &tet, vid);
        if (!inside){
            cout << "no enclosing tet found" << endl;
        }
        // inserting node into the Mesh using Bowyer-Watson algorithm
        //Flip_Insertion_tet(&DT,&vid,tet);
        if ((double) DT.nelems >= (0.8)*((double) ub)){
            cout << "approaching size bound, freeing up space" << endl;
            //delete_tets(&DT);
        }
    }

    for (int n = 0; n<DT.nelems; n++){
        for (int i = 0; i<4; i++){
            if (DT.elems[n][i] >= nv){
                DT.delete_elem[n] = true;
            }
        }
    }

    //delete_tets(&DT);
    DT.coords.resize(nv);
    for (n=0; n<nv; n++){
        DT.coords[n][0] = xs[n][0];
        DT.coords[n][1] = xs[n][1];
        DT.coords[n][2] = xs[n][2];
    }
    cout << "created " << DT.nelems << " tetrahedra from initial points" << endl;
    return DT;
}




vector<double> normal3d(const vector<vector<double>> &xs, int facet){
    vector<vector<int>> Faces = {{0,2,1},{0,1,3},{1,2,3},{2,0,3}};
    vector<double> a = {0,0,0};
    vector<double> b = {0,0,0};
    a[0] = xs[Faces[facet][1]][0] - xs[Faces[facet][0]][0];
    a[1] = xs[Faces[facet][1]][1] - xs[Faces[facet][0]][1];
    a[2] = xs[Faces[facet][1]][2] - xs[Faces[facet][0]][2];
    b[0] = xs[Faces[facet][2]][0] - xs[Faces[facet][0]][0];
    b[1] = xs[Faces[facet][2]][1] - xs[Faces[facet][0]][1];
    b[2] = xs[Faces[facet][2]][2] - xs[Faces[facet][0]][2];

    vector<double> N = {0,0,0};
    N[0] = a[1]*b[2] - a[2]*b[1];
    N[1] = a[2]*b[0] - a[0]*b[2];
    N[2] = a[0]*b[1] - a[1]*b[0];
    return N;
}

bool inside_tet(const Mat &xs, const vector<double> &ps){
    vector<vector<int>> Faces = {{0,2,1},{0,1,3},{1,2,3},{2,0,3}};

    vector<double> N = normal3d(xs,0);
    vector<double> a = {xs[3][0]-xs[Faces[0][0]][0], xs[3][1]-xs[Faces[0][0]][1], xs[3][2]-xs[Faces[0][0]][2]};
    double D0 = N[0]*a[0]+N[1]*a[1]+N[2]*a[2];
    vector<double> p = {ps[0]-xs[Faces[0][0]][0], ps[1]-xs[Faces[0][0]][1], ps[2]-xs[Faces[0][0]][2]};
    double Dp = N[0]*p[0]+N[1]*p[1]+N[2]*p[2];
    bool S1 = signbit(D0) == signbit(Dp);

    N = normal3d(xs,1);
    a = {xs[2][0]-xs[Faces[1][0]][0], xs[2][1]-xs[Faces[1][0]][1], xs[2][2]-xs[Faces[1][0]][2]};
    D0 = N[0]*a[0]+N[1]*a[1]+N[2]*a[2];
    p = {ps[0]-xs[Faces[1][0]][0], ps[1]-xs[Faces[1][0]][1], ps[2]-xs[Faces[1][0]][2]};
    Dp = N[0]*p[0]+N[1]*p[1]+N[2]*p[2];
    bool S2 = signbit(D0) == signbit(Dp);

    N = normal3d(xs,2);
    a = {xs[0][0]-xs[Faces[2][0]][0], xs[0][1]-xs[Faces[2][0]][1], xs[0][2]-xs[Faces[2][0]][2]};
    D0 = N[0]*a[0]+N[1]*a[1]+N[2]*a[2];
    p = {ps[0]-xs[Faces[2][0]][0], ps[1]-xs[Faces[2][0]][1], ps[2]-xs[Faces[2][0]][2]};
    Dp = N[0]*p[0]+N[1]*p[1]+N[2]*p[2];
    bool S3 = signbit(D0) == signbit(Dp);

    N = normal3d(xs,3);
    a = {xs[1][0]-xs[Faces[3][0]][0], xs[1][1]-xs[Faces[3][0]][1], xs[1][2]-xs[Faces[3][0]][2]};
    D0 = N[0]*a[0]+N[1]*a[1]+N[2]*a[2];
    p = {ps[0]-xs[Faces[3][0]][0], ps[1]-xs[Faces[3][0]][1], ps[2]-xs[Faces[3][0]][2]};
    Dp = N[0]*p[0]+N[1]*p[1]+N[2]*p[2];
    bool S4 = signbit(D0) == signbit(Dp);

    return S1 && S2 && S3 && S4;
}

bool find_enclosing_tet(Mesh* DT, int* tet, int vid){
    int v1,v2,v3,v4,i,hfid;
    bool stop;
    int iters = 0;
    i = 0;
    stop = false;
    vector<vector<double>> xs = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
    while (!stop){
        v1 = DT->elems[*tet][0];
        v2 = DT->elems[*tet][1];
        v3 = DT->elems[*tet][2];
        v4 = DT->elems[*tet][3];
        if (DT->delete_elem[*tet]){
            cout << "deleted tet passed ran through in find_enclosing_tet, should not be" << endl;

        }
        xs[0] = DT->coords[v1];
        xs[1] = DT->coords[v2];
        xs[2] = DT->coords[v3];
        xs[3] = DT->coords[v4];
        if (inside_tet(xs,DT->coords[vid])){
            stop = true;
            return true;
        } else {
            double AP[3] = {DT->coords[vid][0]-DT->coords[v1][0], DT->coords[vid][1]-DT->coords[v1][1], DT->coords[vid][2]-DT->coords[v1][2]};
            double BP[3] = {DT->coords[vid][0]-DT->coords[v2][0], DT->coords[vid][1]-DT->coords[v2][1], DT->coords[vid][2]-DT->coords[v2][2]};
            double CP[3] = {DT->coords[vid][0]-DT->coords[v3][0], DT->coords[vid][1]-DT->coords[v3][1], DT->coords[vid][2]-DT->coords[v3][2]};
            double DP[3] = {DT->coords[vid][0]-DT->coords[v4][0], DT->coords[vid][1]-DT->coords[v4][1], DT->coords[vid][2]-DT->coords[v4][2]};
            vector<double> N1 = normal3d(xs,0);
            vector<double> N2 = normal3d(xs,1);
            vector<double> N3 = normal3d(xs,2);
            vector<double> N4 = normal3d(xs,3);
            double S1 = AP[0]*N1[0]+AP[1]*N1[1]+AP[2]*N1[2];
            double S2 = BP[0]*N2[0]+BP[1]*N2[1]+BP[2]*N2[2];
            double S3 = CP[0]*N3[0]+CP[1]*N3[1]+CP[2]*N3[2];
            double S4 = DP[0]*N4[0]+DP[1]*N4[1]+DP[2]*N4[2];
            if ((S1>0)&&(S1>=S2)&&(S1>=S3)&&(S1>=S4)){
                hfid = DT->sibhfs[*tet][0];
                if (hfid != 0){
                    *tet = hfid2eid(hfid)-1;
                } else {
                    stop = true;
                    *tet = elids2hfid(*tet+1,1);
                    return false;
                }
            } else if ((S2>0)&&(S2>=S1)&&(S2>=S3)&&(S2>=S4)) {
                hfid = DT->sibhfs[*tet][1];
                if (hfid != 0){
                    *tet = hfid2eid(hfid)-1;
                } else {
                    stop = true;
                    *tet =  elids2hfid(*tet+1,2);
                    return false;
                }
            } else if ((S3>0)&&(S3>=S1)&&(S3>=S2)&&(S3>=S4)){
                hfid = DT->sibhfs[*tet][2];
                if (hfid != 0){
                    *tet = hfid2eid(hfid)-1;
                } else {
                    stop = true;
                    *tet = elids2hfid(*tet+1,3);
                    return false;
                }
            } else if ((S4>0)&&(S4>=S1)&&(S4>=S2)&&(S4>=S3)){
                hfid = DT->sibhfs[*tet][3];
                if (hfid != 0){
                    *tet = hfid2eid(hfid)-1;
                } else {
                    stop = true;
                    *tet = elids2hfid(*tet+1,4);
                    return false;
                }
            } 
            
        }
        iters++;
    }
    *tet = -1;
    return false;
    
}


void WrtieVtk_tet(const Mesh &msh){
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
            fprintf(fid,"\n%g %g %g",msh.coords[i][0],msh.coords[i][1],msh.coords[i][2]);
        }
    } else {
        for (int i=0; i<nv;i++){
            fprintf(fid,"\n%g %g %g",msh.coords[i][0],msh.coords[i][1],msh.coords[i][2]);
        }
    }

    // write out connectivity header
    fprintf(fid,"\n\nCELLS %i %i", nelems_2, 5*nelems_2);
    for (int i = 0; i<nelems; i++){
        if (!msh.delete_elem[i]){
            fprintf(fid,"\n%d %d %d %d %d",4,msh.elems[i][0],msh.elems[i][1],msh.elems[i][2],msh.elems[i][3]);
        }
    }


    // write out cell types
    fprintf(fid, "\n\nCELL_TYPES %i", nelems_2);
    for (int i = 0; i<nelems; i++){
        if (!msh.delete_elem[i]){
            fprintf(fid,"\n%i",10);
        }
    }

    fclose(fid);
}