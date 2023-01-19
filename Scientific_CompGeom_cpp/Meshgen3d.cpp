#include "GeoComp.hpp"
#include <unistd.h>
using namespace std;


bool find_enclosing_tet(Mesh* DT, int* tet, int vid);
int inside_tet(const Mat &xs, const vector<double> &ps);
void Flip_Insertion_tet(Mesh* DT, int vid, int tet_s);
bool check_jacobians_tet(Mesh* DT);
bool check_jac(Mat xs);
bool check_sibhfs_tet(Mesh* DT);


Mesh GeoComp_Delaunay_Mesh3d(vector<vector<double>> &xs){
    // size checking
    default_random_engine re;
    uniform_real_distribution<double> unif(-1, 1);
    int nv = xs.size();
    int n,i,j;

    Mesh DT;
    int ub = max(2*nv*nv,40);
    DT.coords = Zeros(nv+4,3);
    DT.elems = Zerosi(ub,4);
    DT.sibhfs = Zerosi(ub,4);
    DT.facets = Zerosi(ub,3);
    DT.delete_elem.resize(ub);
    DT.on_boundary.resize(ub);
    vector<double> a = min_array(xs);
    vector<double> b = max_array(xs);
    double dx = (b[0]-a[0]);
    double dy = (b[1]-a[1]);
    double dz = (b[2]-a[2]);
    double d = max(max(dx,dy),dz);

    for (n=0; n<nv; n++){
        DT.coords[n][0] = (xs[n][0]-a[0])/d + 1e-4*unif(re);
        DT.coords[n][1] = (xs[n][1]-a[1])/d + 1e-4*unif(re);
        DT.coords[n][2] = (xs[n][2]-a[2])/d + 1e-4*unif(re);
        
        DT.coords[n][0] = xs[n][0];
        DT.coords[n][1] = xs[n][1];
        DT.coords[n][2] = xs[n][2];
    }
    printMat(DT.coords);



    DT.coords[nv][0] = -10;
    DT.coords[nv][1] = -10;
    DT.coords[nv][2] = -10;
    DT.coords[nv+1][0] = 10;
    DT.coords[nv+1][1] = -10;
    DT.coords[nv+1][2] = -10;
    DT.coords[nv+2][0] = 0;
    DT.coords[nv+2][1] = 10;
    DT.coords[nv+2][2] = -10;
    DT.coords[nv+3][0] = 0;
    DT.coords[nv+3][1] = 0;
    DT.coords[nv+3][2] = 10;

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
        /*
        bool stop = false;
        vector<vector<double>> xs = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
        i = tet;
        while (i>-1 && !stop){
            if(!DT.delete_elem[i]){
                xs[0] = DT.coords[DT.elems[i][0]];
                xs[1] = DT.coords[DT.elems[i][1]];
                xs[2] = DT.coords[DT.elems[i][2]];
                xs[3] = DT.coords[DT.elems[i][3]];
                if (inside_tet(xs, DT.coords[vid])==0){
                    tet = i;
                    stop = true;
                    inside = true;
                }
            }
            i--;
        }
        */
        if (!inside){
            cout << "no enclosing tet found for starting tet " << tet << endl;
        }
        // inserting node into the Mesh using Bowyer-Watson algorithm
        Flip_Insertion_tet(&DT,vid,tet);
        if ((double) DT.nelems >= (0.8)*((double) ub)){
            cout << "approaching size bound, freeing up space" << endl;
            //delete_tets(&DT);
        }
    }

    for (int n = 0; n<DT.nelems; n++){
        for (int i = 0; i<4; i++){
            if (DT.elems[n][i] >= nv){
                //DT.delete_elem[n] = true;
            }
        }
    }

    //delete_tets(&DT);
    //DT.coords.resize(nv);
    for (n=0; n<nv; n++){
        DT.coords[n][0] = xs[n][0];
        DT.coords[n][1] = xs[n][1];
        DT.coords[n][2] = xs[n][2];
    }
    cout << "created " << DT.nelems << " tetrahedra from initial points" << endl;
    bool check = check_jacobians_tet(&DT);
    check = check_sibhfs_tet(&DT);
    return DT;
}

void Flip_Insertion_tet(Mesh* DT, int vid, int tet_s){
    DT->delete_elem[tet_s] = true;
    vector<vector<int>> Faces = {{1,0,2},{0,1,3},{1,2,3},{2,0,3}};
    int hfid,eid,lid;

    vector<int> tet = {DT->elems[tet_s][0], DT->elems[tet_s][1], DT->elems[tet_s][2], DT->elems[tet_s][3]};
    vector<int> sib = {DT->sibhfs[tet_s][0], DT->sibhfs[tet_s][1], DT->sibhfs[tet_s][2], DT->sibhfs[tet_s][3]};
    
    // splitting tetangle and adding subtetangles to the stack if they are not on the boundary
    stack* head = NULL;
    vector<int> eids = {DT->nelems, DT->nelems+1, DT->nelems+2, DT->nelems+3};
    vector<vector<double>> xs = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};

    // tet 1
    DT->elems[eids[0]] = {vid, tet[Faces[0][0]], tet[Faces[0][1]],tet[Faces[0][2]]};
    hfid = sib[0];
    DT->sibhfs[eids[0]] = {elids2hfid(eids[1]+1 ,1), elids2hfid(eids[2]+1 ,1), hfid, elids2hfid(eids[3]+1,1)};
    if (hfid2eid(hfid) > 0){
        DT->sibhfs[hfid2eid(hfid)-1][hfid2lid(hfid)-1] = elids2hfid(eids[0]+1, 3);
        push_stack(&head, elids2hfid(eids[0]+1, 3));
        DT->on_boundary[eids[0]] = false;
    } else {
        DT->on_boundary[eids[0]] = true;
    }
    xs[0] = DT->coords[DT->elems[eids[0]][0]];
    xs[1] = DT->coords[DT->elems[eids[0]][1]];
    xs[2] = DT->coords[DT->elems[eids[0]][2]];
    xs[3] = DT->coords[DT->elems[eids[0]][3]];
    if (!check_jac(xs)){
        cout << "negative jacobian at " << eids[0] << endl;
    }
    cout << "elem " << eids[0] << " oppelems: " << hfid2eid(DT->sibhfs[eids[0]][0])-1 << " " << hfid2eid(DT->sibhfs[eids[0]][1])-1 << " " << hfid2eid(DT->sibhfs[eids[0]][2])-1 << " " << hfid2eid(DT->sibhfs[eids[0]][3])-1 << endl;
    cout << "elem " << eids[0] << " oppfaces: " << hfid2lid(DT->sibhfs[eids[0]][0])-1 << " " << hfid2lid(DT->sibhfs[eids[0]][1])-1 << " " << hfid2lid(DT->sibhfs[eids[0]][2])-1 << " " << hfid2lid(DT->sibhfs[eids[0]][3])-1 << endl;

    // tet 2
    DT->elems[eids[1]] = {vid, tet[Faces[1][0]], tet[Faces[1][1]],tet[Faces[1][2]]};
    hfid = sib[1];
    DT->sibhfs[eids[1]] = {elids2hfid(eids[0]+1, 1), elids2hfid(eids[3]+1, 4), hfid, elids2hfid(eids[2]+1, 2)};
    if (hfid2eid(hfid) > 0){
        DT->sibhfs[hfid2eid(hfid)-1][hfid2lid(hfid)-1] = elids2hfid(eids[1]+1, 3);
        push_stack(&head, elids2hfid(eids[1]+1, 3));
        DT->on_boundary[eids[1]] = false;
    } else {
        DT->on_boundary[eids[1]] = true;
    }
    xs[0] = DT->coords[DT->elems[eids[1]][0]];
    xs[1] = DT->coords[DT->elems[eids[1]][1]];
    xs[2] = DT->coords[DT->elems[eids[1]][2]];
    xs[3] = DT->coords[DT->elems[eids[1]][3]];
    if (!check_jac(xs)){
        cout << "negative jacobian at " << eids[1] << endl;
    }
    cout << "elem " << eids[1] << " oppelems: " << hfid2eid(DT->sibhfs[eids[1]][0])-1 << " " << hfid2eid(DT->sibhfs[eids[1]][1])-1 << " " << hfid2eid(DT->sibhfs[eids[1]][2])-1 << " " << hfid2eid(DT->sibhfs[eids[1]][3])-1 << endl;
    cout << "elem " << eids[1] << " oppfaces: " << hfid2lid(DT->sibhfs[eids[1]][0])-1 << " " << hfid2lid(DT->sibhfs[eids[1]][1])-1 << " " << hfid2lid(DT->sibhfs[eids[1]][2])-1 << " " << hfid2lid(DT->sibhfs[eids[1]][3])-1 << endl;


    // tet 3
    DT->elems[eids[2]] = {vid, tet[Faces[2][0]], tet[Faces[2][1]],tet[Faces[2][2]]};
    hfid = sib[2];
    DT->sibhfs[eids[2]] = {elids2hfid(eids[0]+1, 2), elids2hfid(eids[1]+1, 4), hfid, elids2hfid(eids[3]+1, 2)};
    if (hfid2eid(hfid) > 0){
        DT->sibhfs[hfid2eid(hfid)-1][hfid2lid(hfid)-1] = elids2hfid(eids[2]+1, 3);
        push_stack(&head, elids2hfid(eids[2]+1, 3));
        DT->on_boundary[eids[2]] = false;
    } else {
        DT->on_boundary[eids[2]] = true;
    }
    xs[0] = DT->coords[DT->elems[eids[2]][0]];
    xs[1] = DT->coords[DT->elems[eids[2]][1]];
    xs[2] = DT->coords[DT->elems[eids[2]][2]];
    xs[3] = DT->coords[DT->elems[eids[2]][3]];
    if (!check_jac(xs)){
        cout << "negative jacobian at " << eids[2] << endl;
    }
    cout << "elem " << eids[2] << " oppelems: " << hfid2eid(DT->sibhfs[eids[2]][0])-1 << " " << hfid2eid(DT->sibhfs[eids[2]][1])-1 << " " << hfid2eid(DT->sibhfs[eids[2]][2])-1 << " " << hfid2eid(DT->sibhfs[eids[2]][3])-1 << endl;
    cout << "elem " << eids[2] << " oppfaces: " << hfid2lid(DT->sibhfs[eids[2]][0])-1 << " " << hfid2lid(DT->sibhfs[eids[2]][1])-1 << " " << hfid2lid(DT->sibhfs[eids[2]][2])-1 << " " << hfid2lid(DT->sibhfs[eids[2]][3])-1 << endl;


    // tet 4
    DT->elems[eids[3]] = {vid, tet[Faces[3][0]], tet[Faces[3][1]],tet[Faces[3][2]]};
    hfid = sib[3];
    DT->sibhfs[eids[3]] = {elids2hfid(eids[0]+1, 4), elids2hfid(eids[2]+1, 4), hfid, elids2hfid(eids[1]+1, 2)};
    if (hfid2eid(hfid) > 0){
        DT->sibhfs[hfid2eid(hfid)-1][hfid2lid(hfid)-1] = elids2hfid(eids[3]+1, 3);
        push_stack(&head, elids2hfid(eids[3]+1, 3));
        DT->on_boundary[eids[3]] = false;
    } else {
        DT->on_boundary[eids[3]] = true;
    }
    xs[0] = DT->coords[DT->elems[eids[3]][0]];
    xs[1] = DT->coords[DT->elems[eids[3]][1]];
    xs[2] = DT->coords[DT->elems[eids[3]][2]];
    xs[3] = DT->coords[DT->elems[eids[3]][3]];
    if (!check_jac(xs)){
        cout << "negative jacobian at " << eids[3] << endl;
    }
    cout << "elem " << eids[3] << " oppelems: " << hfid2eid(DT->sibhfs[eids[3]][0])-1 << " " << hfid2eid(DT->sibhfs[eids[3]][1])-1 << " " << hfid2eid(DT->sibhfs[eids[3]][2])-1 << " " << hfid2eid(DT->sibhfs[eids[3]][3])-1 << endl;
    cout << "elem " << eids[3] << " oppfaces: " << hfid2lid(DT->sibhfs[eids[3]][0])-1 << " " << hfid2lid(DT->sibhfs[eids[3]][1])-1 << " " << hfid2lid(DT->sibhfs[eids[3]][2])-1 << " " << hfid2lid(DT->sibhfs[eids[3]][3])-1 << endl;


    DT->nelems+=4;
    int oppeid,opplid;

    // loop through stack and flip edges
    /*
    while (head != NULL){
        hfid = head->hfid;
        pop_stack(&head);
        eid = hfid2eid(hfid)-1;
        lid = hfid2lid(hfid)-1;
        oppeid = hfid2eid(DT->sibhfs[eid][lid])-1;
        opplid = hfid2lid(DT->sibhfs[eid][lid])-1;
        xs[0] = DT->coords[DT->elems[oppeid][0]];
        xs[1] = DT->coords[DT->elems[oppeid][1]];
        xs[2] = DT->coords[DT->elems[oppeid][2]];
        
        if (inside_circumtet(xs, DT->coords[DT->elems[eid][0]])){
            flip_edge(DT,eid,1);

            hfid = DT->sibhfs[oppeid][1];
            if (hfid2eid(hfid) > 0){
                push_stack(&head, elids2hfid(oppeid+1,2));
            }
            hfid = DT->sibhfs[eid][1];
            if (hfid2eid(hfid) > 0){
                push_stack(&head, elids2hfid(eid+1,2));
            }
        }

    }
    */
    return;
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

int inside_tet(const Mat &xs, const vector<double> &ps){
    vector<vector<int>> Faces = {{0,2,1},{0,1,3},{1,2,3},{2,0,3}};

    vector<double> N = normal3d(xs,0);
    vector<double> a = {xs[3][0]-xs[Faces[0][0]][0], xs[3][1]-xs[Faces[0][0]][1], xs[3][2]-xs[Faces[0][0]][2]};
    double D0 = N[0]*a[0]+N[1]*a[1]+N[2]*a[2];
    vector<double> p = {ps[0]-xs[Faces[0][0]][0], ps[1]-xs[Faces[0][0]][1], ps[2]-xs[Faces[0][0]][2]};
    double D1 = N[0]*p[0]+N[1]*p[1]+N[2]*p[2];
    cout << D0 << " " << D1 << endl;
    bool S1 = D1 < 0;

    N = normal3d(xs,1);
    a = {xs[2][0]-xs[Faces[1][0]][0], xs[2][1]-xs[Faces[1][0]][1], xs[2][2]-xs[Faces[1][0]][2]};
    D0 = N[0]*a[0]+N[1]*a[1]+N[2]*a[2];
    p = {ps[0]-xs[Faces[1][0]][0], ps[1]-xs[Faces[1][0]][1], ps[2]-xs[Faces[1][0]][2]};
    double D2 = N[0]*p[0]+N[1]*p[1]+N[2]*p[2];
    cout << D0 << " " << D2 << endl;
    bool S2 = D2 < 0;

    N = normal3d(xs,2);
    a = {xs[0][0]-xs[Faces[2][0]][0], xs[0][1]-xs[Faces[2][0]][1], xs[0][2]-xs[Faces[2][0]][2]};
    D0 = N[0]*a[0]+N[1]*a[1]+N[2]*a[2];
    p = {ps[0]-xs[Faces[2][0]][0], ps[1]-xs[Faces[2][0]][1], ps[2]-xs[Faces[2][0]][2]};
    double D3 = N[0]*p[0]+N[1]*p[1]+N[2]*p[2];
    cout << D0 << " " << D3 << endl;
    bool S3 = D3 < 0;

    N = normal3d(xs,3);
    a = {xs[1][0]-xs[Faces[3][0]][0], xs[1][1]-xs[Faces[3][0]][1], xs[1][2]-xs[Faces[3][0]][2]};
    D0 = N[0]*a[0]+N[1]*a[1]+N[2]*a[2];
    p = {ps[0]-xs[Faces[3][0]][0], ps[1]-xs[Faces[3][0]][1], ps[2]-xs[Faces[3][0]][2]};
    double D4 = N[0]*p[0]+N[1]*p[1]+N[2]*p[2];
    cout << D0 << " " << D4 << endl;
    bool S4 = D4 < 0;

    cout << S1 << " " << S2 << " " << S3 << " " << S4 << endl;
    if (S1 && S2 && S3 && S4){
        return 0;
    } else if(D1>=0 && (D1>=D2) && (D1>=D3) && (D1>=D4)){
        return 1;
    } else if(D2>=0 && (D2>=D1) && (D2>=D3) && (D2>=D4)){
        return 2;
    } else if(D3>=0 && (D3>=D1) && (D3>=D2) && (D3>=D4)){
        return 3;
    } else if(D4>=0 && (D4>=D1) && (D4>=D2) && (D4>=D3)){
        return 4;
    } else {
        cout << "error in inside_tet returning null" << endl;
        assert(false);
        return -1;
    }
}

bool find_enclosing_tet(Mesh* DT, int* tet, int vid){
    int v1,v2,v3,v4,i,hfid;
    bool stop;
    int iters = 0;
    i = 0;
    stop = false;
    vector<vector<double>> xs = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
    int IT;
    while (!stop){
        v1 = DT->elems[*tet][0];
        v2 = DT->elems[*tet][1];
        v3 = DT->elems[*tet][2];
        v4 = DT->elems[*tet][3];
        cout << "elem: " << *tet << " " << v1 << " " << v2 << " " << v3 << " " << v4 << endl;
        if (DT->delete_elem[*tet]){
            cout << "deleted tet passed ran through in find_enclosing_tet, should not be" << endl;

        }
        xs[0] = DT->coords[v1];
        xs[1] = DT->coords[v2];
        xs[2] = DT->coords[v3];
        xs[3] = DT->coords[v4];
        IT = inside_tet(xs,DT->coords[vid]);
        if (IT==0){
            stop = true;
            return true;
        } else if (IT==1){
            hfid = DT->sibhfs[*tet][0];
            if (hfid != 0){
                *tet = hfid2eid(hfid)-1;
            } else {
                stop = true;
                *tet = elids2hfid(*tet+1,1);
                return false;
            }
        } else if (IT==2){
            hfid = DT->sibhfs[*tet][1];
            if (hfid != 0){
                *tet = hfid2eid(hfid)-1;
            } else {
                stop = true;
                *tet = elids2hfid(*tet+1,1);
                return false;
            }
        } else if (IT==3){
            hfid = DT->sibhfs[*tet][2];
            if (hfid != 0){
                *tet = hfid2eid(hfid)-1;
            } else {
                stop = true;
                *tet = elids2hfid(*tet+1,1);
                return false;
            }
        } else if (IT==4){
            hfid = DT->sibhfs[*tet][3];
            if (hfid != 0){
                *tet = hfid2eid(hfid)-1;
            } else {
                stop = true;
                *tet = elids2hfid(*tet+1,1);
                return false;
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

bool check_jac(Mat xs){
    const Mat dphi = {{-1,-1,-1},{1,0,0},{0,1,0},{0,0,1}};
    Mat J = Transpose(xs)*dphi;
    bool check = true;
    double detJ = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0]) + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    if (detJ < 0){;
        check = false;
    }
    return check;
}

bool check_jacobians_tet(Mesh* DT){
    const Mat dphi = {{-1,-1,-1},{1,0,0},{0,1,0},{0,0,1}};
    Mat ps = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
    Mat J;
    double detJ;
    bool check = true;
    for (int i = 0; i<DT->nelems; i++){
        if (!DT->delete_elem[i]){
        ps[0] = DT->coords[DT->elems[i][0]];
        ps[1] = DT->coords[DT->elems[i][1]];
        ps[2] = DT->coords[DT->elems[i][2]];
        ps[3] = DT->coords[DT->elems[i][3]];
        J = Transpose(ps)*dphi;
        detJ = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0]) + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
        if (detJ < 0){
            cout << "negative jacobian at eid: " << i << " Jacobian " << detJ << endl;
            check = false;
        }
        }
    }
    return check;
}

bool Compare_faces(Mesh* DT, int eid, int lid, int oppeid, int opplid){
    vector<vector<int>> Faces = {{0,2,1},{0,1,3},{1,2,3},{2,0,3}};
    vector<int> f1 = {DT->elems[eid][Faces[lid][0]], DT->elems[eid][Faces[lid][1]], DT->elems[eid][Faces[lid][2]]};
    vector<int> f2 = {DT->elems[oppeid][Faces[opplid][0]], DT->elems[oppeid][Faces[opplid][1]], DT->elems[oppeid][Faces[opplid][2]]};

    bool v1 = f1[0] == f2[0] || f1[0]==f2[1] || f1[0]==f2[2];
    bool v2 = f1[1] == f2[0] || f1[1]==f2[1] || f1[1]==f2[2];
    bool v3 = f1[2] == f2[0] || f1[2]==f2[1] || f1[2]==f2[2];
    bool check = v1&&v2&&v3;
    if (!check){
        cout << f1[0] << " " << f1[1] << " " << f1[2] << endl;
        cout << f2[0] << " " << f2[1] << " " << f2[2] << endl;
        cout << DT->elems[eid][0] << " " << DT->elems[eid][1] << " " << DT->elems[eid][2] << " " << DT->elems[eid][3] << endl; 
        cout << DT->elems[oppeid][0] << " " << DT->elems[oppeid][1] << " " << DT->elems[oppeid][2] << " " << DT->elems[oppeid][3] << endl; 
    }
    return v1&&v2&&v3;
}

bool check_sibhfs_tet(Mesh* DT){
    int nelems = DT->nelems;
    int hfid,eid,lid;
    bool check = true;
    for (int i = 0; i<nelems; i++){
        if (!DT->delete_elem[i]){
        for (int j = 0; j<4; j++){
            hfid = DT->sibhfs[i][j];
            if (hfid != 0){
            eid = hfid2eid(hfid)-1;
            lid = hfid2lid(hfid)-1;
            if (eid > nelems || lid > 3 || eid < 0 || DT->delete_elem[eid]){
                cout << "sibhfs is wrong at elem: " << i << " face: " << j << " oppeid: " << eid << " opplid: " << lid << endl;
                check = false;
            }  
            if (!Compare_faces(DT,i,j,eid,lid)){
                cout << "sibhfs is wrong at elem: " << i << " face: " << j << " oppeid: " << eid << " opplid: " << lid << endl;
                check = false;
            }
            }
        }
        }
    }
    return check;
} 