#include "GeoComp.hpp"
using namespace std;

int elids2hfid(int eid, int lid){
    return (eid << 8) + lid;
}

int hfid2eid(int hfid){
    return hfid >> 8;
}

int hfid2lid(int hfid){
    return (hfid & 255);
}

vector<double> min_array(const vector<vector<double>> &xs){
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

vector<double> max_array(const vector<vector<double>> &xs){
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


Triangulation GeoComp_Delaunay_Triangulation(vector<vector<double>> &xs){
    // size checking
    int nv = xs.size();
    int n;
    int ndims = xs[0].size();
    assert(ndims == 2);

    // establishing upper bound
    int upper_bound = 2*nv;

    Triangulation DT;
    DT.coords = Zeros(nv+3,2);
    DT.elems.reserve(upper_bound);
    DT.deletes.reserve(upper_bound);
    DT.sibhfs.reserve(upper_bound);

    vector<double> a = min_array(xs);
    vector<double> b = max_array(xs);

    // reorder points for optimal triangle search O(N*log(N))
    reorder(xs);

    for (n=0; n<nv; n++){
        DT.coords[n] = xs[n];
    }

    DT.coords[nv][0] = a[0];
    DT.coords[nv][1] = a[1];
    DT.coords[nv+1][0] = a[0] + 2*(b[0]-a[0]);
    DT.coords[nv+1][1] = a[1];
    DT.coords[nv+2][0] = a[0];
    DT.coords[nv+2][1] = a[1] + 2*(b[1]-a[1]);

    DT.elems.push_back({nv,nv+1,nv+2});
    DT.deletes.push_back(false);

    // loop inserting each point into triangulation
    int tri = -1;
    DT.nelems = 1;
    for (int n = 0; n<nv; n++){
        for (int ii = 0; ii<DT.nelems; ii++){
            if (!DT.deletes[ii]) {
                if (inside_tri(DT, ii, n)){
                    tri = ii;
                    break;
                }
            }
        }
        // inserting node into the triangulation using Bowyer-Watson algorithm
    }



    return DT;
}

Triangulation Bowyer_watson_insert_point2d(Triangulation DT, int n, int starting_tri){
    vector<vector<int>> edges = {{1,2},{2,3},{3,1}};
    int nbad = 0;
    // find all triangles whos point is inside circumcircle 
}

double detv(const vector<double> &u, const vector<double> &v){
    return (u[0]*v[1] - u[1]*v[0]);
}

bool inside_tri(const Triangulation &DT, int tri, int vertex){
    vector<double> point = DT.coords[vertex];
    vector<double> v0 = DT.coords[DT.elems[tri][0]];
    vector<double> v1 = DT.coords[DT.elems[tri][1]] - v0;
    vector<double> v2 = DT.coords[DT.elems[tri][2]] - v0;
    double D = detv(v1,v2);
    double a = (detv(point,v2) - detv(v0,v2))/D;
    double b = -(detv(point,v1) - detv(v0,v1))/D;
    return (a>0 && b>0 && (a+b)<1);
}

bool inside_circumtri(const Triangulation &DT, int tri, int vid){
    vector<vector<double>> J;
    J.resize(3,vector<double>(3));
    for (int i = 0; i<3; i++){
        for (int j = 0; j<2; j++){
            J[i][j] = DT.coords[DT.elems[tri][i]][j] - DT.coords[vid][j];
        }
        J[i][2] = inner(DT.coords[DT.elems[tri][i]]-DT.coords[vid],DT.coords[DT.elems[tri][i]]-DT.coords[vid]);
    }
    return (det3(J) > 0);
}

double det3(const vector<vector<double>> &J){
    double D = J[0][2]*(J[1][0]*J[2][1] - J[2][0]*J[1][1]) + \
    J[1][2]*(J[2][0]*J[0][1] - J[0][0]*J[2][1]) + \
    J[2][2]*(J[0][0]*J[1][1] - J[1][0]*J[0][1]);
}

vector<double> find_center(const vector<vector<double>> &xs){
    int nv = xs.size();
    int ndims = xs[0].size();
    vector<double> center(ndims);
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
    vector<double> center(ndims);

    center = find_center(xs);

    // measuring some quantity for each point
    vector<double> e = {1,0};
    vector<double> u(ndims);
    for (int n = 0; n<nv; n++){
        u = xs[n] - center;
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





void WriteStl(const Triangulation &msh){
    FILE *fid;
    fid = fopen("test.stl","w");
    fprintf(fid, "solid mesh\n");
    int nelems = msh.elems.size();
    int ndims = msh.coords[0].size();
    for (int ii = 0; ii<nelems; ii++){
        if (ndims == 2){
            fprintf(fid, "\tfacet normal %e %e %e\n",0.0,0.0,1.0);
        }
        fprintf(fid,"\t\touter loop\n");
        for (int jj = 0; jj<msh.elems[0].size(); jj++){
            fprintf(fid,"\t\t\tvertex %e %e %e\n",msh.coords[msh.elems[ii][jj]][0],msh.coords[msh.elems[ii][jj]][1],0.0);
        }
        fprintf(fid, "\t\tendloop\n");
        fprintf(fid, "\tendfacet\n");
    }
    fprintf(fid, "endsolid\n");
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