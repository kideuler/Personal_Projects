#include "GeoComp.hpp"
#include <unistd.h>
using namespace std;

Spline Cubic_spline(const vector<vector<double>> &xs, const vector<bool> &corners);

vector<double> spline_point_segment(Spline* spl, double a, double b, double ratio){
    double arclength=0.0;
    double q[5] = {-(1/3)*sqrt(5+2*sqrt(10/7)), -(1/3)*sqrt(5-2*sqrt(10/7)),0,(1/3)*sqrt(5-2*sqrt(10/7)),(1/3)*sqrt(5+2*sqrt(10/7))};
    double w[5] = {(322-13*sqrt(70))/900.0, (322+13*sqrt(70))/900.0, 125/225, (322+13*sqrt(70))/900.0,(322-13*sqrt(70))/900.0};

    // find total arclength
    for (int k=0; k<5; k++){
        arclength += w[k]*norm(spline_var(spl,q[k]*(b-a)/2 + (a+b)/2, 1));
    }

    // implement full bisection search method later if needed
    if (abs(b-a) > 0.5){
        if (b<a){
            b +=1;
        } else {
            a +=1;
        }
    }
    vector<double> xs = spline_var(spl, (a+b)/2,0);
    return xs;
}

vector<double> spline_var(Spline* spl, double t, int order){
    if (t<0){
        t = 1-t;
    }
    if (t>1){
        t = fmod(t,1);
    }
    int n = (int) floor(double(spl->nv)*t);
    
    if (n==spl->nv){
        n--;
    }
    bool stop = false;
    while (!stop){
        if (n == 0){
            if (spl->params[n+1] >= t){
                stop = true;
            } else {
                n++;
            }
        } else if (n == spl->nv-1){
            if (spl->params[n] <= t){
                stop = true;
            } else {
                n--;
            }
        } else {
            if (spl->params[n] <= t && spl->params[n+1] >= t){
                stop = true;
            } else if (spl->params[n] > t){
                n--;
            } else {
                n++;
            }
        }
    }
    
    // computing function values or derivatives
    double x_j = (t-spl->params[n])/(spl->params[n+1]-spl->params[n]);
    vector<double> xs = {0.0,0.0}; double coeff;
    for (int i = order; i<spl->degree+1; i++){
        coeff = 1.0;
        for (int j = 0; j<order; j++){
            coeff = coeff*double(i-j);
        }
        xs[0] += coeff*pow(x_j,double(i-order))*spl->xweights[n][i];
        xs[1] += coeff*pow(x_j,double(i-order))*spl->yweights[n][i];
    }

    return xs;
}

Spline spline_init(const vector<vector<double>> &xs, const vector<bool> &corners){

    Spline spl;
    
    spl = Cubic_spline(xs, corners);


    double arclength=0.0;
    double q[5] = {-(1/3)*sqrt(5+2*sqrt(10/7)), -(1/3)*sqrt(5-2*sqrt(10/7)),0,(1/3)*sqrt(5-2*sqrt(10/7)),(1/3)*sqrt(5+2*sqrt(10/7))};
    double w[5] = {(322-13*sqrt(70))/900.0, (322+13*sqrt(70))/900.0, 125/225, (322+13*sqrt(70))/900.0,(322-13*sqrt(70))/900.0};

    double a,b,I;
    double* temp = new double[spl.nv+1];
    temp[0] = 0.0;
    for (int i = 0; i<spl.nv; i++){
        a = spl.params[i];
        b = spl.params[i+1];
        I = 0.0;
        for (int j = 0; j<5; j++){
            I += w[j]*norm(spline_var(&spl,q[j]*(b-a)/2 + (a+b)/2, 1));
        }
        I = I*(b-a)/2;
        arclength += I;
        temp[i+1] = arclength;
    }

    for (int i = 1; i<spl.nv+1; i++){
        spl.params[i] = temp[i]/arclength;
    }

    delete temp;
    return spl;
}

Spline Cubic_spline(const vector<vector<double>> &xs, const vector<bool> &corners){

    // Boolean array for which points are corners
    int nv = xs.size();
    vector<bool> flat(nv);
    int ii,jj;
    for (ii = 0; ii<nv; ii++){
        flat[ii] = corners[ii] || flat[ii];
        flat[(ii+nv-1)%nv] = corners[ii] || flat[(ii+nv-1)%nv];
    }

    // Setting up spline data structure
    Spline spl;
    spl.degree = 3;
    spl.nv = nv;
    spl.coords = Zeros(nv,2);
    spl.xweights = Zeros(nv, 4);
    spl.yweights = Zeros(nv, 4);
    spl.params.resize(nv+1);

    vector<vector<double>> A;
    A.resize(nv);
    for (int i = 0; i<nv; i++){
        A[i].resize(nv);
        spl.coords[i] = xs[i];
        spl.params[i] = (double) i / ((double) nv+1);
    }
    spl.params[nv] = 1;
    vector<double> bx;
    bx.resize(nv);
    vector<double> by;
    by.resize(nv);

    ii = 0;
    while (ii < nv){
        if (flat[ii]){
            jj = (ii+1)%nv;
            bx[ii] = xs[jj][0] - xs[ii][0];
            by[ii] = xs[jj][1] - xs[ii][1];
            if (!flat[jj]){
                A[jj][jj] = 2.0;
                A[jj][ii] = 1.0;
                bx[jj] = 2*(xs[jj][0] - xs[ii][0]);
                by[jj] = 2*(xs[jj][1] - xs[ii][1]);
                ii++;
            }

        } else {
            A[ii][ii] = 4.0;
            A[ii][(ii+nv-1)%nv] = 1.0;
            A[ii][(ii+1)%nv] = 1.0;
            bx[ii] = 3*(xs[(ii+1)%nv][0] - xs[(ii+nv-1)%nv][0]);
            by[ii] = 3*(xs[(ii+1)%nv][1] - xs[(ii+nv-1)%nv][1]);
        }
        ii++;
    }

    bool solves;
    vector<double> Dx = LUP_solve(A, bx, solves);
    vector<double> Dy = LUP_solve(A, by, solves);
    if (!solves){
        cout << "spline matrix was not solved properly" << endl;
    }

    for (ii=0; ii<nv; ii++){
        if (flat[ii]){
            spl.xweights[ii][0] = xs[ii][0];
            spl.yweights[ii][0] = xs[ii][1];
            spl.xweights[ii][1] = Dx[ii];
            spl.yweights[ii][1] = Dy[ii];
            spl.xweights[ii][2] = 0.0;
            spl.yweights[ii][2] = 0.0;
            spl.xweights[ii][3] = 0.0;
            spl.yweights[ii][3] = 0.0;
        } else {
            spl.xweights[ii][0] = xs[ii][0];
            spl.yweights[ii][0] = xs[ii][1];
            spl.xweights[ii][1] = Dx[ii];
            spl.yweights[ii][1] = Dy[ii];
            spl.xweights[ii][2] = 3*(xs[(ii+1)%nv][0] - xs[ii][0]) - 2*Dx[ii] - Dx[(ii+1)%nv];
            spl.yweights[ii][2] = 3*(xs[(ii+1)%nv][1] - xs[ii][1]) - 2*Dy[ii] - Dy[(ii+1)%nv];
            spl.xweights[ii][3] = 2*(xs[ii][0] - xs[(ii+1)%nv][0]) + Dx[ii] + Dx[(ii+1)%nv];
            spl.yweights[ii][3] = 2*(xs[ii][1] - xs[(ii+1)%nv][1]) + Dy[ii] + Dy[(ii+1)%nv];
        }
    }

    return spl;
}

