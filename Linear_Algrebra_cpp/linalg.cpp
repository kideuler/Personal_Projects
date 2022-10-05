#include "linalg.hpp"
using namespace std;

typedef vector<double> vec;
typedef vector<vector<double>> Mat;


// list of operators

// copy matrix
Mat copy(const Mat &A){
    int m = A.size();
    int n = A[0].size();
    Mat B(m);
    for (int i = 0; i<m; i++){
        B[i].resize(n);
        for (int j = 0; j<n; j++){
            B[i][j] = A[i][j];
        }
    }
    return B;
}

// mat mat multiplication
Mat operator*(const Mat &A, const Mat &B){
    int i,j,k;
    int m1 = A.size();
    int m2 = B.size();
    int n1 = A[0].size();
    int n2 = B[0].size();

    Mat prod(m1);
    if (n1!=m2){
        cout << "matrices incompatible for multiplication" << endl;
        return prod;
    }

    for (i = 0; i<m1;i++){
        prod[i].resize(n2);
        for (j = 0; j<n2; j++){
            prod[i][j] = 0.0;
            for (k = 0; k<n1; k++){
                prod[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return prod;
}

// mat vec multiplication operator
vec operator*(const Mat &A, const vec &b){
    int i,j;
    int m = A.size();
    int n = A[0].size();
    int x = b.size();

    vec prod(m);
    if(x!=n){
        cout << "matrix and vector have incompatible sizes" << endl;
        return prod;
    }

    for(i=0; i<m; i++){
        prod[i] = 0;
        for(j=0; j<n;j++){
            prod[i] += A[i][j]*b[j];
        }
    }
    return prod;
}

// multiply vector on rhs
vec operator*(const vec &b, const Mat &A){
    int i,j;
    int m = A.size();
    int n = A[0].size();
    int x = b.size();

    vec prod(n);
    if(x!=m){
        cout << "matrix and vector have incompatible sizes" << endl;
        return prod;
    }

    for(i=0; i<n; i++){
        prod[i] = 0;
        for(j=0; j<m;j++){
            prod[i] += A[j][i]*b[j];
        }
    }
    return prod;
}

// add scalar to vector
vec operator+(const vec &u, double a){
    int n = u.size();
    vec w(n);
    for (int i = 0; i<n; i++){
        w[i] = u[i]+a;
    }
    return w;
}

// multiply vector and scalar
vec operator*(double a, const vec &u){
    int n = u.size();
    vec w(n);
    for (int i = 0; i<n; i++){
        w[i] = u[i]*a;
    }
    return w;
}

// add vectors
vec operator+(const vec &u, const vec &v){
    int n1 = u.size();
    int n2 = v.size();
    int n = min(n1,n2);
    if (n1 != n2){
        cout << 'vector addition not possible, different sizes';
    }
    vec w(n);
    for (int i = 0; i<n; i++){
        w[i] = u[i]+v[i];
    }
    return w;
}

// subract vectors
vec operator-(const vec &u, const vec &v){
    int n1 = u.size();
    int n2 = v.size();
    int n = min(n1,n2);
    if (n1 != n2){
        cout << 'vector addition not possible, different sizes';
    }
    vec w(n);
    for (int i = 0; i<n; i++){
        w[i] = u[i]-v[i];
    }
    return w;
}

// matrix addition
Mat operator+(const Mat &A, const Mat &B){
    int m1 = A.size();
    int m2 = B.size();
    int n1 = A[0].size();
    int n2 = B[0].size();
    int m = min(m1,m2);
    int n = min(n1,n2);

    if (n1 != n2 || m1!=m2){
        cout << 'vector addition not possible, different sizes';
    }
    Mat C(m);
    for (int i = 0; i<m; i++){
        C[i].resize(n);
        for (int j = 0; j<n; j++){
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    return C;
}

// matrix subtraction
Mat operator-(const Mat &A, const Mat &B){
    int m1 = A.size();
    int m2 = B.size();
    int n1 = A[0].size();
    int n2 = B[0].size();
    int m = min(m1,m2);
    int n = min(n1,n2);

    if (n1 != n2 || m1!=m2){
        cout << 'vector addition not possible, different sizes';
    }
    Mat C(m);
    for (int i = 0; i<m; i++){
        C[i].resize(n);
        for (int j = 0; j<n; j++){
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

double inner(const vec &u, const vec &v){
    int n1 = u.size();
    int n2 = v.size();
    int n = min(n1,n2);
    if (n1 != n2){
        cout << 'vector dot product not possible, different sizes';
    }
    double val = 0.0;
    for(int i = 0; i<n; i++){
        val = val + u[i]*v[i];
    }
    return val;
}

Mat outer(const vec &u, const vec &v){
    int m = u.size();
    int n = v.size();
    Mat prod(m);
    for (int i = 0; i<m; i++){
        prod[i].resize(n);
        for (int j = 0; j<n; j++){
            prod[i][j] = u[i]*v[j];
        }
    }

    return prod;
}

double norm(const vec &u){
    double val = 0;
    for (int i = 0; i<u.size(); i++){
        val  = val + u[i]*u[i];
    }
    val = sqrt(val);
    return val;
}

double norm_inf(const Mat &A){
    int m = A.size();
    int n = A[0].size();
    double nrm = 0.0;
    double mnrm = 0.0;
    for (int j = 0; j<n; j++){
        nrm = 0.0;
        for (int i = 0; i<m; i++){
            nrm += abs(A[i][j]);
        }
        mnrm = max(mnrm, nrm);
    }
    return mnrm;
}

vec operator/(const vec &u, double a) {
    int n = u.size();
    vec w(n);
    for (int i = 0; i<n; i++){
        w[i] = u[i]/a;
    }
    return w;
}

// Create specific matrices

// random vector
vec rvec(int n, double lower, double upper){
    default_random_engine re;
    vec rand(n);
    uniform_real_distribution<double> unif(lower, upper);
    for(int i = 0; i < n; i++){
        rand[i] = unif(re);
    }
    return rand;
}

// random matrix
Mat rMat(int m, int n, double lower, double upper){
    default_random_engine re(time(0));
    Mat rand(m);
    uniform_real_distribution<double> unif(lower, upper);
    for(int i = 0; i < m; i++){
        rand[i].resize(n);
        for(int j = 0; j < n; j++){
            rand[i][j] = unif(re);
        }
        
    }

    return rand;
}

// identity
Mat Eye(int n){
    Mat I(n);
    for(int i = 0; i<n; i++){
        I[i].resize(n);
        I[i][i] = 1.0;
    }
    return I;
}

Mat Zeros(int m, int n){
    Mat A(m);
    for(int i = 0; i<n; i++){
        A[i].resize(n);
    }
    return A;
}

// identity vector
vec e_i(int n, int i){
    vec e(n);
    if (i > n){
        return e;
    }
    e[i] = 1.0;
    return e;
}

// utilities

Mat Transpose(Mat &A){
    int m = A.size();
    int n = A[0].size();
    double val;

    if (n == m){
        for (int i = 0; i<m; i++){
            for (int j = i; j<n; j++){
                val = A[i][j];
                A[i][j] = A[j][i];
                A[j][i] = val;
            }
        }
    } else {
        cout << "must be a squre matrix for efficient in place transpose";
    }
    return A;
}


void printvec(vec const &A){
    for (int i = 0; i<A.size(); i++){
        cout << A[i] << endl;
    }
    cout << endl;
}

void printMat(Mat const &A){
    for (int i = 0; i<A.size(); i++){
        for (int j = 0; j<A[i].size(); j++){
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

double sign(double x){
    double y = 0.0;
    if (x > 0) {
        y = 1.0;
    } else if (x < 0){
        y = -1.0;
    }
    return y;
}





// Linear algebra

// QR decomposition
void QR(const Mat &A, Mat &Q, Mat &R){
    int m = A.size();
    int n = A[0].size();
    int j, k, i, sz;
    double normx,s,u1,tau, dd;
    assert(m >= n);
    
    Q = Eye(m);
    R = Zeros(m,n);
    R = copy(A);

    for (j = 0; j<n; j++){
        sz = m-j;
        vec *r = new vec(sz);
        vec *w = new vec(sz);
        Mat *Rbuf = new Mat(sz);
        Mat *Qbuf = new Mat(m);
        
        // R matrix buffer
        for (k = 0; k<sz; k++){
            (*r)[k] = R[k+j][j];
            (*Rbuf)[k].resize(n);
            for (i = 0; i<n; i++){
                (*Rbuf)[k][i] = R[k+j][i];
            }
        }
        // Q matrix bufffer
        for (k = 0; k<m; k++){
            (*Qbuf)[k].resize(sz);
            for (i = 0; i<sz; i++){
                (*Qbuf)[k][i] = Q[k][i+j];
            }
        }

        // calculating householder vector
        normx = norm(*r);
        s = -sign(R[j][j]);
        u1 = R[j][j] - s*normx;
        *w = (*r)/u1;
        (*w)[0] = 1.0;
        tau = -s*u1/normx;
        dd = inner(*w,*r);
        *r = *r - tau*dd*(*w);
        
        // updating R
        *Rbuf = outer(tau * (*w), *w * *Rbuf);
        for (k = 0; k<sz; k++){
            for (i = 0; i<n; i++){
                R[k+j][i] = R[k+j][i] - (*Rbuf)[k][i];
            }
        }

        // updating Q
        *Qbuf = outer(*Qbuf * *w, tau * *w);
        for (k = 0; k<m; k++){
            for (i = 0; i<sz; i++){
                Q[k][i+j] = Q[k][i+j] - (*Qbuf)[k][i];
            }
        }
        // deleting buffers
        delete r;
        delete w;
        delete Rbuf;
        delete Qbuf;
    }
}

// function to swap rows in a matrix
void swap_cols(Mat *A, int i, int j){
    int m = (*A).size();
    int n = (*A)[0].size();
    double val;
    assert(i <= n && j <= n);

    for (int k = 0; k<m; k++){
        val = (*A)[k][i];
        (*A)[k][i] = (*A)[k][j];
        (*A)[k][j] = val;
    }
}
// QR_CP decomposition
int QR(const Mat &A, Mat &Q, Mat &R, Mat &P){
    int m = A.size();
    int n = A[0].size();
    int j, k, i, sz, rank, itag;
    double normx,s,u1,tau, dd,val;
    assert(m >= n);
    double tol = 1e-15;
    Q = Eye(m);
    P = Eye(n);
    R = copy(A);

    vec *c = new vec(n);
    for (j = 0; j<n; j++){
        (*c)[j] = 0.0;
        for (k = 0; k<m; k++){
            (*c)[j] += R[k][j]*R[k][j];
        }
    }

    rank = 0;
    while (rank < n){
        itag = rank;
        for (k = rank; k < n; k++){
            if ((*c)[k] > (*c)[itag]){
                itag = k;
            }
        }
        
        if (rank != itag){
            swap_cols(&P,rank,itag);
            swap_cols(&R,rank,itag);
        }
        val = (*c)[rank];
        (*c)[rank] = (*c)[itag];
        (*c)[itag] = val;
        
        sz = m-rank;
        vec *r = new vec(sz);
        vec *w = new vec(sz);
        Mat *Rbuf = new Mat(sz);
        Mat *Qbuf = new Mat(m);


        // R matrix buffer
        for (k = 0; k<sz; k++){
            (*r)[k] = R[k+rank][rank];
            (*Rbuf)[k].resize(n);
            for (i = 0; i<n; i++){
                (*Rbuf)[k][i] = R[k+rank][i];
            }
        }
        // Q matrix bufffer
        for (k = 0; k<m; k++){
            (*Qbuf)[k].resize(sz);
            for (i = 0; i<sz; i++){
                (*Qbuf)[k][i] = Q[k][i+rank];
            }
        }

        // calculating householder vector
        normx = norm(*r);
        if (normx > tol){
            s = -sign(R[rank][rank]);
            u1 = R[rank][rank] - s*normx;
            *w = (*r)/u1;
            (*w)[0] = 1.0;
            tau = -s*u1/normx;
            dd = inner(*w,*r);
        } else {
            tau = 0.0;
            dd = 0.0;
        }
        *r = *r - tau*dd*(*w);
        
        // updating R
        *Rbuf = outer(tau * (*w), *w * *Rbuf);
        for (k = 0; k<sz; k++){
            for (i = 0; i<n; i++){
                R[k+rank][i] = R[k+rank][i] - (*Rbuf)[k][i];
            }
        }

        // updating Q
        *Qbuf = outer(*Qbuf * *w, tau * *w);
        for (k = 0; k<m; k++){
            for (i = 0; i<sz; i++){
                Q[k][i+rank] = Q[k][i+rank] - (*Qbuf)[k][i];
            }
        }
        // deleting buffers
        delete r;
        delete w;
        delete Rbuf;
        delete Qbuf;
        rank++;
        for (k = rank; k<n; k++){
            (*c)[k] -= A[rank][k]*A[rank][k];
        }
        
    }

    delete c;
    for (j = 0; j<n; j++){
        for (i = j+1; i<m; i++){
            R[i][j] = 0;
        }
    }

    for (j = 0; j<n; j++){
        if (abs(R[j][j]) < abs(R[1][1])*tol){
            R[j][j] = 0;
            rank--;
        }
    }
    return rank;
}