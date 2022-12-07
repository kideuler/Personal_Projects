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

// copy vector
vec copy(const vec &x){
    int n = x.size();
    vec y(n);
    for (int i = 0; i<n; i++){
        y[i] = x[i];
    }
    return y;
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

Mat operator*(double a, const Mat &A){
    int i,j;
    int m = A.size();
    int n = A[0].size();
    
    Mat prod(m);
    for (i = 0; i<m; i++){
        prod[i].resize(n);
        for (j = 0; j<n; j++){
            prod[i][j] = a*A[i][j];
        }
    }
    return prod;
}

Mat operator/(const Mat &A,double a){
    int i,j;
    int m = A.size();
    int n = A[0].size();
    
    Mat prod(m);
    for (i = 0; i<m; i++){
        prod[i].resize(n);
        for (j = 0; j<n; j++){
            prod[i][j] = A[i][j]/a;
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

vec operator-(vec &u){
    int n = u.size();
    for (int i = 0; i<n; i++){
        u[i] = -u[i];
    }
    return u;
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

vec add_i(vec x, double a, int place){
    int n = x.size();
    vec u(n);
    for (int i = 0; i<n; i++){
        u[i] = x[i];
    }
    u[place] += a;
    return u;
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
    for(int i = 0; i<m; i++){
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

// get row vector from matrix
vec row(const Mat &A, int i){
    int m = A.size();
    assert(i<m);
    int n = A[0].size();
    vec row_(n);
    for (int j = 0; j<n; j++){
        row_[j] = A[i][j];
    }
    return row_;
}

// get column from matrix
vec col(const Mat &A, int j){
    int m = A.size();
    int n = A[0].size();
    assert(j<n);
    vec col_(m);
    for (int i = 0; i<m; i++){
        col_[i] = A[i][j];
    }
    return col_;
}

// utilities

Mat Transpose(Mat &A){
    int m = A.size();
    int n = A[0].size();
    double val;

   
    Mat B = Zeros(n,m);
    for (int i = 0; i<m; i++){
        for (int j = 0; j<n; j++){
            B[j][i] = A[i][j];
        }
    }
    return B;
}


void printvec(vec const &A){
    for (int i = 0; i<A.size(); i++){
        cout << A[i] << endl;
    }
    cout << endl;
}

void printMat(Mat const &A){
    cout << endl;
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
    double tol = 1e-12;
    Q = Eye(m);
    P = Eye(n);
    R = copy(A);

    // square of column norms
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
        if (abs((*c)[itag]) < tol){
            break;
        }
        
        // pivoting 
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
        
        // updating column norms
        for (k = rank; k<n; k++){
            (*c)[k] -= R[rank][k]*R[rank][k];
        }
        rank++;
    }

    delete c;
    for (j = 0; j<n; j++){
        for (i = j+1; i<m; i++){
            R[i][j] = 0;
        }
    }

    return rank;
}

vec upp_tri_inv(const Mat &U, const vec &b, int rank){
    // initilizing values
    int m = U.size();
    int n = U[0].size();
    int i,j;
    double val;
    assert(m >= n);
    vec x = vec(n);
    // checking the matrix is upper triangular
    if (rank < n) {
        for (i = rank-1; i<n; i++){
            assert(abs(U[i][i]) < 1e-12);
        }
    }

    // doing back substitution
    x[rank-1] = b[rank-1]/U[rank-1][rank-1];
    for (i = rank-1; i>-1; i--){
        val = 0.0;
        for (j = rank-1; j>i; j--){
            val += x[j]*U[i][j];
        }
        x[i] = (b[i] - val)/U[i][i];
    }
    return x;
}

vec low_tri_inv(const Mat &L, const vec &b, int rank){
    // initilizing values
    int m = L.size();
    int n = L[0].size();
    int i,j;
    double val;
    assert(m >= n);
    vec x = vec(n);

    // doing back substitution
    x[0] = b[0]/L[0][0];
    for (i = 1; i<rank; i++){
        val = 0.0;
        for (j = 0; j<i; j++){
            val += x[j]*L[i][j];
        }
        x[i] = (b[i] - val)/L[i][i];
    }
    return x;
}


vec QR_solve(const Mat &A, const vec &b, bool &solves) {
    int i,j;
    double val;
    int m = A.size();
    int n = A[0].size();
    assert(m == b.size());
    int rank;

    vec x = vec(0);
    if (m == 2 && n == 2){
        // explicit formula for 2x2 system
        double det = A[0][0]*A[1][1] - A[1][0]*A[0][1];
        Mat A_inv = Zeros(2,2);
        if (abs(det) < 1e-12){
            solves = false;
            return A_inv*b;
        } else {
            A_inv[0][0] = A[1][1]/det;
            A_inv[1][1] = A[0][0]/det;
            A_inv[0][1] = -A[0][1]/det;
            A_inv[1][0] = -A[1][0]/det;
            return A_inv*b;
        }
    }

    // creating buffer matrices
    Mat *Q = new Mat(0); 
    Mat *Q1 = new Mat(0);
    Mat *R = new Mat(0);
    Mat *R1 = new Mat(0);
    Mat *P = new Mat(0);
    Mat *P1 = new Mat(0);
    // overdetermined or square system
    if (m >= n){
        // finding QR factorization
        rank = QR(A,*Q,*R,*P);
        if (m > n || rank < n){
            *Q1 = Zeros(m,rank);
            for (i = 0; i<m; i++){
                for (j = 0; j<rank; j++){
                    (*Q1)[i][j] = (*Q)[i][j];
                }
                for (j = rank; j<n; j++){
                    (*Q)[i][j] = 0.0;
                }
            }

            *P1 = Zeros(n,rank);
            for (i = 0; i<n; i++){
                for (j = 0; j<rank; j++){
                    (*P1)[i][j] = (*P)[i][j];
                }
            }

            *R1 = Zeros(rank,rank);
            for (i = 0; i<rank; i++){
                for (j = 0; j<rank; j++){
                    (*R1)[i][j] = (*R)[i][j];
                }
            }
            for (i = rank; i<n; i++){
                for (j = 0; j<n; j++){
                    (*R)[i][j] = 0.0;
                }
            }
            // back-sub
            vec x = upp_tri_inv(*R1, Transpose(*Q1)*b,rank);
            x = (*P1)*x;
            vec err = (*R1)*Transpose(*P1)*x - Transpose(*Q1)*b;
            if (norm(err) < abs((*R1)[0][0])*1e-12) {
                solves = true;
            } else {
                solves = false;
            }  
        } else {
            vec x = upp_tri_inv(*R, Transpose(*Q)*b,rank);
            x = (*P)*x;
            vec err = (*Q)*(*R)*Transpose(*P)*x - b;
            if (norm(err) < abs((*R)[0][0])*1e-12) {
                solves = true;
            } else {
                solves = false;
            }
            printvec(err);
            return x;
        }
        
        // deleting buffers
        delete Q;
        delete Q1;
        delete R;
        delete R1;
        delete P;
        delete P1;

        return x;
    } else if (m < n) { // underdetermined system
        Mat *At = new Mat(0);
        *At = copy(A);
        int m1 = n;
        int n1 = m;
        *At = Transpose(*At);
        rank = QR(*At,*Q,*R,*P);

        *R = Transpose(*R);
        *Q1 = Zeros(m1,rank);
        for (i = 0; i<m1; i++){
            for (j = 0; j<rank; j++){
                (*Q1)[i][j] = (*Q)[i][j];
            }
            for (j = rank; j<n1; j++){
               (*Q)[i][j] = 0.0;
            }
        }

        *P1 = Zeros(n1,rank);
        for (i = 0; i<n1; i++){
            for (j = 0; j<rank; j++){
                (*P1)[i][j] = (*P)[i][j];
            }
        }

        *R1 = Zeros(rank,rank);
        for (i = 0; i<rank; i++){
            for (j = 0; j<rank; j++){
                (*R1)[i][j] = (*R)[i][j];
            }
        }
        for (i = rank; i<n1; i++){
            for (j = 0; j<n1; j++){
                (*R)[i][j] = 0.0;
            }
        }
        // forw-sub
        vec x = low_tri_inv(*R1, Transpose(*P)*b,rank);
        x = (*Q1)*x;
        vec err = (*P)*(*R1)*Transpose(*Q1)*x - b;
        if (norm(err) < abs((*R1)[0][0])*1e-12) {
            solves = true;
        } else {
            solves = false;
        }  

        delete At;
    }
}


Mat QRinv(const Mat &A){
    int m = A.size();
    int n = A[0].size();
    int i, j;

    // ensure A is a square system
    assert(m == n);
    Mat A_inv = Zeros(m,n);

    // creating buffer matrices
    Mat *Q = new Mat(0); 
    Mat *R = new Mat(0);
    Mat *P = new Mat(0);
    vec *x = new vec(m);
    int rank = QR(A,*Q,*R,*P);
        if (rank == n){
            for (i = 0; i<m; i++){
                *x = upp_tri_inv(*R, row(*Q,i),rank);
                for (j = 0; j<n; j++){
                    A_inv[j][i] = (*x)[j];
                }
            }
            A_inv = (*P)*A_inv;
        } else {
            cout << endl << "rank deficient matrix, returning zeros" << endl;
        }
        
        // deleting buffers
        delete Q;
        delete R;
        delete P;
        delete x;
        return A_inv;
}

vec LU_solve(Mat A, vec b, bool &solves) { // LU solver without pivoting
    solves = true;
    int m = A.size();
    int n = A[0].size();
    assert(m == n);
    double r;
    int i,j,k;
    vec x(n);

    // forward elimination
    for (k = 0; k<n-1; k++){
        if (abs(A[k][k]) < 1e-14){
            solves = false;
        }
        for (i = k+1; i<n; i++){
            r = A[i][k]/A[k][k];
            for (j = 0; j<n; j++){
                A[i][j] = A[i][j] - r*A[k][j];
            }
            b[i] = b[i] - r*b[k];
        }
    }
    // backward solve
    x[n-1] = b[n-1]/A[n-1][n-1];
    for (i = n-1; i>-1; i--){
        x[i] = b[i];
        for (j=i+1; j<n; j++){
            x[i] = x[i] - A[i][j]*x[j];
        }
        x[i] = x[i]/A[i][i];
    }

    return x;
}

vec LUP_solve(Mat A, vec b, bool &solves) { // LU solver with partial pivoting
    solves = true;
    int m = A.size();
    int n = A[0].size();
    assert(m == n);
    double r, r_max;
    int i,j,k,temp;
    vec x(n);
    vector<int> *P = new vector<int>(n);
    vec s(n);

    for (i = 0; i<n; i++){
        s[i] = 0.0;
        for (j = 0; j<n; j++){
            s[i] = max(s[i],abs(A[i][j]));
        }
        (*P)[i] = i; // initialize row pointer
    }

    for (k = 0; k<n-1; k++){
        r_max = 0.0;
        for (i = k; i<n; i++){
            r = abs(A[(*P)[i]][k]/s[(*P)[i]]);
            if (r > r_max) {
                r_max = r;
                j = i;
            }
        }
        temp = (*P)[k];
        (*P)[k] = (*P)[j];
        (*P)[j] = temp;

        for (i = k+1; i<n; i++){
            A[(*P)[i]][k] = A[(*P)[i]][k] / A[(*P)[k]][k];
            for (j = k+1; j<n; j++){
                A[(*P)[i]][j] = A[(*P)[i]][j] - A[(*P)[i]][k]*A[(*P)[k]][j];
            }
        }
    }

    // Forward Elimination
    for (k=0; k<n-1; k++){ 
        for (i=k+1; i<n; i++){
            b[(*P)[i]] = b[(*P)[i]] - A[(*P)[i]][k]*b[(*P)[k]];
        }
    } 

    // backward solve
    x[n-1] = b[n-1]/A[n-1][n-1];
    for (i = n-1; i>-1; i--){
        x[i] = b[(*P)[i]];
        for (j=i+1; j<n; j++){
            x[i] = x[i] - A[(*P)[i]][j]*x[j];
        }
        if (abs(A[(*P)[i]][i]) < 1e-12){
            solves = false;
        }
        x[i] = x[i]/A[(*P)[i]][i];
    }
    delete P;
    return x;
}