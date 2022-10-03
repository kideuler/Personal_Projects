#include "linalg.hpp"
using namespace std;

typedef vector<double> vec;
typedef vector<vector<double>> Mat;


// list of operators

// multiplication operator
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
    default_random_engine re;
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

// identity vector
vec e_i(int n, int i){
    vec e(n);
    if (i > n){
        return e;
    }
    e[i] = 1.0;
    return e;
}

void Transpose(Mat &A){
    int m = A.size();
    int n = A[0].size();
    cout << m << ' '<< n << endl;;
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
}


void printvec(vec const &A){
    for (int i = 0; i<A.size(); i++){
        cout << A[i] << endl;
    }
}

void printMat(Mat const &A){
    for (int i = 0; i<A.size(); i++){
        for (int j = 0; j<A[i].size(); j++){
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}
