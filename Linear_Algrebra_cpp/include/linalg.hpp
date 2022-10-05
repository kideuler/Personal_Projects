#include <vector>
#include <iostream>
#include <random>
#include <cassert>
#include <cmath>
using namespace std;

#ifndef LINALG_H
#define LINALG_H
typedef vector<double> vec;
typedef vector<vector<double>> Mat;

// operators
Mat copy(const Mat &A);
Mat operator+(const Mat &A, const Mat &B);
Mat operator-(const Mat &A, const Mat &B);
Mat operator*(const Mat &A, const Mat &B);
vec operator*(const Mat &A, const vec &b);
vec operator*(const vec &b, const Mat &A);
vec operator/(const vec &u, double a);
vec operator*(double a, const vec &u);

vec operator+(const vec &u, const vec &v);
vec operator+(const vec &u, double a);
vec operator-(const vec &u, const vec &v);

double inner(const vec &u, const vec &v);
Mat outer(const vec &u, const vec &v);
double norm(const vec &u);
double norm_inf(const Mat &A);

// Matrix and vector initilizers
vec rvec(int n, double lower=0.0, double upper=1.0);

Mat rMat(int m, int n, double lower=0.0, double upper=1.0);

Mat Eye(int n);

Mat Zeros(int m, int n);

vec e_i(int n, int i = 1);


// Utilities
void printvec(vec const &A);

void printMat(Mat const &A);

Mat Transpose(Mat &A);


// Linear algebra algorithms
void QR(const Mat &A, Mat &Q, Mat &R);
int QR(const Mat &A, Mat &Q, Mat &R, Mat &P);

#endif