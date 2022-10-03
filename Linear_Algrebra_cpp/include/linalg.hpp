#include <vector>
#include <iostream>
#include <random>
using namespace std;

#ifndef LINALG_H
#define LINALG_H
typedef vector<double> vec;
typedef vector<vector<double>> Mat;

// operators
vec operator*(const Mat &A, const vec &b);


// Matrix and vector initilizers
vec rvec(int n, double lower=0.0, double upper=1.0);

Mat rMat(int m, int n, double lower=0.0, double upper=1.0);

Mat Eye(int n);

vec e_i(int n, int i = 1);


// Utilities
void printvec(vec const &A);

void printMat(Mat const &A);

void Transpose(Mat &A);

#endif