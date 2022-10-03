#include<linalg.hpp>

int main(){
    int n = 20;
    vec u = e_i(n,2);
    Mat A = rMat(n,n);

    printvec(u);
    printMat(A);
    Transpose(A);
    printMat(A);
    return 0;
}