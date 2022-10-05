#include<linalg.hpp>


void Test_QR_decomp(){
    bool pass = true;
    int rows[4] = {5,20,40,100};
    int cols[4] = {4,20,35,50};
    for(int i = 0; i<4; i++){
        int m = rows[i];
        int n = cols[i];
        Mat Q = Mat(0);
        Mat R = Mat(0);

        Mat A = rMat(m,n);
        QR(A,Q,R);
        double nrm = norm_inf(A-Q*R);
        if (nrm < 1e-16*double(m)*double(n)){
            cout << "QR test 1: " << "pass" << endl;
        } else {
            cout << "QR test 1: " << "fails" << endl;
        }
    }
}


int main(){
    int m = 10;
    int n = 10;
    Mat A = rMat(m,n);
    for (int i = 0; i<m; i++){
        A[i][2] = 0*A[i][1];
    }
    Mat B = rMat(m,n);
    Mat Q = Mat(0);
    Mat P = Mat(0);
    Mat R = Zeros(m,n);
    int rank = QR(A, Q, R, P);
    Test_QR_decomp();
    return 0;
}