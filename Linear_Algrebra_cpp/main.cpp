#include<linalg.hpp>

#define FAIL "\033[0;31mFAIL\033[0m"
#define PASS "\033[0;32mPASS\033[0m"

void Test_QR_decomp(){
    bool pass = true;
    int rows[4] = {6,20,40,100};
    int cols[4] = {6,20,35,50};
    for(int i = 0; i<4; i++){
        int m = rows[i];
        int n = cols[i];
        Mat Q = Mat(0);
        Mat R = Mat(0);

        Mat A = rMat(m,n);
        
        QR(A,Q,R);
        double nrm = norm_inf(A-Q*R);
        if (nrm < 1e-16*double(m)*double(n)){
            cout << "QR test " << i+1 <<": " << PASS << ", ";
        } else {
            cout << "QR test " << i+1 <<": " << FAIL << ", ";
        }
    }
    cout << endl;
}

void Test_QRCP_decomp(){
    bool pass = true;
    int rows[4] = {12,20,40,100};
    int cols[4] = {10,20,35,50};
    for(int i = 0; i<4; i++){
        int m = rows[i];
        int n = cols[i];
        Mat Q = Mat(0);
        Mat R = Mat(0);
        Mat P = Mat(0);

        Mat A = rMat(m,n);
        for (int j = 0; j<m; j++){
            A[j][2] = A[j][3];
        }
        int rank = QR(A,Q,R,P);
        double nrm = norm_inf(A-Q*R*Transpose(P));
        if (nrm < 1e-16*double(m)*double(n) && rank == n-1){
            cout << "QRCP test " << i+1 <<": " << PASS << ", ";
        } else {
            cout << "QRCP test " << i+1 <<": " << FAIL << ", ";
        }
    }
    cout << endl;
}

void Test_QR_solve(){
    bool pass = true;
    int rows[4] = {10,20,40,100};
    int cols[4] = {10,20,35,50};
    for(int i = 0; i<4; i++){
        int m = rows[i];
        int n = cols[i];
        Mat A = rMat(m,n);
        vec b = rvec(m);
        for (int j = 0; j<m; j++){
            //A[j][3] = A[j][6];
        }
        vec x = QR_solve(A, b, pass);
        if (pass){
            cout << "QRCP solve test " << i+1 <<": " << PASS << ", ";
        } else {
            cout << "QRCP solve test " << i+1 <<": " << FAIL << ", ";
        }
    }
    cout << endl;
}

int main(){
    int m = 6;
    int n = 6;
    Mat A = rMat(m,n);
    for (int j = 0; j<m; j++){
            A[j][2] = A[j][4];
        }
    Test_QR_solve();
    
    return 0;
}