#include<linalg.hpp>

#define FAIL "\033[0;31mFAIL\033[0m"
#define PASS "\033[0;32mPASS\033[0m"

void Test_QR_decomp(){
    bool pass = true;
    int rows[4] = {10,20,40,100};
    int cols[4] = {10,20,35,50};
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
    int rows[6] = {10,20, 50,100, 35,50};
    int cols[6] = {10,20, 35,50, 50,100};
    for(int i = 0; i<6; i++){
        int m = rows[i];
        int n = cols[i];
        Mat A = rMat(m,n);
        vec b = rvec(m);
        for (int j = 0; j<m; j++){
            A[j][3] = A[j][6];
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

void Test_QR_inv(){
    bool pass = true;
    int rows[4] = {10,50,100,200};
    int cols[4] = {10,50,100,200};
    for(int i = 0; i<4; i++){
        int m = rows[i];
        int n = cols[i];
        Mat A = rMat(m,n);
        vec b = rvec(m);
        Mat A_inv = QRinv(A);
        vec err = A*A_inv*b - b;
        if (norm(err) < 1e-12){
            cout << "QRCP inverse test " << i+1 <<": " << PASS << ", ";
        } else {
            cout << "QRCP inverse test " << i+1 <<": " << FAIL << ", ";
        }
    }
    cout << endl;
}

void Test_LU_solve(){
    bool pass = true;
    int rows[4] = {10,50,100,200};
    int cols[4] = {10,50,100,200};
    for(int i = 0; i<4; i++){
        int m = rows[i];
        int n = cols[i];
        Mat A = rMat(m,n);
        vec b = rvec(m);
        vec x = LU_solve(A,b,pass);
        vec err = A*x-b;
        if (norm(err) < 1e-10 && pass){
            cout << "LU solve test " << i+1 <<": " << PASS << ", ";
        } else {
            cout << "LU solve test " << i+1 <<": " << FAIL << ", ";
        }
    }
    cout << endl;
}

void Test_LUP_solve(){
    bool pass = true;
    int rows[4] = {10,50,100,200};
    int cols[4] = {10,50,100,200};
    for(int i = 0; i<4; i++){
        int m = rows[i];
        int n = cols[i];
        Mat A = rMat(m,n);
        vec b = rvec(m);
        vec x = LUP_solve(A,b,pass);
        vec err = A*x-b;
        if (norm(err) < 1e-10 && pass){
            cout << "LUP solve test " << i+1 <<": " << PASS << ", ";
        } else {
            cout << "LUP solve test " << i+1 <<": " << FAIL << ", ";
        }
    }
    cout << endl;
}


int main(){
    Test_QR_decomp();
    Test_QRCP_decomp();
    Test_QR_solve();
    Test_QR_inv();
    Test_LU_solve();
    Test_LUP_solve();
    return 0;
}