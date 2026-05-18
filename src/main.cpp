#include <iostream>
#include "matrix.h"

using namespace std;

int main() {
    cout << "--- Linear System Solver (Ax = b) Test ---" << endl;

    // 1. 初始化系数矩阵 A (3x3)
    Matrix<double> A(3, 3);
    A(0,0) = 2; A(0,1) = 1; A(0,2) = 1;
    A(1,0) = 4; A(1,1) = 3; A(1,2) = 3;
    A(2,0) = 8; A(2,1) = 7; A(2,2) = 9;
    
    // 2. 初始化结果列向量 b (3x1)
    Matrix<double> b(3, 1);
    b(0,0) = 4;
    b(1,0) = 10;
    b(2,0) = 24;

    cout << "Matrix A:\n" << A;
    cout << "Vector b:\n" << b;

    // 3. 求解 Ax = b
    try {
        Matrix<double> x = A.solve(b);
        cout << "\nSolution vector x:\n" << x;
        
        // 4. 验算: A * x 是否等于 b
        cout << "\nVerification (A * x should equal b):\n" << A * x;
        
    } catch (const exception& e) {
        cerr << "Error solving system: " << e.what() << endl;
    }

    return 0;
}