#include <iostream>
#include <cmath>
#include <functional>
#include "matrix.h"

using namespace std;

// Nonlinear system F(x, y) = 0
//   f1: x^2 + y^2 - 5 = 0
//   f2: x^2 - y - 3   = 0
Matrix<double> nonlinear_system(const Matrix<double>& X) {
    Matrix<double> F(2, 1);
    double x = X(0, 0);
    double y = X(1, 0);

    F(0, 0) = x * x + y * y - 5.0;
    F(1, 0) = x * x - y - 3.0;

    return F;
}

int main() {
    cout << "--- Newton-Raphson Solver Test ---" << endl;

    // Initial guess near (2.5, 1.5)
    Matrix<double> guess(2, 1);
    guess(0, 0) = 2.5;
    guess(1, 0) = 1.5;

    cout << "Initial guess:\n" << guess;

    try {
        Matrix<double> root = Matrix<double>::newton_raphson(nonlinear_system, guess);

        cout << "\nRoot found:\n" << root;

        // Verify that F(root) is close to zero
        cout << "\nF(root) should be near zero:\n" << nonlinear_system(root);

    } catch (const exception& e) {
        cerr << "Solver failed: " << e.what() << endl;
    }

    return 0;
}
