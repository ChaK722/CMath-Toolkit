# CMath-Toolkit

A C++20 matrix library providing basic linear algebra and numerical methods.

## Requirements

- C++20 (concepts, `std::is_arithmetic_v`)
- GCC (MinGW-w64) or any C++20 compiler

## Build

```bash
make
./cmath_test
```

## Supported Types

`float` and `double` are instantiated by default. To add more types, append to `src/matrix.cpp`:

```cpp
template class Matrix<int>;
```

## API Reference

### Construction

```cpp
Matrix<double> A(3, 3);           // 3x3 zero matrix
Matrix<double> B(A);              // copy
Matrix<double> C = Matrix<double>::identity(3);  // 3x3 identity
```

### Element Access & Properties

```cpp
A(0, 1) = 5.0;                   // bounds-checked write
double v = A(0, 1);              // bounds-checked read
int r = A.get_row();             // number of rows
int c = A.get_col();             // number of columns
```

### Arithmetic

```cpp
A + B;   A - B;   A * B;         // matrix-matrix
A + 2;   A * 3;   A / 4.0;       // matrix-scalar
2 * A;                            // scalar-matrix (friend)
A += B;  A -= B;  A *= B;        // in-place
```

### Comparison & Output

```cpp
A == B;   A != B;
std::cout << A;                   // formatted output
```

### Linear Algebra

```cpp
A.transpose();                    // returns transposed matrix

auto [L, U] = A.lu_decomposition();  // Doolittle LU decomposition
double d = A.determinant();           // product of U diagonal
Matrix<double> x = A.solve(b);        // solve Ax = b (forward + back substitution)
```

### Numerical Methods

```cpp
// Numerical Jacobian via central finite difference
Matrix<double> J = Matrix<double>::jacobian(func, x);

// Newton-Raphson root finding
Matrix<double> root = Matrix<double>::newton_raphson(func, guess);
// Optional: newton_raphson(func, guess, tol, max_iter)
```

## Example: Nonlinear System Solver

```cpp
#include <iostream>
#include "matrix.h"

Matrix<double> system(const Matrix<double>& X) {
    double x = X(0, 0), y = X(1, 0);
    Matrix<double> F(2, 1);
    F(0, 0) = x * x + y * y - 5.0;   // x^2 + y^2 = 5
    F(1, 0) = x * x - y - 3.0;       // x^2 - y = 3
    return F;
}

int main() {
    Matrix<double> guess(2, 1);
    guess(0, 0) = 2.5;   // initial guess: x ~ 2.5
    guess(1, 0) = 1.5;   //                 y ~ 1.5

    Matrix<double> root = Matrix<double>::newton_raphson(system, guess);
    std::cout << "Root: " << root << "\n";
    std::cout << "F(root) ~ 0: " << system(root) << "\n";
}
```

## Project Structure

```
CMath-Toolkit/
├── include/
│   └── matrix.h        # Template declarations
├── src/
│   ├── matrix.cpp      # Template implementations + explicit instantiations
│   └── main.cpp         # Example usage
├── Makefile
└── README.md
```
