#ifndef MATRIX_H
#define MATRIX_H

#include <cmath>
#include <concepts>
#include <functional>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

template <typename T>
concept Number = std::is_arithmetic_v<T>;

template <Number T>
class Matrix {
  int row, col;
  std::vector<T> data;

 public:
  // Construction and destruction
  Matrix(int row, int col);
  Matrix(const Matrix& other);
  Matrix(Matrix&& other) noexcept;

  // Assignment
  Matrix& operator=(const Matrix& other);
  Matrix& operator=(Matrix&& other) noexcept;

  // In-place arithmetic
  Matrix& operator+=(const Matrix& other);
  Matrix& operator-=(const Matrix& other);
  Matrix& operator*=(const Matrix& other);
  Matrix& operator*=(const T& scalar);

  // Arithmetic (non-modifying)
  Matrix operator+(const Matrix& other) const;
  Matrix operator-(const Matrix& other) const;
  Matrix operator*(const Matrix& other) const;
  Matrix operator*(const T& scalar) const;
  Matrix operator/(const T& scalar) const;

  // Scalar * Matrix
  friend Matrix operator*(const T& scalar, const Matrix& matrix) {
    Matrix result(matrix);
    return result *= scalar;
  }

  // Static factory
  static Matrix identity(int n);

  // Accessors
  int get_row() const;
  int get_col() const;
  T& operator()(int i, int j);
  const T& operator()(int i, int j) const;

  // Comparison
  bool operator==(const Matrix& other) const;
  bool operator!=(const Matrix& other) const;

  // Linear algebra
  Matrix transpose() const;
  std::pair<Matrix, Matrix> lu_decomposition() const;
  T determinant() const;
  Matrix solve(const Matrix& b) const;

  // Numerical methods
  static Matrix jacobian(const std::function<Matrix(const Matrix&)>& func,
                         const Matrix& x);
  static Matrix newton_raphson(
      const std::function<Matrix(const Matrix&)>& func,
      const Matrix& x0, T tol = T(1e-6), int max_iter = 100);

  // Output
  friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
    for (int i = 0; i < matrix.row; i++) {
      for (int j = 0; j < matrix.col; j++) {
        os << matrix.data[i * matrix.col + j] << " ";
      }
      os << std::endl;
    }
    return os;
  }
};

#endif  // MATRIX_H
