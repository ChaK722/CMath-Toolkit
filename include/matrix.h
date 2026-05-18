#ifndef MATRIX_H
#define MATRIX_H

#include <cmath>
#include <concepts>
#include <functional>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

template <typename T>
concept Number = std::is_arithmetic_v<T>;

template <Number T>
class Matrix {
  int row, col;
  std::vector<T> data;

 public:
  Matrix(int row, int col) : row(row), col(col) { data.resize(row * col, T{}); }

  Matrix(const Matrix& other)
      : row(other.row), col(other.col), data(other.data) {};
  Matrix(Matrix&& other) noexcept
      : row(other.row), col(other.col), data(std::move(other.data)) {};

  Matrix& operator=(const Matrix& other) {
    if (this == &other) return *this;
    row = other.row;
    col = other.col;
    data = other.data;
    return *this;
  }

  Matrix& operator=(Matrix&& other) noexcept {
    if (this == &other) return *this;
    row = other.row;
    col = other.col;
    data = std::move(other.data);
    return *this;
  }

  Matrix& operator+=(const Matrix& other) {
    if (row != other.row || col != other.col)
      throw std::invalid_argument("Matrix dimensions do not match");
    for (int i = 0; i < row * col; i++) {
      data[i] += other.data[i];
    }
    return *this;
  }

  Matrix& operator-=(const Matrix& other) { return *this += (other * -1); }

  Matrix operator+(const Matrix& other) const {
    Matrix result(*this);
    return result += other;
  }

  Matrix operator-(const Matrix& other) const {
    Matrix result(*this);
    return result -= other;
  }

  Matrix& operator*=(const Matrix& other) {
    if (col != other.row)
      throw std::invalid_argument("Matrix dimensions do not match");
    std::vector<T> newData(row * other.col, T{});
    for (int i = 0; i < row; i++) {
      for (int j = 0; j < other.col; j++) {
        for (int k = 0; k < col; k++) {
          newData[i * other.col + j] +=
              data[i * col + k] * other.data[k * other.col + j];
        }
      }
    }
    data = std::move(newData);
    col = other.col;
    return *this;
  }

  Matrix& operator*=(const T& scalar) {
    for (int i = 0; i < row * col; i++) {
      data[i] *= scalar;
    }
    return *this;
  }

  friend Matrix operator*(const T& scalar, const Matrix& matrix) {
    Matrix result(matrix);
    return result *= scalar;
  }

  Matrix operator*(const Matrix& other) const {
    Matrix result(*this);
    return result *= other;
  }

  Matrix operator*(const T& scalar) const {
    Matrix result(*this);
    return result *= scalar;
  }

  Matrix operator/(const T& scalar) const {
    if (scalar == 0) {
      throw std::invalid_argument("Division by zero!");
    }
    Matrix result(row, col);
    for (size_t i = 0; i < data.size(); ++i) {
      result.data[i] = data[i] / scalar;
    }
    return result;
  }

  static Matrix identity(int n) {
    Matrix result(n, n);
    for (int i = 0; i < n; i++) {
      result.data[i * n + i] = 1;
    }
    return result;
  }

  int get_row() const { return row; }
  int get_col() const { return col; }

  friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
    for (int i = 0; i < matrix.row; i++) {
      for (int j = 0; j < matrix.col; j++) {
        os << matrix.data[i * matrix.col + j] << " ";
      }
      os << std::endl;
    }
    return os;
  }

  bool operator==(const Matrix& other) const {
    if (row != other.row || col != other.col) {
      return false;
    }
    for (size_t i = 0; i < data.size(); ++i) {
      if (data[i] != other.data[i]) {
        return false;
      }
    }
    return true;
  }

  bool operator!=(const Matrix& other) const { return !(*this == other); }

  T& operator()(int i, int j) {
    if (i >= row || j >= col || i < 0 || j < 0)
      throw std::out_of_range("Index out of range");
    return data[i * col + j];
  }

  const T& operator()(int i, int j) const {
    if (i >= row || j >= col || i < 0 || j < 0)
      throw std::out_of_range("Index out of range");
    return data[i * col + j];
  }

  Matrix transpose() const {
    Matrix result(col, row);
    for (int i = 0; i < row; ++i) {
      for (int j = 0; j < col; ++j) {
        result(j, i) = data[i * col + j];
      }
    }
    return result;
  }

  std::pair<Matrix, Matrix> lu_decomposition() const {
    if (row != col)
      throw std::invalid_argument("Matrix must be square for LU decomposition");
    Matrix l(row, col);
    Matrix u(row, col);

    for (int i = 0; i < row; i++) {
      l(i, i) = 1;
      for (int j = i; j < col; j++) {
        u(i, j) = (*this)(i, j);
        for (int k = 0; k < i; k++) {
          u(i, j) -= l(i, k) * u(k, j);
        }
      }
      for (int j = i + 1; j < col; j++) {
        l(j, i) = (*this)(j, i);
        for (int k = 0; k < i; k++) {
          l(j, i) -= l(j, k) * u(k, i);
        }
        l(j, i) /= u(i, i);
      }
    }
    return std::make_pair(l, u);
  }

  T determinant() const {
    if (row != col)
      throw std::invalid_argument("Matrix must be square to compute determinant");
    T det = 1;
    auto [L, U] = lu_decomposition();
    for (int i = 0; i < row; ++i) {
      det *= U(i, i);
    }
    return det;
  }

  Matrix solve(const Matrix& b) const {
    if (row != col)
      throw std::invalid_argument("Matrix must be square to solve linear system");
    if (b.get_row() != row || b.get_col() != 1)
      throw std::invalid_argument("b must be a column vector matching matrix dimensions");
    std::pair<Matrix, Matrix> lu = lu_decomposition();
    Matrix L = lu.first;
    Matrix U = lu.second;

    Matrix y(col, 1);
    for (int i = 0; i < col; i++) {
      y(i, 0) = b(i, 0);
      for (int j = 0; j < i; j++) {
        y(i, 0) -= L(i, j) * y(j, 0);
      }
    }

    Matrix x(col, 1);
    for (int i = col - 1; i >= 0; i--) {
      x(i, 0) = y(i, 0);
      for (int j = i + 1; j < col; j++) {
        x(i, 0) -= U(i, j) * x(j, 0);
      }
      x(i, 0) /= U(i, i);
    }
    return x;
  }

  static Matrix jacobian(const std::function<Matrix(const Matrix&)>& func,
                         const Matrix& x) {
    int n = x.get_row();
    Matrix fx = func(x);
    int m = fx.get_row();
    Matrix J(m, n);
    T h = T(1e-6);

    for (int j = 0; j < n; j++) {
      Matrix x_plus = x;
      Matrix x_minus = x;
      x_plus(j, 0) += h;
      x_minus(j, 0) -= h;

      Matrix f_plus = func(x_plus);
      Matrix f_minus = func(x_minus);

      for (int i = 0; i < m; i++) {
        J(i, j) = (f_plus(i, 0) - f_minus(i, 0)) / (T(2) * h);
      }
    }
    return J;
  }

  static Matrix newton_raphson(
      const std::function<Matrix(const Matrix&)>& func,
      const Matrix& x0, T tol = T(1e-6), int max_iter = 100) {
    Matrix x = x0;

    for (int iter = 0; iter < max_iter; iter++) {
      Matrix f = func(x);
      T norm = 0;
      for (int i = 0; i < f.get_row(); i++) norm += f(i, 0) * f(i, 0);
      norm = std::sqrt(norm);
      if (norm < tol) break;

      Matrix J = jacobian(func, x);
      Matrix dx = J.solve(f * T(-1));
      x = x + dx;
    }
    return x;
  }
};

#endif  // MATRIX_H
