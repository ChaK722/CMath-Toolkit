#ifndef MATRIX_H
#define MATRIX_H

#include <concepts>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

template <typename T>
concept Number = std::is_arithmetic_v<T>;

template <Number T>
class Matrix {
  int row, col;
  vector<T> data;

 public:
  Matrix(int row, int col) : row(row), col(col) { data.resize(row * col, T{}); }

  Matrix(const Matrix& other)
      : row(other.row), col(other.col), data(other.data) {};
  Matrix(Matrix&& other) noexcept
      : row(other.row), col(other.col), data(move(other.data)) {};

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
    data = move(other.data);
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
    vector<T> newData(row * other.col, T{});
    for (int i = 0; i < row; i++) {
      for (int j = 0; j < other.col; j++) {
        for (int k = 0; k < col; k++) {
          newData[i * other.col + j] +=
              data[i * col + k] * other.data[k * other.col + j];
        }
      }
    }
    data = move(newData);
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

  int get_row() const { return row; }
  int get_col() const { return col; }

  friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
      for(int i = 0; i < matrix.row; i++) {
        for(int j = 0; j < matrix.col; j++) {
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

  bool operator!=(const Matrix& other) const {
    return !(*this == other);
  }

  T& operator()(int i, int j) {
    if(i >= row || j >= col || i < 0 || j < 0) throw std::out_of_range("Index out of range");
    return data[i * col + j];
  }

  const T& operator()(int i, int j) const {
    if(i >= row || j >= col || i < 0 || j < 0) throw std::out_of_range("Index out of range");
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


  
}

#endif  // MATRIX_H