#include <vector>
#include <cmath>
#include <iostream>

class Matrix {
public:
    Matrix(int rows, int cols) : rows_(rows), cols_(cols), data_(rows* cols) {}

    int rows() const { return rows_; }
    int cols() const { return cols_; }

    double& operator()(int i, int j) { return data_[i * cols_ + j]; }
    const double& operator()(int i, int j) const { return data_[i * cols_ + j]; }

    Matrix operator*(const Matrix& other) const {
        if (cols_ != other.rows()) {
            throw std::runtime_error("Incompatible matrix dimensions for multiplication");
        }

        Matrix result(rows_, other.cols());
        for (int i = 0; i < result.rows(); ++i) {
            for (int j = 0; j < result.cols(); ++j) {
                for (int k = 0; k < cols_; ++k) {
                    result(i, j) += (*this)(i, k) * other(k, j);
                }
            }
        }

        return result;
    }

    Matrix operator+(const Matrix& other) const {
        if (rows_ != other.rows() || cols_ != other.cols()) {
            throw std::runtime_error("Incompatible matrix dimensions for addition");
        }

        Matrix result(rows_, cols_);
        for (int i = 0; i < result.rows(); ++i) {
            for (int j = 0; j < result.cols(); ++j) {
                result(i, j) = (*this)(i, j) + other(i, j);
            }
        }

        return result;
    }

    Matrix operator*(double scalar) const {
        Matrix result(rows_, cols_);
        for (int i = 0; i < result.rows(); ++i) {
            for (int j = 0; j < result.cols(); ++j) {
                result(i, j) = (*this)(i, j) * scalar;
            }
        }

        return result;
    }

    Matrix transpose() const {
        Matrix result(cols_, rows_);
        for (int i = 0; i < result.rows(); ++i) {
            for (int j = 0; j < result.cols(); ++j) {
                result(i, j) = (*this)(j, i);
            }
        }

        return result;
    }

    double dot(const Matrix& other) const {
        if (rows_ != 1 || other.rows() != 1 || cols_ != other.cols()) {
            throw std::runtime_error("Incompatible matrix dimensions for dot product");
        }

        double result = 0;
        for (int i = 0; i < cols_; ++i) {
            result += (*this)(0, i) * other(0, i);
        }
        return result;
    }

    Matrix operator*(double scalar) const {
        Matrix result(rows_, cols_);
        for (int i = 0; i < result.rows(); ++i) {
            for (int j = 0; j < result.cols(); ++j) {
                result(i, j) = (*this)(i, j) * scalar;
            }
        }

        return result;
    }

    double maxOffDiagonal() const {
        if (rows_ != cols_) {
            throw std::runtime_error("Matrix must be square for off-diagonal elements");
        }

        double maxVal = std::numeric_limits<double>::lowest();
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                if (i != j) {
                    maxVal = std::max(maxVal, (*this)(i, j));
                }
            }
        }

        return maxVal;
    }

private:
    int rows_;
    int cols_;
    std::vector<double> data_;
    };