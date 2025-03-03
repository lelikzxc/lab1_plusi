#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <random>

using namespace std;

template <typename T>
class Matrix
{
    vector<T> data;
    unsigned rows, cols;

public:
    Matrix(unsigned rows, unsigned cols, T value = T()) : rows(rows), cols(cols), data(rows * cols, value) {}

    static Matrix Identity(unsigned size)
    {
        Matrix identity(size, size);
        for (unsigned i = 0; i < size; ++i)
        {
            identity(i, i) = 1;
        }
        return identity;
    }

    unsigned numRows() const { return rows; }
    unsigned numCols() const { return cols; }

    Matrix transpose() const
    {
        Matrix transposed(cols, rows);
        for (unsigned i = 0; i < rows; ++i)
        {
            for (unsigned j = 0; j < cols; ++j)
            {
                transposed(j, i) = (*this)(i, j);
            }
        }
        return transposed;
    }

    T &operator()(unsigned row, unsigned col)
    {
        return data[row * cols + col];
    }

    T operator()(unsigned row, unsigned col) const
    {
        return data[row * cols + col];
    }

    pair<Matrix, T> triangularForm(const vector<T> &b) const
    {
        Matrix augmented(*this);
        vector<T> modifiedB(b);
        T determinantFactor = 1;

        for (unsigned k = 0; k < rows; ++k)
        {
            unsigned pivotRow = k;
            for (unsigned i = k + 1; i < rows; ++i)
            {
                if (std::abs(augmented(i, k)) > std::abs(augmented(pivotRow, k)))
                {
                    pivotRow = i;
                }
            }

            if (std::abs(augmented(pivotRow, k)) < 1e-9)
            {
                continue;
            }

            if (pivotRow != k)
            {
                for (unsigned j = 0; j < cols; ++j)
                {
                    swap(augmented(k, j), augmented(pivotRow, j));
                }
                swap(modifiedB[k], modifiedB[pivotRow]);
                determinantFactor = -determinantFactor;
            }

            for (unsigned i = k + 1; i < rows; ++i)
            {
                T factor = augmented(i, k) / augmented(k, k);
                for (unsigned j = k; j < cols; ++j)
                {
                    augmented(i, j) -= factor * augmented(k, j);
                }
                modifiedB[i] -= factor * modifiedB[k];
            }
        }

        return {augmented, determinantFactor};
    }

    T determinant() const
    {
        if (rows != cols)
        {
            throw invalid_argument("Matrix must be square to calculate determinant.");
        }

        auto result = triangularForm(vector<T>(rows, 0));
        Matrix<T> augmented = result.first;
        T determinantFactor = result.second;

        T det = determinantFactor;
        for (unsigned i = 0; i < rows; ++i)
        {
            det *= augmented(i, i);
        }

        return det;
    }

    vector<T> solve(const vector<T> &b) const
    {
        if (b.size() != rows)
        {
            throw invalid_argument("Size of b does not match number of rows.");
        }

        auto result = triangularForm(b);
        Matrix<T> augmented = result.first;
        vector<T> modifiedB = b;

        vector<T> x(rows, 0);
        for (int i = rows - 1; i >= 0; --i)
        {
            if (std::abs(augmented(i, i)) < 1e-9)
            {
                if (std::abs(modifiedB[i]) > 1e-9)
                {
                    throw runtime_error("No solutions.");
                }
                continue;
            }

            x[i] = modifiedB[i];
            for (unsigned j = i + 1; j < cols; ++j)
            {
                x[i] -= augmented(i, j) * x[j];
            }
            x[i] /= augmented(i, i);
        }

        return x;
    }

    void print() const
    {
        for (unsigned i = 0; i < rows; ++i)
        {
            for (unsigned j = 0; j < cols; ++j)
            {
                cout << (*this)(i, j) << " ";
            }
            cout << endl;
        }
    }
};

template <typename T>
Matrix<T> generateLowerTriangularMatrixWithDeterminant(unsigned size, T determinant)
{
    Matrix<T> matrix(size, size);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<T> dis(1.0, 2.0);

    T product = 1.0;
    for (unsigned i = 0; i < size; ++i)
    {
        for (unsigned j = 0; j <= i; ++j)
        {
            if (i == j)
            {
                T value = dis(gen);
                matrix(i, j) = value;
                product *= value;
            }
            else
            {
                matrix(i, j) = dis(gen);
            }
        }
    }

    T scale = pow(determinant / product, 1.0 / size);
    for (unsigned i = 0; i < size; ++i)
    {
        matrix(i, i) *= scale;
    }

    return matrix;
}

#endif // MATRIX_HPP
