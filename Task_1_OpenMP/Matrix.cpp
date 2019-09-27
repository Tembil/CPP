#include "pch.h"
#include "Matrix.h"

inline double & Matrix::set(int x, int y)
{
	return mat[y * size + x];
}

inline double Matrix::get(int x, int y) const
{
	return mat[y * size + x];
}

int Matrix::getSize() const
{
	return size;
}

Matrix::Matrix()
{
}

Matrix::Matrix(int size)
{
	this->size = size;
	mat = new double[size * size];
}

void Matrix::Mul(const Matrix & a, const Matrix & b, Matrix & result)
{
	const int N = result.getSize();
	for (int x = 0; x < N; x++) {
		for (int y = 0; y < N; y++) {
			double sum = 0;
			for (int k = 0; k < N; k++) {
				sum += a.get(x, k) * b.get(k, y);
			}
			result.set(x, y) = sum;
		}
	}
}

Matrix::~Matrix()
{
	delete[] mat;
}