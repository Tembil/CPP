#include "pch.h"
#include "StrassenOMP.h"

void StrassenOMP::Add(const Matrix & a, const Matrix & b, Matrix & result)
{
	const int size = result.getSize();
	for (int x = 0; x < size; x++) {
		for (int y = 0; y < size; y++) {
			result.set(x, y) = a.get(x, y) + b.get(x, y);
		}
	}
}

void StrassenOMP::Sub(const Matrix & a, const Matrix & b, Matrix & result)
{
	const int size = result.getSize();
	for (int x = 0; x < size; x++) {
		for (int y = 0; y < size; y++) {
			result.set(x, y) = a.get(x, y) - b.get(x, y);
		}
	}
}

void StrassenOMP::Strassen(const Matrix & a, const Matrix & b, Matrix & result)
{
	const int size = result.getSize();
	if (size <= 64)
	{
		Matrix::Mul(a, b, result);
		return;
	}

	const int halfSize = size / 2;

	Matrix a11 = Matrix(halfSize);
	Matrix a12 = Matrix(halfSize);
	Matrix a21 = Matrix(halfSize);
	Matrix a22 = Matrix(halfSize);

	Matrix b11 = Matrix(halfSize);
	Matrix b12 = Matrix(halfSize);
	Matrix b21 = Matrix(halfSize);
	Matrix b22 = Matrix(halfSize);

	Matrix r11 = Matrix(halfSize);
	Matrix r12 = Matrix(halfSize);
	Matrix r21 = Matrix(halfSize);
	Matrix r22 = Matrix(halfSize);

	Matrix p1 = Matrix(halfSize);
	Matrix p2 = Matrix(halfSize);
	Matrix p3 = Matrix(halfSize);
	Matrix p4 = Matrix(halfSize);
	Matrix p5 = Matrix(halfSize);
	Matrix p6 = Matrix(halfSize);
	Matrix p7 = Matrix(halfSize);

	// разделяем каждую матрицу на 4 матрицы
#pragma omp parallel for default(shared)
	for (int x = 0; x < halfSize; x++) {
		for (int y = 0; y < halfSize; y++) {
			a11.set(x, y) = a.get(x, y);
			a12.set(x, y) = a.get(x, y + halfSize);
			a21.set(x, y) = a.get(x + halfSize, y);
			a22.set(x, y) = a.get(x + halfSize, y + halfSize);

			b11.set(x, y) = b.get(x, y);
			b12.set(x, y) = b.get(x, y + halfSize);
			b21.set(x, y) = b.get(x + halfSize, y);
			b22.set(x, y) = b.get(x + halfSize, y + halfSize);
		}
	}

	// Calculating p1 to p7:
#pragma omp parallel default(shared)
#pragma omp single
	{
#pragma omp task
		{
			Matrix aResult = Matrix(halfSize);
			Matrix bResult = Matrix(halfSize);
			Add(a11, a22, aResult); // a11 + a22
			Add(b11, b22, bResult); // b11 + b22
			Strassen(aResult, bResult, p1); // p1 = (a11+a22) * (b11+b22)
		}

#pragma omp task
		{
			Matrix aResult = Matrix(halfSize);
			Add(a21, a22, aResult); // a21 + a22
			Strassen(aResult, b11, p2); // p2 = (a21+a22) * (b11)
		}

#pragma omp task
		{
			Matrix bResult = Matrix(halfSize);
			Sub(b12, b22, bResult); // b12 - b22
			Strassen(a11, bResult, p3); // p3 = (a11) * (b12 - b22)
		}

#pragma omp task
		{
			Matrix bResult = Matrix(halfSize);
			Sub(b21, b11, bResult); // b21 - b11
			Strassen(a22, bResult, p4); // p4 = (a22) * (b21 - b11)
		}

#pragma omp task
		{
			Matrix aResult = Matrix(halfSize);
			Add(a11, a12, aResult); // a11 + a12
			Strassen(aResult, b22, p5); // p5 = (a11+a12) * (b22)
		}

#pragma omp task
		{
			Matrix aResult = Matrix(halfSize);
			Matrix bResult = Matrix(halfSize);
			Sub(a21, a11, aResult); // a21 - a11
			Add(b11, b12, bResult); // b11 + b12
			Strassen(aResult, bResult, p6); // p6 = (a21-a11) * (b11+b12)
		}

#pragma omp task
		{
			Matrix aResult = Matrix(halfSize);
			Matrix bResult = Matrix(halfSize);
			Sub(a12, a22, aResult); // a12 - a22
			Add(b21, b22, bResult); // b21 + b22
			Strassen(aResult, bResult, p7); // p7 = (a12-a22) * (b21+b22)
		}

#pragma omp taskwait

#pragma omp task
		{
			Add(p3, p5, r12); // c12 = p3 + p5
			Add(p2, p4, r21); // c21 = p2 + p4
		}

#pragma omp task
		{
			Matrix aResult = Matrix(halfSize);
			Matrix bResult = Matrix(halfSize);
			Add(p1, p4, aResult); // p1 + p4
			Add(aResult, p7, bResult); // p1 + p4 + p7
			Sub(bResult, p5, r11); // c11 = p1 + p4 - p5 + p7
		}

#pragma omp task
		{
			Matrix aResult = Matrix(halfSize);
			Matrix bResult = Matrix(halfSize);
			Add(p1, p3, aResult); // p1 + p3
			Add(aResult, p6, bResult); // p1 + p3 + p6
			Sub(bResult, p2, r22); // c22 = p1 + p3 - p2 + p6
		}

#pragma omp taskwait
	}

	// соединяем матрицы результата
#pragma omp parallel for default(shared)
	for (int x = 0; x < halfSize; x++) {
		for (int y = 0; y < halfSize; y++) {
			result.set(x, y) = r11.get(x, y);
			result.set(x, y + halfSize) = r12.get(x, y);
			result.set(x + halfSize, y) = r21.get(x, y);
			result.set(x + halfSize, y + halfSize) = r22.get(x, y);
		}
	}
}

StrassenOMP::StrassenOMP()
{
}

StrassenOMP::~StrassenOMP()
{
}