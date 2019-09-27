// Task_1_OpenMP.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <typeinfo>
#include <omp.h>
#include "StrassenOMP.h"

using namespace std;

void printMatrix(Matrix& matrix) { // печатать матрицу
	int size = matrix.getSize();
	for (int y = 0; y < size; y++) {
		for (int x = 0; x < size; x++) {
			double v = matrix.get(x, y);
			printf("%10.5G  ", v);
		}
		cout << endl;
	}
	cout << endl;
}

Matrix createRandomMatrix(int n) { // создать рандомную матрицу
	Matrix matrix(n);
	int size = matrix.getSize();
	for (int y = 0; y < size; y++) {
		for (int x = 0; x < size; x++) {
			matrix.set(x, y) = rand();
		}
	}
	return matrix;
}

bool equals(Matrix& a, Matrix& b) { // сравнить матрицы
	int N = a.getSize();
	for (int y = 0; y < N; y++) {
		for (int x = 0; x < N; x++) {
			double v1 = a.get(x, y);
			double v2 = b.get(x, y);
			if (abs(v2 - v1) > 0.01) {
				cout << "Error at (" << x << "," << y << ") ";
				cout << "v1 = " << v1 << " v2 = " << v2 << " ";
				cout << "delta = " << (v1 - v2) << endl;
				return false;
			}
		}
	}
	return true;
}

int main(int argc, char** argv) {
	int N = 16;
	srand(time(NULL));

	const Matrix a = createRandomMatrix(N);
	const Matrix b = createRandomMatrix(N);
	Matrix result_strassen_omp(N);

	int start = clock();
	StrassenOMP* strassen = new StrassenOMP();
	strassen->Strassen(a, b, result_strassen_omp);
	int result = clock() - start;
	cout << "it works";

	system("pause");
	return 0;
}