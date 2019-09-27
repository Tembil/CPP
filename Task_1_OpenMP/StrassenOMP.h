#pragma once

#include <omp.h>
#include "iostream"
#include "Matrix.h"

class StrassenOMP
{
private:
	void Add(const Matrix& a, const Matrix& b, Matrix& result);
	void Sub(const Matrix& a, const Matrix& b, Matrix& result);

public:
	void Strassen(const Matrix& a, const Matrix& b, Matrix& result);
	StrassenOMP();
	~StrassenOMP();
};
