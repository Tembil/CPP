#include <omp.h>
#include "iostream"

struct Matrix
{
private:
	int size;
	double* mat;

public:
	inline double& set(int x, int y);
	inline double get(int x, int y) const;
	int getSize() const;
	static void Mul(const Matrix& a, const Matrix& b, Matrix& result);
	Matrix();
	Matrix(int size);
	~Matrix();
};