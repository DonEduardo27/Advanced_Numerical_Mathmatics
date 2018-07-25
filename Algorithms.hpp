#include <iostream>
#include "Matrix.hpp"

class Algorithms
{
public:
	Matrix thomas(Matrix mat, Matrix vec);
	Matrix gauss(Matrix mat, Matrix vec);
	Matrix fixpoint(Matrix mat, Matrix vec, int iterations);
	Matrix singleStep(Matrix mat, Matrix vec, int iterations);
};