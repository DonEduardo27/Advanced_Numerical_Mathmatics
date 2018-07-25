#ifndef Matrix_H
#define Matrix_H

#include <vector>

class Matrix
{
public:
	Matrix(int n, int m); //Empty
	void makeBand3(double left, double middle, double right);
	Matrix operator*(Matrix matR);
	Matrix operator+(Matrix matR);


	const void printMat();

	double getElement(int i, int j);
	void   setElement(int i, int j, double a);
	int getDimensionN();
	int getDimensionM();
	bool checkIfBand3();

private:
	std::vector<std::vector<double> > val;
	int dimN; //y 
	int dimM; //x
};

#endif