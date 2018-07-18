#include "Matrix.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>


Matrix::Matrix(int n, int m): dimN(n),dimM(m) //Empty
{
	val.resize(n,std::vector<double>(m));
	for (int i = 0; i < dimN; ++i)
	{
		for (int j = 0; j < dimM; ++j)
		{
			if(i == j)
			{
				val[i][j] = 1;	
			}
		}
	}

}

void Matrix::makeBand3(double left, double middle, double right)//zeros = anzahl der 0en in der ersten zeile bevor wert left kommt
{
	val[0][0] = middle;
	for (int i =1; i < dimN; ++i)
	{
		val[i-1][i] = right;
		val[i][i] = middle;
		val[i][i-1] = left;
	}
}

void Matrix::printMat()
{
	std::stringstream str;

	for (int i = 0; i < dimN; ++i)
	{
		for (int j = 0; j < dimM; ++j)
		{
			str << val[i][j] ;	

			str<< "--";
		}
		str<<"\n";
	}

	std::cout<< str.str()<<std::endl;
}

double Matrix::getElement(int i, int j)
{
	return val[i-1][j-1];
}
void   Matrix::setElement(int i, int j, float a)
{
	val[i-1][j-1] = a;
}
int Matrix::getDimensionN()
{
	return dimN;
}
int Matrix::getDimensionM()
{
	return dimM;
}
bool Matrix::checkIfBand3()
{
	for (int i = 0; i < dimN; ++i)
	{
		for (int j = 0; j < dimM; ++j)
		{
			if(i - j > 1 or j - i > 1)
			{
				if(val[i][j] > 0.001)return false;
			}	
		}
	}
	return true;
}

