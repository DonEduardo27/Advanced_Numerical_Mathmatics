#include "Matrix.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>


Matrix::Matrix(int n, int m): dimN(n),dimM(m) //Unity atrix
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
void Matrix::makeBand3(double left, double middle, double right)
{
	val[0][0] = middle;
	for (int i =1; i < dimN; ++i)
	{
		val[i-1][i] = right;
		val[i][i] = middle;
		val[i][i-1] = left;
	}
}
Matrix Matrix::operator*(Matrix matR)
{
	Matrix matL = *this;
	if(matL.getDimensionM() != matR.getDimensionN())
	{
		std::cout<<"To do a Matrix multiplication Dimensions have to match. (Now: left has "<<matL.getDimensionM()<< " colums and right has "<<matR.getDimensionN()<<" rows)\n";
		Matrix m(0,0);
		return m;
	}

	Matrix res(matL.getDimensionN(),matR.getDimensionM());


	for (int rowLeft = 1; rowLeft <= matL.getDimensionN(); ++rowLeft)
	{
		for (int colRight = 1; colRight <= matR.getDimensionM(); ++colRight)
		{
			double sum = 0;
			for (int cumulated = 1; cumulated <= matR.getDimensionN(); ++cumulated)
			{
				sum += matL.getElement(rowLeft, cumulated) * matR.getElement(cumulated, colRight);
			}
			res.setElement(rowLeft,colRight,sum);
		}
	}
	return res;
}
Matrix Matrix::operator+(Matrix matR)
{
	Matrix matL = *this;
	if(matL.getDimensionN() != matR.getDimensionN() or matL.getDimensionM() != matR.getDimensionM())
	{
		std::cout<<"To do a Matrix addition Dimensions have to be equal. (Now: left has "<<matL.getDimensionM()<< " colums and right has "<<matR.getDimensionM()<<" colums,\nleft has "<<matL.getDimensionN()<< " rows and right has "<<matR.getDimensionN()<<" rows)\n";
		Matrix m(0,0);
		return m;
	}	
	Matrix res( matR.getDimensionN(), matR.getDimensionM());
	for (int j = 1; j <= matR.getDimensionN(); ++j)
	{
		for (int i = 1; i <= matR.getDimensionM(); ++i)
		{
			double sum = matR.getElement(j,i) + matL.getElement(j,i);
			res.setElement(j,i,sum);		
		}
	}
	return res;
}

const void Matrix::printMat()
{
	std::stringstream str;

	for (int i = 0; i < dimN; ++i)
	{
		for (int j = 0; j < dimM; ++j)
		{
			str << val[i][j] ;	

			str<< "   ";
		}
		str<<"\n";
	}

	std::cout<< str.str()<<std::endl;
}

double Matrix::getElement(int i, int j)
{
	return val[i-1][j-1];
}
void   Matrix::setElement(int i, int j, double a)
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
				if(val[i][j] > 0.001)return false;//Quick and dirty. Besser wäre den Wert adaptiv zu machen
			}	
		}
	}
	return true;
}

