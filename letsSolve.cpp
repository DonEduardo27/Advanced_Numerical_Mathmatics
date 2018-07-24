#include <iostream>
#include <string>
#include <sstream>

#include <vector>

#include "Algorithms.hpp"
#include "Matrix.hpp"
#include "UserInput.hpp"






int main()
{
	std::cout<<"Matrix Solver \n";

	//UserInput UInp;
	Algorithms algo;

	Matrix left(3,3);
	Matrix right(3,3);

	left.setElement(1,1, 1);
	left.setElement(2,1, 2);
	left.setElement(3,1, 3);
	
	left.setElement(1,2, 3);
	left.setElement(2,2, 1);
	left.setElement(3,2, 2);
	
	left.setElement(1,3, 1);
	left.setElement(2,3, 3);
	left.setElement(3,3, 2);

	right.setElement(1,1, 3);
	right.setElement(2,1, 2);
	right.setElement(3,1, 3);
	
	right.setElement(1,2, 1);
	right.setElement(2,2, 2);
	right.setElement(3,2, 1);
	
	right.setElement(1,3, 1);
	right.setElement(2,3, 3);
	right.setElement(3,3, 3);
	left.multiply(left,right).printMat();

/*
	Matrix m = UInp.make_matrix();
	m.printMat();

	double f = UInp.calculate_f();
	double g = UInp.calculate_g();

	Matrix solvec(0,0);
	Matrix v = UInp.calculate_rhs(0,solvec);//initial	
	std::cout<<"\n Using Chase-Method, doing "<< UInp.getIterations()<<" Iterations, vector u is:\n";
	for (int i = 1; i <= UInp.getIterations(); ++i)
	{
		solvec = algo.thomas(m,v);
		v = UInp.calculate_rhs(i+1,solvec);
	}
	solvec.printMat();
	std::cout<<"\n Using Gauss elimination:\n";

	for (int i = 1; i <= UInp.getIterations(); ++i)
	{
		solvec = algo.gauss(m,v);
		v = UInp.calculate_rhs(i+1,solvec);
	}
	solvec.printMat();
*/


	std::cout<<std::endl;
	return 0;
}
