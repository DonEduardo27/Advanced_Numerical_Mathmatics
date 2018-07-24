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
	UserInput UInp;
	Algorithms algo;

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

	std::cout<<std::endl;
	return 0;
}
