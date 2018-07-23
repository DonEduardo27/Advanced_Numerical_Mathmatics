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

	std::cout<<"\n---\n";


	double f = UInp.calculate_f();
	double g = UInp.calculate_g();

	Matrix solvec(0,0);


	Matrix v = UInp.calculate_rhs(0,m);//initial	
	for (int i = 1; i <= UInp.getIterations(); ++i)
	{
		
		solvec = algo.thomas(m,v);
		std::cout<<"Iteration "<<i <<":\n";
		solvec.printMat();
		
		v = UInp.calculate_rhs(i,solvec);

	}


	std::cout<<std::endl;
	return 0;
}
