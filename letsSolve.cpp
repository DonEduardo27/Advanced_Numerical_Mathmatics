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

	Matrix v = UInp.calculate_rhs();
	v.printMat();
/*

	Matrix m(4,4);
	m.makeBand3(1,2,3);
	m.printMat();

	Matrix v(4,1);
	v.setElement(1,1,8);
	v.setElement(2,1,8);
	v.setElement(3,1,10);
	v.setElement(4,1,5);
	v.printMat();	
	algo.thomas(m,v).printMat();
	algo.gauss(m,v).printMat();
	//printSolution(thomas(m,v), UInp)

*/

	double f = UInp.calculate_f();
	double g = UInp.calculate_g();



	algo.thomas(m,v).printMat();		




	std::cout<<std::endl;
	return 0;
}
