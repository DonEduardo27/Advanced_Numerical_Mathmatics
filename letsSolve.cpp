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
/*

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
	left = left + right;
	left.printMat();
*/

	Matrix m = UInp.make_matrix();
	m.printMat();

	double f = UInp.calculate_f();
	double g = UInp.calculate_g();

	Matrix solvec(0,0);
	Matrix solfix(0,0);
	Matrix solgau(0,0);
	Matrix v = UInp.calculate_rhs(0,solvec);//initial - Thomas
	Matrix v2= UInp.calculate_rhs(0,solfix);//initial - Fixpoint
	Matrix v3= UInp.calculate_rhs(0,solgau);//initial - Gauss
	std::cout<<"V DEBUG init\n";
	v.printMat();

	std::cout<<"\n Using Chase Method:\n";
	for (int i = 1; i <= UInp.getIterations(); ++i)
	{
		solvec = algo.thomas(m,v);
		std::cout<<"solvec DEBUG\n"<<i<<"\n";
		solvec.printMat();
		
		v = UInp.calculate_rhs(i+1,solvec);
		std::cout<<"V DEBUG\n"<<i<<"\n";
		v.printMat();
	}
	solvec.printMat();
	/*
	std::cout<<"\n Using Fixpoint Iteration:\n";
	for (int i = 0; i < UInp.getIterations(); ++i)
	{
		solfix = algo.fixpoint(m,v2,10);
		if(solfix.getDimensionM() == 0)break;//fixpoint scheitert
		v2= UInp.calculate_rhs(i+1,solfix);
	}
	solfix.printMat();
	std::cout<<"\n Using Gauss elimination:\n";
	for (int i = 0; i < UInp.getIterations(); ++i)
	{
		solgau = algo.gauss(m,v3);
		v3= UInp.calculate_rhs(i+1,solgau);
	}
	solgau.printMat();
*/



	std::cout<<std::endl;
	return 0;
}
