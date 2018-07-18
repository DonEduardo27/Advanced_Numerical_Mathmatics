#include <iostream>
#include <string>
#include <sstream>

#include <vector>

#include "Matrix.hpp"
#include "UserInput.hpp"


Matrix thomas(Matrix mat, Matrix vec)
{
	//ERROR Messages
	if(mat.getDimensionM() != mat.getDimensionN())
	{
		std::cout<<"Error in Thomas algorithm: Only possible square matrices!\n";
		Matrix m(0,0);
		return m;
	}
	if (!mat.checkIfBand3())
	{
		std::cout<<"Error in Thomas algorithm: Only possible on 3 Band matrices!\n";
		Matrix m(0,0);
		return m;
	}
	if(vec.getDimensionM() != 1)
	{
		std::cout<<"Error in Thomas algorithm: Left hand side has to be a column vector! (Did you try to transpose?)\n";
		Matrix m(0,0);
		return m;		
	}
	if(vec.getDimensionN() != mat.getDimensionN())
	{
		std::cout<<"Error in Thomas algorithm: vector- and Matrix dimension have to match (Did you try to transpose?)\n";
		Matrix m(0,0);
		return m;		
	}
	//Actual algorithm
	//(1)
	Matrix alpha(vec.getDimensionN()+1,1);
	Matrix beta (vec.getDimensionN()+1,1);
	alpha.setElement(1,1,0);
	beta.setElement(1,1,0);
	//(2)
	for (int i = 1; i <= vec.getDimensionN(); ++i)
	{
		double alphajplus1 = ((-1) * mat.getElement(i,i+1)) /
												(mat.getElement(i,i) + alpha.getElement(i,1) * mat.getElement(i,i-1)); 
		double betajplus1  = (vec.getElement(i,1) - beta.getElement(i,1) * mat.getElement(i,i-1)) /
												(mat.getElement(i,i) + alpha.getElement(i,1) * mat.getElement(i,i-1));
		alpha.setElement(i+1,1,alphajplus1);
		 beta.setElement(i+1,1,betajplus1);
	}
	//(3)
	Matrix solVec (vec.getDimensionN(),1);
	solVec.setElement(vec.getDimensionN(),1, beta.getElement(vec.getDimensionN()+1,1));
	//(4)
	for (int i = vec.getDimensionN() - 1; i > 0; i--)
	{
		solVec.setElement(i,1, alpha.getElement(i+1,1) * solVec.getElement(i+1,1) + beta.getElement(i+1,1));
	}
	return solVec;

}


int main()
{
	std::cout<<"Matrix Solver \n";
	UserInput UInp;

	double f = UInp.calculate_f();
	double g = UInp.calculate_g();

	Matrix m(4,4);
	m.makeBand3(g,f,g);
	m.printMat();

	Matrix v(4,1);
	v.setElement(1,1,8);
	v.setElement(2,1,12);
	v.setElement(3,1,17);
	v.setElement(4,1,20);
	v.printMat();

	thomas(m,v).printMat();		




	std::cout<<std::endl;
	return 0;
}
