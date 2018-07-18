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

	//The c_i's
	double newc1 = mat.getElement(1,2) / mat.getElement(1,2);
	mat.setElement(1,2,newc1);
	for (int i = 2; i <= mat.getDimensionM(); ++i)
	{
		double newci = (mat.getElement(i,i+1)/*ci*/) / ( mat.getElement(i,i)/*bi*/ - mat.getElement(i-1,i)/*ci'-1*/ * mat.getElement(i,i-1)  );
		mat.setElement(i,i+1, newci);
	}
	//The d_i's
	double newd1 = vec.getElement(1,1) / mat.getElement(1,2);
	vec.setElement(1,1,newd1);
	for (int i = 2; i <= vec.getDimensionN(); ++i)
	{
		double newdi = ( vec.getElement(1,i)/*di*/ - vec.getElement(1,i-1)/*di-1*/  *  mat.getElement(i,i-1)/*ai*/ ) /
									 ( mat.getElement(i,i)/*bi*/ - mat.getElement(i-1,i)/*c'i*/   *  mat.getElement(i,i-1)/*ai*/ );
		vec.setElement(1,i, newdi);
	}

	Matrix solVec(1, mat.getDimensionN());
	solVec.setElement(1,solVec.getDimensionM(), vec.getElement(1,vec.getDimensionM()));

	for (int i = mat.getDimensionN() - 1; i > 0; --i)
	{
		double newxi = vec.getElement(1,i) - mat.getElement(i-1,i) * solVec.getElement(1,i+1);
		solVec.setElement(1,i,newxi);
	}
	return solVec;

}


int main()
{
	std::cout<<"Wie zur Hoelle ging C++ nochmal? \n";
	Matrix m(5,5);
	m.printMat();
	
	m.setElement(1,1,9);
	m.printMat();

	m.setElement(5,5,8);
	m.printMat();

	std::cout<<std::endl;

	m.makeBand3(1,2,3);
	m.printMat();

	UserInput u;

	//u.userDialog_init();
	//std::cout<<"g: "<<u.calculate_g()<<"f: "<<u.calculate_f()<<std::endl;




	std::cout<<std::endl;
	return 0;
}
