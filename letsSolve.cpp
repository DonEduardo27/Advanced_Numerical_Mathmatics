#include <iostream>
#include <string>
#include <sstream>

#include <vector>
#include <iomanip>

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
	Matrix solfix(0,0);
	Matrix solgau(0,0);
	Matrix solsig(0,0);
	Matrix ChaseVec = UInp.calculate_rhs(0,solvec);//initial - Thomas
	Matrix FixptVec= UInp.calculate_rhs(0,solfix);//initial - Fixpoint
	Matrix GaussVev= UInp.calculate_rhs(0,solgau);//initial - Gauss
	Matrix SStepVec= UInp.calculate_rhs(0,solsig);//initial - S.Step

	std::cout<<"\n Using Chase Method:\n";
	for (int i = 1; i <= UInp.getIterations(); ++i)
	{
		solvec = algo.thomas(m,ChaseVec);
		ChaseVec = UInp.calculate_rhs(i+1,solvec);
	}
	solvec.printMat();
	
	std::cout<<"\n Using Fixpoint Iteration:\n";
	for (int i = 0; i < UInp.getIterations(); ++i)
	{
		solfix = algo.fixpoint(m,FixptVec,10);
		if(solfix.getDimensionM() == 0)break;//fixpoint scheitert
		FixptVec= UInp.calculate_rhs(i+1,solfix);
	}
	solfix.printMat();
	std::cout<<"\n Using Gauss elimination:\n";
	for (int i = 0; i < UInp.getIterations(); ++i)
	{
		solgau = algo.gauss(m,GaussVev);
		GaussVev= UInp.calculate_rhs(i+1,solgau);
	}
	solgau.printMat();

	std::cout<<"\n Using Single step Iteration:\n";
	for (int i = 0; i < UInp.getIterations(); ++i)
	{
		solsig = algo.singleStep(m,SStepVec,10);
		if(solsig.getDimensionM() == 0)break;//S.Step scheitert
		SStepVec= UInp.calculate_rhs(i+1,solsig);
	}
	solsig.printMat();


		std::cout <<"     Chase / Thomas |    Fixpoint   |      Gauss     |     Single step \n";
		for (int i = 1; i < m.getDimensionM(); ++i)
		{
			/* code */
        	std::cout <<std::setw(4)<<"u"<<i<<std::setw(13)<<solvec.getElement(i,1)<<"  |  "<< std::setw(11)<< solfix.getElement(i,1) <<"  |  " <<  std::setw(11)<<solgau.getElement(i,1)<<"   |      "<<solsig.getElement(i,1)<<std::endl;
		}
    
	

	std::cout<<std::endl;
	return 0;
}
