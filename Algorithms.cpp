#include "Algorithms.hpp"
#include <string>
#include <math.h>

bool basicCheck(Matrix mat, Matrix vec, std::string algo)
{
		if(mat.getDimensionM() != mat.getDimensionN())
	{
		std::cout<<"Error in "<<algo<<" algorithm: Only possible square matrices!\n";
		Matrix m(0,0);
		return false;
	}
	if(vec.getDimensionM() != 1)
	{
		std::cout<<"Error in "<<algo<<" algorithm: Left hand side has to be a column vector! (Did you try to transpose?)\n";
		Matrix m(0,0);
		return false;		
	}
	if(vec.getDimensionN() != mat.getDimensionN())
	{
		std::cout<<"Error in "<<algo<<" algorithm: vector- and Matrix dimension have to match (Did you try to transpose?)\n";
		Matrix m(0,0);
		return false;		
	}
	return true;
}

Matrix Algorithms::thomas(Matrix mat, Matrix vec)
{
	//ERROR Messages
	if(!basicCheck(mat,vec,"Thomas"))
	{
		Matrix m(0,0);
		return m;
	}
	if (!mat.checkIfBand3())
	{
		std::cout<<"Error in Thomas algorithm: Only possible on 3 Band matrices!\n";
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

double abs(double a)
{
	if(a<0)return -a;
	return a;
}

Matrix  Algorithms::gauss(Matrix matInp, Matrix vec)
{
	if(!basicCheck(matInp,vec,"Gauss"))
	{
		Matrix m(0,0);
		return m;
	}

	Matrix mat(matInp.getDimensionN(),matInp.getDimensionN()+1);

	for (int i = 1; i <= matInp.getDimensionN(); ++i)
	{
		for (int j = 1; j <= matInp.getDimensionN(); ++j)
		{
			mat.setElement(i,j,matInp.getElement(i,j));
		}
	}

	for (int i = 1; i <= matInp.getDimensionM(); ++i)
	{
		mat.setElement(i,matInp.getDimensionN()+1,vec.getElement(i,1));
	}


//Actual Algorithm

	for (int i = 1; i <= matInp.getDimensionN(); ++i)
	{
		double max  = abs(mat.getElement(i,i));
		int maxRow = i;
		for (int k = i + 1; k <= matInp.getDimensionN(); ++k)//Maximum
		{
			if(abs(mat.getElement(k,i)) > max)
			{
				max = mat.getElement(k,i);
				maxRow = k;
			}
		}
		//Swap function
		for (int k = i; k <= matInp.getDimensionN()+1; ++k)//swap
		{
			double temp = mat.getElement(maxRow,k);
			mat.setElement(maxRow,k, mat.getElement(i,k));
			mat.setElement(i,k, temp);
		}	
		//0er für die dreiecksmatrix
		for (int k = 1+i; k <= matInp.getDimensionN(); ++k)
		{
			double subtractMat = -((mat.getElement(k,i))/(mat.getElement(i,i)));

			for (int j = i; j <= matInp.getDimensionN()+1; ++j)
			{
				if(i==j) mat.setElement(k,j, 0.0);
				else mat.setElement(k,j, mat.getElement(k,j) + mat.getElement(i,j) * subtractMat);
			}
		}
	}

	Matrix vecSol(matInp.getDimensionN(),1);
	for (int i = matInp.getDimensionN(); i > 0; i--)
	{
		double devidedElements = mat.getElement(i,matInp.getDimensionN() + 1) / mat.getElement(i,i);
		vecSol.setElement(i,1, devidedElements);
		for (int k = i-1; k > 0; k--)
		{
			mat.setElement(k,matInp.getDimensionN()+1,  mat.getElement(k,matInp.getDimensionN()+1) - mat.getElement(k,i)* vecSol.getElement(i,1));
		}
	}
	return vecSol;
}

Matrix Algorithms::fixpoint(Matrix mat, Matrix vec, int iterations)
{

	if(!basicCheck(mat,vec, "Fixpoint"))
	{
		Matrix m(0,0);
		return m;
	}

	double f = mat.getElement(1,1);
	double g = mat.getElement(1,2);
	double norm = 2*(g/f);
	if(norm < 0)norm = -norm;
	if(norm > 1)
	{
		std::cout<<"Error in Fixpoint algorithm: Won't converge as norm of C is "<<norm<<" > 1!\n";
		Matrix m(0,0);
		return m;
	}
	//Actual Fixpoint algorithm

	Matrix B(mat.getDimensionN(),mat.getDimensionM());
	B.makeBand3(-(g/f),0,-(g/f));

	int a;

	Matrix C(mat.getDimensionN(),1);
	for (int i = 1; i <= mat.getDimensionN(); ++i)
	{
		C.setElement(i,1,vec.getElement(i,1)/f);
	}
	//u = B u + C
	

	for (int i = 0; i < iterations; ++i)
	{
		vec = (B * vec );
		vec = vec + C;
	}
	return vec;

}
Matrix Algorithms::singleStep(Matrix mat, Matrix vec, int iterations)
{
	if(!basicCheck(mat,vec, "Single Step"))
	{
		Matrix m(0,0);
		return m;
	}

	double f = mat.getElement(1,1);
	double g = mat.getElement(1,2);
	double norm = 2*(g/f);
	if(norm < 0)norm = -norm;
	if(norm > 1)
	{
		std::cout<<"Error in Single Step method: Won't converge as norm of C is "<<norm<<" > 1!\n";
		Matrix m(0,0);
		return m;
	}
	//Actual SS algorithm

	//D
	Matrix D(mat.getDimensionN(),mat.getDimensionM());
	Matrix L(mat.getDimensionN(),mat.getDimensionM());
	Matrix R(mat.getDimensionN(),mat.getDimensionM());

	for (int i = 1; i <= mat.getDimensionN(); ++i)
	{
		D.setElement(i,i, mat.getElement(i,i));
		R.setElement(i,i, 0);
		L.setElement(i,i, 0);
	}
	//L und R


	for (int i = 1; i <= mat.getDimensionN()-1; ++i)
	{
		for (int j = i+1; j <= mat.getDimensionN(); ++j)
		{
			R.setElement(i,j, mat.getElement(i,j));
			L.setElement(mat.getDimensionN()-(i-1),mat.getDimensionN()-(j-1),mat.getElement(i,j));
		}
	}
	D.printMat();
	L.printMat();
	R.printMat();
	//(d+l) invertieren
	Matrix DLInverse(mat.getDimensionN(),mat.getDimensionM());
	double d = D.getElement(1,1);
	double l = L.getElement(2,1); 
	for (int i = 1; i <= mat.getDimensionN(); ++i)
	{
		double diagonalval = (pow(l,i-1) / pow(d,i)) * pow(-1,i-1);
		for (int j = 1; j <= mat.getDimensionN()-(i-1); ++j)
		{
			DLInverse.setElement(j+(i-1),j, diagonalval);
		}
	}
	Matrix test = (D + L)*DLInverse;
	test.printMat();

	/*Matrix B(mat.getDimensionN(),mat.getDimensionM());
	B.makeBand3(-(g/f),0,-(g/f));

	int a;

	Matrix C(mat.getDimensionN(),1);
	for (int i = 1; i <= mat.getDimensionN(); ++i)
	{
		C.setElement(i,1,vec.getElement(i,1)/f);
	}
	//u = B u + C
	

	for (int i = 0; i < iterations; ++i)
	{
		vec = (B * vec );
		vec = vec + C;
	}*/
	return vec;
}
