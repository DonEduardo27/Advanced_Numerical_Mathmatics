#include "UserInput.hpp"
#include <iostream>
#include <math.h>

UserInput::UserInput()
{
	userDialog_init(0);
}

void UserInput::userDialog_init(bool user)
{
	if(user)
	{
		double tTotal = 0;
		std::cout<<"Initializing Values. \nValue a of heat equation: ";
		std::cin>>a_heat;
		std::cout<<"\n\nSetting parameters.\nLength of rod: ";
		std::cin>>L_length;
		std::cout<<"\nNumber of measuring points n: ";
		
		std::cin>>n_steps;
		while (n_steps < 1)
		{
			std::cout<<"n has to be greater than 0\n";
			std::cin>>n_steps;	
		}

		h_step = L_length / (float)(n_steps-1) ;
		std::cout<<"\nStep size is: "<< h_step <<"\nTotal Time: ";
		std::cin>>tTotal;
		std::cout<<"Number of discrete timesteps:";
		std::cin>>noTimeSteps;
		Tau_step = tTotal / (float)noTimeSteps;
		std::cout<<"Tau = "<<Tau_step;
		std::cout<<"\n\nBoundary condition.\nSet A: ";
		std::cin>>A_bound;
		std::cout<<"\nThanks. Set omega: ";
		std::cin>>omega_bound;
		std::cout<<"\n\n\nComplete.\n\n\n";
	}
	else
	{
		a_heat = 12;
		L_length = 10;
		n_steps = 20;
		h_step = L_length / (float)(n_steps-1) ;
		noTimeSteps = 10;
		Tau_step = 0.0001;
		A_bound = 1;
		omega_bound =1;
	}
}


double UserInput::calculate_f()
{
	return ((1 / Tau_step) -/*war mal +*/ ((2*a_heat*a_heat)/(h_step*h_step))); //uebeltaeter?
}
double UserInput::calculate_g()
{
	return ((a_heat*a_heat) / (h_step*h_step));/*war mal -*/
}
int UserInput::getIterations()
{
	return noTimeSteps;
}
Matrix UserInput::calculate_rhs(int step, Matrix u)
{

	Matrix rhs(n_steps-2,1);//1, 2, 3, ..., n-2, x
	if(step == 0)
	{
		for (int i = 1; i < n_steps-2; ++i)
		{
			rhs.setElement(i,1, (1/Tau_step)  * ( ((float)i * A_bound)/((float)(n_steps-1)) ) );
		}
		double specialvalue = (1/Tau_step)  *        (((n_steps-2)*A_bound)  /  (n_steps-1))  -  calculate_g() * A_bound * cos(omega_bound * Tau_step)   ;
				
	 	rhs.setElement(n_steps-2,1,   specialvalue);
 	}
 	else
 	{
		for (int i = 1; i < n_steps-2; ++i)
		{
			rhs.setElement(i,1,(1/Tau_step)  *u.getElement(i,1) );
		}

		double specialvalue = (1/Tau_step)  * u.getElement(n_steps-2,1) -  calculate_g() * A_bound * cos(omega_bound * Tau_step*(float)(step))   ;
				
	 	rhs.setElement(n_steps-2,1,   specialvalue); 		
 	}
	return rhs;
}
Matrix UserInput::make_matrix()
{   
	if(n_steps > 2)
	{
		Matrix mat(n_steps-2,n_steps-2);
		mat.makeBand3(calculate_g(),calculate_f(),calculate_g());
		return mat;
	}
	else 
	{
		Matrix ERRmat(1,1);
		return ERRmat;
	}
	
}
//~UserInput();
