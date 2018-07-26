#include "UserInput.hpp"
#include <iostream>
#include <math.h>
#include <limits>

UserInput::UserInput()
{	
	char inp;
	do
	{
		std::cout<<"Use default values which guarantee functionality of iterative methods? (j/n)\n";
		std::cin>>inp;
		if(inp == 'j')userDialog_init(0);
		if(inp == 'n')userDialog_init(1);
	}while(inp != 'j' and inp != 'n');
	
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
			std::cout<<"n has to be greater than 0\nNumber of measuring points n: ";
			std::cin>>n_steps;	
			std::cin.clear();
        	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}
		h_step = L_length / (double)(n_steps-1) ;
		std::cout<<"\nStep size is: "<< h_step <<"\nTotal Time: ";
		std::cin>>tTotal;
		
		char inp;
		do 
		{
			std::cout<<"Number of discrete timesteps:";
			std::cin>>noTimeSteps;
			while (noTimeSteps < 1)
			{
				std::cout<<"t has to be greater than 0\nNumber of discrete timesteps:";
				std::cin>>n_steps;	
				std::cin.clear();
        		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			}
			Tau_step = tTotal / (double)noTimeSteps;
			std::cout<<"Tau = "<<Tau_step;

			if((double)Tau_step >= (h_step*h_step)/(4*(a_heat*a_heat)) )
			{
				std::cout<<"\n\nNOTE: Resulting Matrix will be badly conditioned!\nIts recommended to choose more timesteps (t > "<<((4*(a_heat*a_heat))* tTotal)/(h_step*h_step)<<") !\n\nContinue with t = "<<noTimeSteps<<"? (j/n)\n";
				std::cin>>inp;
			}
			else 
			{
				inp = 'j';
			}
		}while (inp != 'j');

		std::cout<<"\n\nBoundary condition.\nSet A: ";
		std::cin>>A_bound;
		std::cout<<"\nThanks. Set omega: ";
		std::cin>>omega_bound;
		std::cout<<"\n\n\nComplete.\n\n\n";

	}
	else
	{
		//Default values
		a_heat = 12;
		L_length = 10;
		n_steps = 20;
		h_step = L_length / (double)(n_steps-1) ;
		noTimeSteps = 10;
		Tau_step = 0.0001;
		A_bound = 1;
		omega_bound =1;
	}
}


double UserInput::calculate_f()
{
	return ((1 / Tau_step) - ((2*a_heat*a_heat)/(h_step*h_step)));
}
double UserInput::calculate_g()
{
	return ((a_heat*a_heat) / (h_step*h_step));
}
int UserInput::getIterations()
{
	return noTimeSteps;
}
//U_j-1 jeweils, endweder von den initioal conditions oder verhergehend berechnete werte
Matrix UserInput::calculate_rhs(int step, Matrix u)
{

	Matrix rhs(n_steps-2,1);
	if(step == 0)
	{
		for (int i = 1; i < n_steps-2; ++i)
		{
			rhs.setElement(i,1, (1/Tau_step)  * ( ((double)i * A_bound)/((double)(n_steps-1)) ) );
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

		double specialvalue = (1/Tau_step)  * u.getElement(n_steps-2,1) -  calculate_g() * A_bound * cos(omega_bound * Tau_step*(double)(step))   ;
				
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
