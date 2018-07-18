#include "UserInput.hpp"
#include <iostream>
#include <math.h>

UserInput::UserInput()
{
	userDialog_init();
}

void UserInput::userDialog_init()
{
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

	h_step = L_length / (float)n_steps;
	std::cout<<"\nStep size is: "<< h_step <<"\nTime step Tau: ";
	std::cin>>Tau_step;
	std::cout<<"\n\nBoundary condition.\nSet A: ";
	std::cin>>A_bound;
	std::cout<<"\nThanks. Set omega: ";
	std::cin>>omega_bound;
	std::cout<<"\n\n\nComplete.\n\n\n";
}

double UserInput::calculate_f()
{
	return (1 / Tau_step) +/*war mal +*/ ((2*a_heat*a_heat)/h_step*h_step);
}
double UserInput::calculate_g()
{
	return -((a_heat*a_heat) / (h_step*h_step));/*war mal -*/
}
Matrix UserInput::calculate_rhs()
{
	Matrix rhs(n_steps-2,1);//1, 2, 3, ..., n-2, x

	for (int i = 1; i < n_steps-2; ++i)
	{
		rhs.setElement(i,1, (1/Tau_step)  * ( ((float)i * A_bound)/((float)(n_steps-1)) ) );
	}
	double specialvalue = (1/Tau_step)  *  (((n_steps-2)*A_bound)  /  (n_steps-1));
			specialvalue = specialvalue - calculate_g() * A_bound * cos(omega_bound * Tau_step);
 	rhs.setElement(n_steps-2,1,   specialvalue);

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
