#include "UserInput.hpp"
#include <iostream>
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
	return (1 / Tau_step) + ((2*a_heat*a_heat)/h_step*h_step);
}
double UserInput::calculate_g()
{
	return -((a_heat*a_heat) / (h_step*h_step));
}
//~UserInput();
