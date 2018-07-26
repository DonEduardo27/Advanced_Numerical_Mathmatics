#include "Matrix.hpp"

class UserInput
{
public:
	UserInput();

	void userDialog_init(bool user);

	double calculate_f();
	double calculate_g();
	int getIterations();
	Matrix calculate_rhs(int step, Matrix u);
	Matrix make_matrix();
private:
double a_heat;

double L_length;
double Tau_step;
double h_step;
int    n_steps;
int    noTimeSteps;

double A_bound;
double omega_bound;
	
double f_Matrix;
double g_Matrix;
};