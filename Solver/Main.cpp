#include <iostream>
#include "Cplex.h"
#include "RungeKutta.h"

double f(double x, double y) {
	return (10 - x);
}

int main()
{
	using RungeKutta::ButcherTable;
	using namespace std;
	
	//RungeKutta::rk4(f, -5, 5, -63.5, 0.25);
	//cout << "\n\nButcher Tablau:\n\n";

	//// TODO: Eulers method has large error when using Butcher table
	//ButcherTable Euler{ {0},{{0}},{1} };
	//ButcherTable Heun{ {0,1},{{0,0},{1,0}},{1 /(double)2, 1 /(double)2} };
	//ButcherTable RK4{
	//	{0, 0.5, 0.5, 1},
	//	{
	//		{ 0, 0, 0, 0},
	//		{.5, 0, 0, 0},
	//		{ 0,.5, 0, 0},
	//		{ 0, 0, 1, 0},
	//	},
	//	{(1 / (double)6), (1 / (double)3), (1 / (double)3), (1 / (double)6)}
	//};

	//RungeKutta::expl(RK4, f, -5, 5, -63.5, 0.25);

	//Cplex::test();
	Cplex::solveLP();
	return 0;
}