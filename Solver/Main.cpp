#include <iostream>
#include "Cplex.h"
#include "RungeKutta.h"
#include "Solve.h"
#include "LinearProblem.h"

double f(double x, double y) {
	return (10 - x);
}

int main()
{	
	//RungeKutta::rk4(f, -5, 5, -63.5, 0.25);
	//RungeKutta::expl(RK4, f, -5, 5, -63.5, 0.25);

	//Cplex::solveLP();

	auto Fc = [](double t) { return 0.0; };
	auto Fy = [](double t) { return 0.4; };
	auto Fu = [](double t) { return -1.0; };
	// auto u = Solve::linear(Fc, Fy, Fu, 0, 2.1, 0.01);

	LinearProblem lin(40, 4, 0.0, 1.0);
	lin.parameterize(Fu, Fy);
	lin.solve();

	std::cout << "\n\nObjective:\n" << lin;
	std::cout << "\nControl:\n" <<= lin;

	return 0;
}