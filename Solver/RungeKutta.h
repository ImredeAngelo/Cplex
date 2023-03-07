#pragma once
#include <vector>
#include <ilcplex/ilocplex.h>

#define Euler	{{0},{{0}},{1}}
#define Heun	{{0,1},{{0,0},{1,0}},{1 /(double)2, 1 /(double)2}}
#define RK4		{{0,0.5,0.5,1}, {{ 0, 0, 0, 0},{.5, 0, 0, 0},{ 0,.5, 0, 0},{ 0, 0, 1, 0},},{(1 / (double)6), (1 / (double)3), (1 / (double)3), (1 / (double)6)}}

namespace RungeKutta
{
	typedef std::vector<double> v;
	typedef std::vector<v> vv;

	struct ButcherTable
	{
		const size_t order;

		const vv a;
		const v b;
		const v c;

		// TODO: Error handling -> dimensions
		// TODO: Construct from matrix
		ButcherTable(v c, vv a, v b) : a(a), b(b), c(c), order(c.size()) {}

		//v operator[](int i) const { return a[i]; }
	};

	// Eulers method (first order)
	void euler(double (*f)(double, double), const double& t0, const double& t1, const double y0, const double h = 0.1);

	// Standard fourth order method
	void rk4(double (*f)(double, double), const double& t0, const double& t1, const double y0, const double h = 0.1);

	// Solves y'(x) = f(x, y(x)) using an explicit Butcher table
	void expl(ButcherTable butcher, double (*f)(double, double), const double& t0, const double& t1, const double y0, const double h = 0.1);
	
	// Discretizes linear function and converts to CPLEX variables
	void discretize(IloModel model, IloNumVarArray& u, IloNumVarArray& y, ButcherTable butcher = RK4);

	/*
	template <typename T, typename... Args>
	void discretize(uint32_t n, double t0, double t1, T y, T u) {};
	*/
};
