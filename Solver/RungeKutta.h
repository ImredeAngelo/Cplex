#pragma once
#include <vector>

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
		ButcherTable(v c, vv a, v b) : a(a), b(b), c(c), order(c.size()) {}

		//v operator[](int i) const { return a[i]; }
	};

	// Eulers method (first order)
	void euler(double (*f)(double, double), const double& t0, const double& t1, const double y0, const double h = 0.1);

	// Standard fourth order method
	void rk4(double (*f)(double, double), const double& t0, const double& t1, const double y0, const double h = 0.1);

	// Solves y'(x) = f(x, y(x)) using an explicit Butcher table
	void expl(ButcherTable butcher, double (*f)(double, double), const double& t0, const double& t1, const double y0, const double h = 0.1);

	// Transform to constrained optimization problem
	void discretize(double (*f)(double, double), const double& t0, const double& t1, const double y0, const double h = 0.1);
};

