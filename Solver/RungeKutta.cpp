#include "RungeKutta.h"
#include <iostream>
#include <iomanip>

// DEBUG
double y_debug(double x, double y) {
	return (10 * x - x * x / 2 - 1);
}

void RungeKutta::euler(double (*f)(double, double), const double& t0, const double& t1, const double y0, const double h)
{
	std::vector<double> y = { y0 };

	for(auto tn = t0, yn = y0; tn < t1; tn += h)
	{
		yn += h * f(tn, yn);
		y.push_back(yn);
	}

#if _DEBUG
	std::cout << "t" << "\t\t" << "eval" << "\t\t" << "true" << "\t\t" << "error" << "\n";
	std::cout << "-------------------------------------------------------\n";

	for (auto i = 0; i < y.size(); i++) {
		std::cout << std::setprecision(4) <<
			(t0 + i * h) << "\t\t" <<
			y[i] << "\t\t" <<
			y_debug(t0 + i * h, y[i]) << "\t\t" <<
			y[i] - y_debug(t0 + i * h, y[i]) << "\n";
	}
#endif
}

void RungeKutta::rk4(double (*f)(double, double), const double& t0, const double& t1, const double y0, const double h)
{
	std::vector<double> y = { y0 };
	
	for (auto tn = t0, yn = y0; tn < t1; tn += h)
	{
		double k1 = f(tn,		yn);
		double k2 = f(tn + h/2, yn + h/2 * k1);
		double k3 = f(tn + h/2, yn + h/2 * k2);
		double k4 = f(tn + h,	yn + h * k3);

		yn += h * (k1 + 2 * k2 + 2 * k3 + k4)/6;
		y.push_back(yn);
	}

#if _DEBUG
	double error = 0;

	std::cout << "t" << "\t\t" << "eval" << "\t\t" << "true" << "\t\t" << "error" << "\n";
	std::cout << "-------------------------------------------------------\n";

	for (auto i = 0; i < y.size(); i++) {
		std::cout << std::setprecision(4) <<
			(t0 + i * h) << "\t\t" <<
			y[i] << "\t\t" <<
			y_debug(t0 + i * h, y[i]) << "\t\t" <<
			y[i] - y_debug(t0 + i * h, y[i]) << "\n";

		auto estimate = y[i];
		auto actual = y_debug(t0 + i * h, y[i]);
		error += (estimate - actual) * (estimate > actual) + (actual - estimate) * (estimate < actual);
	}

	std::cout << "\n\nTotal Error: " << error << "\n\n";
#endif
}

void RungeKutta::expl(ButcherTable butcher, double (*f)(double, double), const double& t0, const double& t1, const double y0, const double h)
{
	std::vector<double> y = { y0 };

	for (double tn = t0, yn = y0, bk = 0; tn < t1; tn += h, bk = 0)
	{
		std::vector<double> k;

		// Get all k_i
		for (auto i = 0; i < butcher.order; i++)
		{
			double ti = tn + butcher.c[i] * h;
			double yi = yn;

			// Explicit only (j < i)
			for (auto j = 0; j < i; j++)
				yi += butcher.a[i][j] * k[j];

			k.emplace_back(h * f(ti, yi));
		}

		// Get next y_i
		for (auto i = 0; i < butcher.order; i++)
			bk += butcher.b[i] * k[i];

		yn += bk;
		y.emplace_back(yn);
	}

#if _DEBUG
	double error = 0;

	std::cout << "t" << "\t\t" << "eval" << "\t\t" << "true" << "\t\t" << "error" << "\n";
	std::cout << "-------------------------------------------------------\n";

	for (auto i = 0; i < y.size(); i++) {
		std::cout << std::setprecision(4) <<
			(t0 + i * h) << "\t\t" <<
			y[i] << "\t\t" <<
			y_debug(t0 + i * h, y[i]) << "\t\t" <<
			y[i] - y_debug(t0 + i * h, y[i]) << "\n";

		auto estimate = y[i];
		auto actual = y_debug(t0 + i * h, y[i]);
		error += (estimate - actual) * (estimate > actual) + (actual - estimate) * (estimate < actual);
	}

	std::cout << "\n\nTotal Error: " << error << "\n\n";
#endif
}

void RungeKutta::discretize(IloModel model, IloNumVarArray& u, IloNumVarArray& y, ButcherTable butcher)
{
	
}