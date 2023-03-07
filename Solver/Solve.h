#pragma once
#include <vector>
#include <ilcplex/ilocplex.h>

namespace Solve
{
	/// <summary>
	/// Solves problem of the form f = Fc(t) + Fy(t)*y + Fu(t)*u
	/// </summary>
	/// <param name="Fc">Constant coefficient</param>
	/// <param name="Fy">Coefficient</param>
	/// <param name="Fu">Coefficient</param>
	/// <param name="t0">RK start time</param>
	/// <param name="t1">RK end time</param>
	/// <param name="dt">Time step</param>
	/// <returns>Discrete approximation of the control function u(t)</returns>
	std::vector<double> linear(double(*Fc)(double), double(*Fy)(double), double(*Fu)(double), double t0, double t1, double dt = 0.1);
};

