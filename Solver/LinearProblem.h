#pragma once
#include <ilcplex/ilocplex.h>
#include <vector>
#include <iostream>

/// <summary>
/// min F = integral(y)
/// <para>s.t. 0 = c*u(t) + d*y(t)</para> 
/// <para> y'(t) = a*u(t) + b*y(t)</para> 
/// <para> lb &#x2264; u(t) &#x2264; ub</para> 
/// </summary>
class LinearProblem
{
public:
	LinearProblem(uint32_t steps, double time, double lb, double ub);

	/// <summary>
	/// Complete parameterization of y'
	/// </summary>
	/// <param name="a"></param>
	/// <param name="b"></param>
	void parameterize(double (*a)(double), double (*b)(double));
	/// <summary>
	/// Solve the problem
	/// </summary>
	void solve();

	/*
	IloNumVarArray getU();
	IloNumVarArray getY();
	void addConstraint();
	*/

private:
	double t_end;
	uint32_t n;

	IloEnv env;
	IloModel model;
	IloCplex cplex;

	IloNumVarArray u;
	IloNumVarArray y;

public:
	friend std::ostream& operator<< (std::ostream& os, const LinearProblem& lp);
	friend std::ostream& operator<<= (std::ostream& os, const LinearProblem& lp);
};
