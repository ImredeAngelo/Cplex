#pragma once
#include <ilcplex/ilocplex.h>
#include <vector>
#include <iostream>

/// <summary>
/// min F = integral(y) s.t.
/// <para> y'(t) = a*u(t) + b*y(t) + c</para> 
/// <para> lb &#x2264; u(t) &#x2264; ub</para> 
/// </summary>
class LinearProblem
{
public:
	LinearProblem(uint32_t steps, double time, double ulb, double uub, double ylb = DBL_MIN, double yub = DBL_MAX);

	/// <summary>
	/// Complete parameterization of y' = a(t)*u(t) + b(t)*y(t)
	/// </summary>
	/// <param name="a"></param>
	/// <param name="b"></param>
	void parameterize(double (*a)(double), double (*b)(double), double (*c)(double));
	/// <summary>
	/// Solve the problem
	/// </summary>
	void solve();

	// For defining custom constraints
	IloNumVarArray getControl();
	IloNumVarArray getConstraint();
	void addConstraint(IloExtractable constraint);

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
