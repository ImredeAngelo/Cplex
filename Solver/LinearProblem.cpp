#include "LinearProblem.h"

LinearProblem::LinearProblem(uint32_t steps, double t_end, double ulb, double uub, double ylb, double yub)
	: n(steps), t_end(t_end), model(env), u(env, n, ulb, uub), y(env, n, ylb, yub)
{
	// TODO - Objective function as argument?
	// Objective = Integral of y
	IloNumExprArg obj = u[n - 1];
	for (auto i = 0; i < n; i++)
		obj = obj + (y[i] * y[i] / 2);

	model.add(IloMinimize(env, obj));
	model.add(y[0] == 1);
}

// TODO: Butcher table
void LinearProblem::parameterize(double(*a)(double), double(*b)(double), double (*c)(double))
{
	const auto dt = t_end/n;
	for(auto i = 0; i < n - 1; i++) {
		auto t = i * dt;
		model.add(y[i + 1] == y[i] + dt * (a(t) * u[i] + b(t) * y[i] + c(t)));
	}
}

void LinearProblem::solve()
{
	try
	{
		cplex = IloCplex(model);
		cplex.solve();
	}
	catch (IloException& e) 
	{
		std::cerr << "Concert exception caught: " << e << std::endl;
	}
	catch (...) 
	{
		std::cerr << "An unknown error occured.";
	}
}

IloNumVarArray LinearProblem::getControl()
{
	return u;
}

IloNumVarArray LinearProblem::getConstraint()
{
	return y;
}

void LinearProblem::addConstraint(IloExtractable constraint)
{
	model.add(constraint);
}

std::ostream& operator<<(std::ostream& os, const LinearProblem& lp)
{
	for (auto i = 0; i < lp.n; i++)
		os << lp.cplex.getValue(lp.y[i]) << (i == (lp.n - 1) ? "\n" : ", ");
	return os;
}

std::ostream& operator<<=(std::ostream& os, const LinearProblem& lp)
{
	for (auto i = 0; i < lp.n; i++)
		os << lp.cplex.getValue(lp.u[i]) << (i == (lp.n - 1) ? "\n" : ", ");
	return os;
}