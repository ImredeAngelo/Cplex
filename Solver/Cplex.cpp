#include "Cplex.h"
#include <ilcplex/ilocplex.h>
#include <vector>

void Cplex::test()
{
	using namespace std;
	
	IloEnv env;
	try {
		IloModel model(env);
		IloNumVarArray vars(env);
		vars.add(IloNumVar(env, 0.0, 40.0));
		vars.add(IloNumVar(env));
		vars.add(IloNumVar(env));
		model.add(IloMaximize(env, vars[0] + 2 * vars[1] + 3 * vars[2]));
		model.add(-vars[0] + vars[1] + vars[2] <= 20);
		model.add(vars[0] - 3 * vars[1] + vars[2] <= 30);
		
		IloCplex cplex(model);
		if (!cplex.solve()) {
			env.error() << "Failed to optimize LP." << endl;
			throw(-1);
		}
		IloNumArray vals(env);
		env.out() << "Solution status = " << cplex.getStatus() << endl;
		env.out() << "Solution value = " << cplex.getObjValue() << endl;
		cplex.getValues(vals, vars);
		env.out() << "Values = " << vals << endl;
	}
	catch (IloException& e) {
		cerr << "Concert exception caught: " << e << endl;
	}
	catch (...) {
		cerr << "Unknown exception caught" << endl;
	}
}

//
//   minimize y*y/2
//   subject to  Ax = b
//               l <= x <= u
//   where
//
//   A = ( [1 + h*fy]		     -1				 0		h * fu			 0			 0 )	b = ( -h * f0 )
//       (		    0	 [1 + h*fy]				-1			 0		h * fu			 0 )		( -h * f1 )
//       (		    0  		      0		[1 + h*fy]	   	 	 0			 0		h * fu )		( -h * f2 )
// 
//   l = (  0  0  0  0  0  0 )
//   u = (  n  n  n  1  1  1 ) -> n = infinity (TODO)
//
//
// x = (y0, y1, y2, u0, u1, u2)
// 
void Cplex::solveLP()
{
	using namespace std;
	cout << "\n\nSolving LP:\n\n";
	
	IloEnv env;
	IloModel model(env);

	const auto dt = 0.25;
	const auto Fo = 0.0, Fy = 0.0, Fu = -0.5;
	const auto t0 = 0.0, t1 = 3.0;
	
	int n = ((t1 - t0) / dt) + 1;
	int i = 0;

	cout << "N: " << n << "\n";

	try {
		IloNumVarArray x(env, 2*n, 0.0, 1.0);

		model.add(IloMinimize(env, x[0]*x[0]/2 + x[1]*x[1]/2 + x[2]*x[2]/2));

		// Initial condition y(0) = 1
		model.add(x[0] == 1);

		//// Constraint on each point u_i
		//model.add(x[1] == x[0] + dt * (Fo + Fy * x[0] + Fu * x[3]));
		//model.add(x[2] == x[1] + dt * (Fo + Fy * x[1] + Fu * x[4]));

		for (auto t = t0; t < t1; t += dt) {
			model.add(x[i + 1] == x[i] + dt * (Fo + Fy * x[i] + Fu * x[n + i]));
			i++;
		}

		IloCplex cplex(model);
		//cplex.setParam(IloCplex::Param::Simplex::Display, 2);
		//cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Network);
		cplex.solve();
		cplex.out() << "\n\nAfter optimization, values are:" << endl;
			//<< "(" << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3] << ", " << x[4] << ", " << x[5] << ")" << endl;

		for (int i = 0; i < 2*n - 1; i++) {
			cout << cplex.getValue(x[i]) << ", ";
			//cout << std::round(cplex.getValue(x[i])) << ", ";
		}
	}
	catch (IloException& e) {
		cerr << "Concert exception caught: " << e << endl;
	}
	catch (...) {
		cerr << "Unknown exception caught" << endl;
	}

	env.end();
	
	std::cout << "\n\nComplete\n\n";
}

//// Set up system to be solved
//IloNumVar u(env, 0.0, 1.0, ILOFLOAT);	// Control
//IloNumVar y(env);						// Dynamics
//// IloNumLinExprTerm dy = -u;			// Dynamics PDE -> RK to get constraints on y

//model.add(IloMinimize(env, y * y / 2));

//IloRange c1(y >= 0);
//model.add(c1);

// solve