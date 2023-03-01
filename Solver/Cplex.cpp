#include "Cplex.h"
#include <ilcplex/ilocplex.h>
#include <vector>

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
	const auto Fo = 0.0, Fy = 0.0, Fu = -1.0;
	const auto t0 = 0.0, t1 = 3.0;
	
	int n = ((t1 - t0) / dt) + 1;
	int i = 0;

	cout << "N: " << n << "\n";

	try {
		IloNumVarArray x(env, 2*n, 0.0, 1.0);

		model.add(IloMinimize(env, x[0]*x[0]/2 + x[1]*x[1]/2 + x[2]*x[2]/2 + x[2*n - 1]));

		// Initial condition y(0) = 1
		model.add(x[0] == 1);

		// Constraint on each point u_i (Eulers method)
		for (auto t = t0; t < t1; t += dt, i++) {
			model.add(x[i + 1] == x[i] + dt * (Fo + Fy * x[i] + Fu * x[n + i]));
		}

		IloCplex cplex(model);
		cplex.solve();

		// Output
		// cplex.out() << "\n\nAfter optimization, values are:" << endl;

		cout << "\nObjective: " << setw(2);
		for (int i = 0; i < n; i++) {
			cout << cplex.getValue(x[i]) << ", " << setw(5);
		}
		cout << "\nControl: " << setw(4);
		for (int i = 0; i < n; i++) {
			cout << cplex.getValue(x[i + n]) << ", " << setw(5);
		}

		cout << endl;
	}
	catch (IloException& e) {
		cerr << "Concert exception caught: " << e << endl;
	}
	catch (...) {
		cerr << "Unknown exception caught" << endl;
	}

	env.end();
}

//// Set up system to be solved
//IloNumVar u(env, 0.0, 1.0, ILOFLOAT);	// Control
//IloNumVar y(env);						// Dynamics
//// IloNumLinExprTerm dy = -u;			// Dynamics PDE -> RK to get constraints on y

//model.add(IloMinimize(env, y * y / 2));

//IloRange c1(y >= 0);
//model.add(c1);

// solve