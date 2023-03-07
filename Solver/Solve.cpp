#include "Solve.h"
#include "RungeKutta.h"

// TODO: Algebraic constraints
std::vector<double> Solve::linear(double(*Fc)(double), double(*Fy)(double), double(*Fu)(double), double t0, double t1, double dt)
{
    int n = ((t1 - t0) / dt); // n as parameter instead of dt?

    std::vector<double> control(n);

    IloEnv env;
    IloModel model(env);
    IloNumVarArray u(env, n, 0.0, 1.0);
    IloNumVarArray y(env, n, 0.0, 10.0); // TODO: Proper constraints

    try {    
        // Constraint on each point u_i (Eulers method)
        for(auto i = 0; i < n - 1; i++) {
            auto t = i * dt + t0;
            model.add(y[i + 1] == y[i] + dt * (Fc(t) + Fy(t) * y[i] + Fu(t) * u[i]));
        }
        //RungeKutta::discretize(model, u, y);

        // TODO - Objective function as argument
        // Objective = Integral of y
        IloNumExprArg obj = u[n - 1];
        for (auto i = 0; i < n; i++) {
            obj = obj + (y[i]/*y[i]*/ / 2);
        }
        
        model.add(IloMinimize(env, obj));
        model.add(y[0] == 1);

        // Solve
        IloCplex cplex(model);
        cplex.solve();

        for (auto i = 0; i < n; i++) {
            control[i] = cplex.getValue(u[i]);
        }

#if _DEBUG
        std::cout << std::endl;
        std::cout << std::setw(14) << "Time";
        std::cout << std::setw(14) << "Objective";
        std::cout << std::setw(14) << "Control";
        std::cout << "\n      ----------------------------------------";
        std::cout << std::endl << std::setprecision(4);

        for (int i = 0; i < n; i++) {
            std::cout << std::setw(14) << dt*i + t0;
            std::cout << std::setw(14) << cplex.getValue(y[i]);
            std::cout << std::setw(14) << cplex.getValue(u[i]);
            std::cout << std::endl;
        }
#endif
    }
    catch (IloException& e) {
        std::cerr << "Concert exception caught: " << e << std::endl;
    }
    catch(...) {
        std::cerr << "An unknown error occured.";
    }

    return control;
}
