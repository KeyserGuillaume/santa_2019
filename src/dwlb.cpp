#include "dwlb.h"
#include "Presets.h"
#include "day_sub_problem.h"

ILOSTLBEGIN

void test_stuff(Presets &presets, const std::vector<unsigned int> &initial_solution, const double& cost) {
    std::vector<std::vector<unsigned int>> selected_choices(0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        for (unsigned int k = 0; k < K_MAX; k++)
            if (initial_solution[i] == presets.get_family_data(i, k))
                selected_choices.push_back(std::vector<unsigned int>(1, k));
    Column col(selected_choices, cost);
    std::vector<Column> columns(1, col);

    IloEnv env;
    IloModel model;

    while(true) {
        env = IloEnv();
        model = IloModel(env);

        // variables
        unsigned int n = columns.size();
        IloNumVarArray lambda(env, n, 0, 1);

        // objective
        IloExpr obj(env);
        for (unsigned int m = 0; m < n; m++)
            obj += columns[m].get_cost() * lambda[m];
        model.add(IloMinimize(env, obj));
        obj.end();

        // constraints
        IloExpr exp(env);
        for (unsigned int m = 0; m < n; m++)
            exp += lambda[m];
        IloRange convexity(env, exp);
        convexity.setBounds(1., 1.);
        model.add(convexity);
        exp.end();

        IloRangeArray family_constraints(env, NB_FAMILIES);
        for (unsigned int i = 0; i < NB_FAMILIES; i++) {
            IloExpr exp(env);
            for (unsigned int m = 0; m < n; m++)
                exp += lambda[m] * int(columns[m].get_family_use_nb(i));
            family_constraints[i] = (1. <= exp <= 1.);
            exp.end();
        }
        model.add(family_constraints);

        // solve
        IloCplex cplex(model);
        std::cout << std::endl;
        std::cout << "solving with Cplex. Cplex says:" << std::endl;
        std::cout << "===================================================================================="
                  << std::endl;
//    cplex.setOut(env.getNullStream());
        cplex.solve();
        std::cout << "===================================================================================="
                  << std::endl;
        std::cout << "solve status: " << cplex.getStatus() << std::endl;
        std::cout << "Objective value: " << cplex.getObjValue() << std::endl;

        std::vector<double> lam(n, 0);
        double sum_lam = 0;
        for (unsigned int m = 0; m < n; m++) {
            lam[m] = cplex.getValue(lambda[m]);
            sum_lam += lam[m];
        }

        // extract dual stuff
        std::vector<double> duals;
        for (unsigned int i = 0; i < NB_FAMILIES; i++) {
            duals.push_back(-cplex.getDual(family_constraints[i]));
            if (duals[i] != 0) std::cout << "family " << i << " will be given a lambda equal to " << duals[i] << std::endl;
        }

        double eta = cplex.getDual(convexity);

        cplex.end();
        family_constraints.end();
        convexity.end();
        model.end();
        env.end();

        // solve the sub-problems
        RLLB rllb(presets, duals);

        // check whether the minimum reduced cost is negative
        double reduced_cost = rllb.get_lb() - eta;
        if (reduced_cost >= 0)
            break;

        columns.push_back(Column(rllb.get_selected_choices(presets), rllb.get_lb()));
    }

}
