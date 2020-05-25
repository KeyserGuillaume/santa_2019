#include "dwlb.h"
#include "Presets.h"
#include "day_sub_problem.h"

ILOSTLBEGIN

void test_stuff(Presets &presets, const std::vector<unsigned int> &initial_solution, const double& cost) {
//    std::vector<std::vector<char>> selected_choices(0);
//    for (unsigned int i = 0; i < NB_FAMILIES; i++)
//        for (unsigned int k = 0; k < K_MAX; k++)
//            if (initial_solution[i] == presets.get_family_data(i, k))
//                selected_choices.push_back(std::vector<char>(1, k));
//    Column col(selected_choices, 10*cost);
//    std::vector<Column> columns(1, col);

    std::vector<Column> columns(0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        columns.push_back(Column(i));

    IloEnv env;

    // define useful array of one
    IloNumArray ones(env, NB_FAMILIES);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) ones[i] = 1;

    IloModel model(env);
    IloCplex cplex(model);

    // objective
    IloObjective objective = IloAdd(model, IloMinimize(env));
    // constraints: first, each family has at most one choice
    IloRangeArray at_most_once = IloAdd(model, IloRangeArray(env, -IloInfinity, ones));
    IloRange sum_to_one = IloAdd(model, IloRange(env, 1, 1));

    // define the variables called lambda in the document "Decomposition de Dantzig-Wolfe"
    IloNumVarArray cols(env);

    // add whatever columns previously defined
    for (unsigned int l = 0; l < columns.size(); l++)
        cols.add(IloNumVar(objective(columns[l].get_cost() + LARGE_UPPER_BOUND * (NB_FAMILIES - columns[l].get_sum_family_use_nb())) +
                           at_most_once(columns[l].get_family_use_nb(env)) +
                           sum_to_one(1), 0));

    // add as many columns as necessary in order to make the constraints matrix a full rank matrix
    /*Presets tmp_presets = presets;
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        tmp_presets.assign_family(i, initial_solution[i], false, false);
    tmp_presets.compute_feasibility();
    tmp_presets.compute_all_bounds();
    std::vector<unsigned int> occupancies(NB_DAYS, 0);
    for (unsigned int j = 0; j < NB_DAYS; j++)
        occupancies[j] = tmp_presets.get_presets_occupancy(j);
    double prev_day_cost = get_day_cost(occupancies);
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        for (unsigned int k = 0; k < K_MAX; k++){
            if (selected_choices[i][0] == k) continue;
            unsigned int j = presets.get_family_data(i, k);
            if (tmp_presets.get_presets_occupancy(j) + presets.get_family_size(i) <= MAX_NB_PEOPLE_PER_DAY){
                selected_choices[i].push_back(k);
                occupancies[j] += presets.get_family_size(i);
                double new_cost = columns[0].get_cost() + get_day_cost(occupancies) - prev_day_cost +
                        CONSTANT_COST[k] + MARGINAL_COST[k] * presets.get_family_size(i);
                columns.push_back(Column(selected_choices, new_cost));
                cols.add(IloNumVar(objective(columns[i + 1].get_cost()) +
                                   at_most_once(columns[i + 1].get_family_use_nb(env)) +
                                   sum_to_one(1)));
                selected_choices[i].pop_back();
                occupancies[j] -= presets.get_family_size(i);
                break;
            }
            if (k == K_MAX - 1){
                char k_orig = selected_choices[i][0];
                j = presets.get_family_data(i, k_orig);
                if (tmp_presets.get_presets_occupancy(j) - presets.get_family_size(i) >= MIN_NB_PEOPLE_PER_DAY){
                    selected_choices[i].pop_back();
                    occupancies[j] -= presets.get_family_size(i);
                    double new_cost = columns[0].get_cost() + get_day_cost(occupancies) - prev_day_cost +
                                      CONSTANT_COST[k] + MARGINAL_COST[k] * presets.get_family_size(i);
                    columns.push_back(Column(selected_choices, new_cost));
                    cols.add(IloNumVar(objective(columns[i + 1].get_cost()) +
                                       at_most_once(columns[i + 1].get_family_use_nb(env)) +
                                       sum_to_one(1)));
                    selected_choices[i].push_back(k_orig);
                    occupancies[j] += presets.get_family_size(i);
                } else {
                    throw std::logic_error("I wish this could not happen");
                }
            }
        }
    }*/

    for (;;) {

        unsigned int n = columns.size();

        // solve
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

//        std::vector<double> lam(n, 0);
//        double sum_lam = 0;
//        for (unsigned int m = 0; m < n; m++) {
//            lam[m] = cplex.getValue(cols[m]);
//            sum_lam += lam[m];
//        }

        // extract dual stuff
        std::vector<double> duals;
        for (unsigned int i = 0; i < NB_FAMILIES; i++) {
            duals.push_back(-cplex.getDual(at_most_once[i]));
            if (duals[i] != 0) std::cout << "family " << i << " will be given a lambda equal to " << duals[i] << std::endl;
        }

        double eta = cplex.getDual(sum_to_one);

        // solve the sub-problems
        RLLB rllb(presets, duals);

        // check whether the minimum reduced cost is negative
        double reduced_cost = rllb.get_lb() - eta;
        if (reduced_cost >= 0)
            break;

        columns.push_back(Column(rllb.get_selected_choices(presets), rllb.get_cost(presets)));
        std::cout << "adding one column to the model with a cost of " << columns[n].get_cost() << " and a sum of family uses of " << columns[n].get_sum_family_use_nb() << std::endl;
        cols.add(IloNumVar(objective(columns[n].get_cost() + LARGE_UPPER_BOUND * (NB_FAMILIES - columns[n].get_sum_family_use_nb()))  +
                           at_most_once(columns[n].get_family_use_nb(env)) +
                           sum_to_one(1)));
    }

}

void tiny_test() {
    IloEnv env;
    IloModel model(env);

    // variables
    IloNumVar x(env);
    IloNumVar y(env);

    // objective
    model.add(IloMinimize(env, x + y));

    // constraints
    IloRange c1(env, 2*x + y);
    c1.setLB(1);
    model.add(c1);
    IloRange c2(env, x + 2*y);
    model.add(c2);
    c2.setLB(1);
    IloRange c3(env, 3*x + 3*y);
    model.add(c3);
    c3.setLB(2);

    IloCplex cplex(model);
    cplex.solve();
    std::cout << "solve status: " << cplex.getStatus() << std::endl;
    std::cout << "Objective value: " << cplex.getObjValue() << std::endl;
    std::cout << cplex.getDual(c1) << std::endl;
    std::cout << cplex.getDual(c2) << std::endl;
    std::cout << cplex.getDual(c3) << std::endl;
}

//         variables
//        unsigned int n = columns.size();
//        IloNumVarArray lambda(env, n, 0, 1);
//
//        // objective
//        IloExpr obj(env);
//        for (unsigned int m = 0; m < n; m++)
//            obj += columns[m].get_cost() * lambda[m];
//        model.add(IloMinimize(env, obj));
//        obj.end();
//
//        // constraints
//        IloExpr exp(env);
//        for (unsigned int m = 0; m < n; m++)
//            exp += lambda[m];
//        IloRange convexity(env, exp);
//        convexity.setBounds(1., 1.);
//        model.add(convexity);
//        exp.end();
//
//        IloRangeArray family_constraints(env, NB_FAMILIES);
//        for (unsigned int i = 0; i < NB_FAMILIES; i++) {
//            IloExpr exp(env);
//            for (unsigned int m = 0; m < n; m++)
//                exp += lambda[m] * int(columns[m].get_family_use_nb(i));
//            family_constraints[i] = (1. <= exp <= 1.);
//            exp.end();
//        }
//        model.add(family_constraints);
//        IloCplex cplex(model);