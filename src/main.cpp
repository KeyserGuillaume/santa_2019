#include "tools.h"
#include "Flow.h"
#include "LocalSearch.h"
#include "Presets.h"
#include "day_sub_problem.h"
#include "BranchAndBound.h"
#include "dwlb.h"

void compare_solutions(const Assignment &A, const Assignment &B) {
    for (unsigned int i = 0; i < NB_DAYS; i++)
        if (A.get_ith_day(i)->get_N() != B.get_ith_day(i)->get_N())
            std::cout << "Day " << i << ": " << A.get_ith_day(i)->get_N() << " -> " << B.get_ith_day(i)->get_N() << std::endl;
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        unsigned int k_A = A.get_ith_family(i)->get_k();
        unsigned int k_B = B.get_ith_family(i)->get_k();
        unsigned int day_A = A.get_ith_family(i)->get_assigned_day()->get_id();
        unsigned int day_B = B.get_ith_family(i)->get_assigned_day()->get_id();
        if (k_B != k_A) {
            std::cout << "Family " << i << " (" << A.get_ith_family(i)->get_nb_people() << "): ";
            std::cout << day_A << " -> " << day_B << " (" << k_A << " -> " << k_B << ")" << std::endl;
        }
    }
}

void try_local_search(Assignment &A) {
    LocalSearch LS(&A);
    LS.run(10000000, 1000000);
    //LS.run(1000, 100);
    //LS.run_on_time_limit(begin + 60*60*CLOCKS_PER_SEC, 1000000);

    std::cout << A.get_cost() << std::endl;
    A.stats();
    LS.stats();

    A.write_solution("./solutions/alocal_search_solution_" + std::to_string(A.get_cost()) + "_.csv");
}

std::vector<bool> choose_families_2(const Assignment &A,const  Presets &presets, const unsigned int& n1, const unsigned int& n2, const unsigned int& n3) {
    // choose n1 days
    std::vector<bool> focused_days = random_day_selection(n1);
    // choose n2 families on the focused days which would rather be on another focused day
    std::vector<unsigned int> coveting_families(0);
    std::vector<bool> res(NB_FAMILIES, false);
    for (unsigned int i = 0; i < NB_DAYS; i++){
        if (!focused_days[i]) continue;
        Day *current_day = A.get_ith_day(i);
        unsigned int n_families = current_day->get_nb_families();
        for (unsigned int j = 0; j < n_families; j++) {
            Family *current_family = current_day->get_ith_family(j);
            for (unsigned int k = 0; k < current_family->get_k(); k++)
                if (focused_days[current_family->get_ith_preferred_day(k)->get_id()])
                    coveting_families.push_back(current_family->get_id());
        }
    }
    std::vector<bool> coveting_selection = random_selection(n2, coveting_families.size());
    for (unsigned int m = 0; m < coveting_selection.size(); m++)
        if (coveting_selection[m])
            res[coveting_families[m]] = true;
    // choose the n3 - n2 families left among those verifying one criterion : interest in the focused days
    std::vector<unsigned int> families_interested(0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        if (res[i]) continue;
        for (unsigned int k = 0; k < K_MAX; k++)
            if (focused_days[presets.get_family_data(i, k)])
                families_interested.push_back(i);
    }
    families_interested = get_unique_vector(families_interested);
    std::vector<bool> family_selection = random_selection(n3 - n2, families_interested.size());
    for (unsigned int m = 0; m < families_interested.size(); m++)
        if (family_selection[m])
            res[families_interested[m]] = true;
    return res;
}

std::vector<bool> choose_families(const Assignment &A,const  Presets &presets, const unsigned int& n1, const unsigned int& n2) {
    std::vector<bool> res(NB_FAMILIES, false);
    std::vector<bool> day_possesses_chosen_family(NB_DAYS, false);

    // choose n2 families of rank 3 and all families of rank 4
    std::vector<unsigned int> rank_3(0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        unsigned int k = A.get_ith_family(i)->get_k();
        if (k == 4) {
            res[i] = true;
            day_possesses_chosen_family[A.get_ith_family(i)->get_assigned_day()->get_id()] = true;
        }
        if (k > 0)
            rank_3.push_back(i);
    }
    std::vector<bool> rank_3_selection = random_selection(n2, rank_3.size());
    for (unsigned int m = 0; m < rank_3.size(); m++) {
        if (rank_3_selection[m]) {
            res[rank_3[m]] = true;
            day_possesses_chosen_family[A.get_ith_family(rank_3[m])->get_assigned_day()->get_id()] = true;
        }
    }

    // choose n1 - n2 families among families which are interested in one of the days with a previously chosen day
    // assigned to it
    std::vector<unsigned int> families_interested(0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        if (res[i]) continue;
        for (unsigned int k = 0; k < K_MAX; k++)
            if (day_possesses_chosen_family[presets.get_family_data(i, k)])
                families_interested.push_back(i);
    }
    families_interested = get_unique_vector(families_interested);
    std::vector<bool> family_selection = random_selection(n1 - n2, families_interested.size());
    for (unsigned int m = 0; m < families_interested.size(); m++)
        if (family_selection[m])
            res[families_interested[m]] = true;
    return res;
}

// takes 50 families, finds the optimal solution when all other families are set, then repeat
void farm(const Assignment &A, Presets &presets, const std::vector<bool> &useful_families, const unsigned int &nb_times) {
    srand((unsigned int) time(0));
    for (unsigned int t = 0; t < nb_times; t++) {
        std::vector<bool> family_selection = choose_families_2(A, presets, 5, 25, 50);
        //std::vector<bool> family_selection = choose_families(A, presets, 15, 10);
//        std::vector<bool> day_selection = random_day_selection(10);
        for (unsigned int i = 0; i < NB_FAMILIES; i++) {
            unsigned int k = A.get_ith_family(i)->get_k();
            if (!family_selection[i] && k <= 3 && !useful_families[i])
            //if (i != 740 && i != 805 && i != 433 && i != 975 && i != 718 && i != 2373)
                presets.assign_family(i, k, false);
        }
        for (unsigned int i = 0; i < NB_DAYS; i++) {
            if (A.get_ith_day(i)->is_Friday_like())
                presets.prescribe_occupancy_ub(i, MIN_NB_PEOPLE_PER_DAY);
        }
        for (unsigned int i = 0; i < NB_DAYS; i++) {
            presets.prescribe_occupancy_lb(i, A.get_ith_day(i)->get_N() - 15);
            presets.prescribe_occupancy_ub(i, A.get_ith_day(i)->get_N() + 15);
        }
        presets.compute_all_bounds();
        presets.compute_feasibility();

        std::vector<unsigned int> reference_costs(NB_DAYS, 0);
        for (unsigned int i = 0; i < NB_DAYS; i++)
            reference_costs[i] = A.get_ith_day(i)->get_assignments_costs();// + ((day_selection[i]) ? 272 : 0);

        for (int i = NB_DAYS - 1; i >= 0; i--) {
            std::vector<unsigned int> assignations(0), counter_assignations(0);
            if (!anticipate_on_day(presets, i, reference_costs[i], assignations, counter_assignations))
                throw std::logic_error("how ?");
            for (unsigned int &family_index: assignations) {
                presets.assign_family(family_index, i, false, false);
            }
            for (unsigned int &family_index: counter_assignations) {
                presets.forbid_assignment(family_index, i, false, false);
            }
            if (assignations.size() > 0 || counter_assignations.size() > 0) {
                 presets.compute_all_bounds();
                presets.compute_feasibility();
            }
        }

        Graph G = Graph(presets);
        G.compute_max_flow_min_cost();
        std::cout << "flow is " << G.get_current_flow() << " with current cost " << G.get_flow_cost()
                  << " and with true cost " << G.get_true_flow_cost() << std::endl;
        G.check_flow();
        G.show_schedule();

        unsigned int nb_nodes = 0;
        clock_t t0 = clock();
        clock_t time_graph = 0;
        clock_t time_bounds = 0;
        clock_t time_anticipating = 0;
        std::cout << "Best solution seen is "
                  << brute_force(presets, nb_nodes, time_graph, time_bounds, time_anticipating,
                                 reference_costs, G.get_distribution()) << std::endl;
        std::cout << "Number of nodes is " << nb_nodes << std::endl;
        std::cout << "Time spent computing flows is  " << time_graph << std::endl;
        std::cout << "Time spent computing bounds is " << time_bounds << std::endl;
        std::cout << "Time spent anticipating is     " << time_anticipating << std::endl;
        std::cout << clock() - t0 << std::endl;

        for (unsigned int i = 0; i < NB_FAMILIES; i++)
            if (presets.is_family_alr_assigned(i))
                presets.deassign_family(i ,false);
    }
}


void try_anticipating_stuff(const Assignment &A, Presets &presets) {
    for (int i = NB_DAYS - 1; i >= 1; i--) {
        Day *day = A.get_ith_day(i);
        unsigned int cost = day->get_assignments_costs();
        std::cout << "We try assigning to day " << i << " of cost " << cost << std::endl;
        std::vector<unsigned int> assignations(0), counter_assignations(0);
        anticipate_on_day(presets, i, cost, assignations, counter_assignations);
        for (unsigned int& family_index: assignations){
            presets.assign_family(family_index, i, false, false);
        }
        for (unsigned int& family_index: counter_assignations){
            presets.forbid_assignment(family_index, i, false, false);
        }
        presets.compute_all_bounds(); presets.compute_feasibility();
    }
}

int main() {
    clock_t begin = clock();

    std::vector<unsigned int> initial_solution = read_solution("./solutions/rllb_lb_69153.773438");
//    std::vector<unsigned int> initial_solution = read_solution("./solutions/rllb_lb_71223.734375");
//    std::vector<unsigned int> initial_solution = read_solution("./solutions/flow_solution_71441_8150286.csv");
    //std::vector<unsigned int> initial_solution = read_solution("./solutions/local_search_solution_71480_.csv");
    //std::vector<unsigned int> initial_solution = read_solution("./solutions/alocal_search_solution_74253_.csv");
    std::vector<std::vector<unsigned int>> family_data = read_instance(INSTANCE_PATH);

    Assignment A(family_data, initial_solution);A.stats();

//    std::vector<unsigned int> initial_solution2 = read_solution("./solutions/flow_solution_71441_59883784.csv");
//    Assignment B(family_data, initial_solution2);
//    compare_solutions(A, B);
//    return 0;

//    try_local_search(A); return 0;

    Presets presets (&family_data);//, "./limited_presets.txt");

//    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
//        unsigned int k = A.get_ith_family(i)->get_k();
//        presets.assign_family(i, k, false);
//    }

//    tiny_test(); return 0;
    test_stuff(presets, initial_solution, 69153.77343); return 0;

    presets.compute_all_bounds();
    std::vector<uint_pair> true_occupancy_bounds = read_bounds("./true_occupancy_bounds.txt");
    for (unsigned int j = 0; j < NB_DAYS; j++){
        presets.prescribe_occupancy_lb(j, true_occupancy_bounds[j].first);
        presets.prescribe_occupancy_ub(j, true_occupancy_bounds[j].second);
    }
//    std::vector<uint_pair> assignations = read_assignations("./true_assignations.txt");
//    for (unsigned int m = 0; m < assignations.size(); m++){
//        presets.assign_family(assignations[m].first, assignations[m].second, false);
//    }
    presets.compute_all_bounds();
    RLLB rllb(presets);
    rllb.read_lambda(presets, "./lambda_2.txt");
//    presets.prescribe_occupancy_lb(57, 126);
//    rllb.compute_true_day_occupancy_bounds(presets, 68889); return 0;
//    rllb.suggest_best_branching_family(presets);return 0;
//    std::cout << rllb.get_lb() << std::endl; return 0;

//    rllb.carry_out_tests(presets);return 0;

//    compute_true_day_occupancy_bounds(presets, rllb);return 0;

//    compute_compulsory_assignations(presets, rllb);
//    compute_forbidden_assignations(presets, rllb);
//    presets.write_presets("./limited_presets.txt");
//    rllb.write_lambda("./lambda_4.txt");
//    return 0;

    // std::vector<bool> useful_families(NB_FAMILIES, false);
    // get_differences_between_solutions(initial_solution, read_solution("./solutions/flow_solution_71441_79041694.csv"), useful_families);
    // farm(A, presets, useful_families, 5000); return 0;
//    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
//        unsigned int k = A.get_ith_family(i)->get_k();
//        if ((i < 4500) &&  k <= 1)
////        if (i %10 != 9 && k <= 1)
//        //if (k <= 1)
//            presets.assign_family(i, k, false);
//    }
//    for (unsigned int i = 0; i < NB_DAYS; i++){
//        if (A.get_ith_day(i)->is_Friday_like())
//            presets.prescribe_occupancy_ub(i, MIN_NB_PEOPLE_PER_DAY);
//    }
//    for (unsigned int i = 0; i < NB_DAYS; i++){
//        unsigned int delta = 20;
//        if (i >= 25 && i <= 5*NB_DAYS/8) delta = 200;
//        if (i == 21) delta = 200;
//        //if (i == 48) delta = 20; //30, 33, 36 (30: 68775; 33: 68709; 36 very false: > 69200)
//        // 39: 68793; 42: 68807; 45: 68680; 46, 47, 48: ~68750
//        presets.prescribe_occupancy_lb(i, (unsigned int)(std::max(125, int(A.get_ith_day(i)->get_N() - delta))));
//        presets.prescribe_occupancy_ub(i, A.get_ith_day(i)->get_N() + delta);
//    }


//    try_anticipating_stuff(A, presets);
//    return 0;

//    Graph G = Graph(presets);
//    G.compute_max_flow_min_cost();
//    std::cout << "flow is " << G.get_current_flow() << " with current cost " << G.get_flow_cost() <<" and with true cost " << G.get_true_flow_cost() << std::endl;
//    G.check_flow();
//    G.show_schedule();

//    std::vector<unsigned int> reference_costs(NB_DAYS, 0);
//    for (unsigned int i = 0; i < NB_DAYS; i++)
//        reference_costs[i] = A.get_ith_day(i)->get_assignments_costs();
    unsigned int nb_nodes = 0;
    clock_t t0 = clock();
    clock_t time_graph = 0;
    clock_t time_bounds = 0;
//    clock_t time_anticipating = 0;
    float* best_cost = new float(BEST_SOLUTION);
    std::cout << "Best solution seen is " << oriented_branch_and_bound(presets, nb_nodes, time_graph, time_bounds, //time_anticipating,
                                                         rllb, best_cost) << std::endl;
//                                                         reference_costs, G.get_distribution()) << std::endl;
    std::cout << "Number of nodes is " << nb_nodes << std::endl;
    std::cout << "Time spent computing lagrangian bound is  " << time_graph << std::endl;
    //std::cout << "Time spent computing flows is  " << time_graph << std::endl;
    std::cout << "Time spent computing bounds is " << time_bounds << std::endl;
    //std::cout << "Time spent anticipating is     " << time_anticipating << std::endl;
    std::cout << clock() - t0 << std::endl;

    delete[] best_cost;
    return 0;
}


