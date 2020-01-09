#include "tools.h"
#include "Flow.h"
#include "LocalSearch.h"
#include "Presets.h"
#include "day_sub_problem.h"

std::vector<uint_pair> do_assignments_by_default(Presets& presets){
    std::vector<uint_pair> assignments_by_default = presets.get_assignments_by_default();
    for (unsigned int m = 0; m < assignments_by_default.size(); m++)
        presets.assign_family(assignments_by_default[m].first, assignments_by_default[m].second, false);
    if (assignments_by_default.size() > 0)
        presets.compute_all_bounds(); presets.compute_feasibility();
    return assignments_by_default;
}

void undo_assignments_by_default(Presets& presets, const std::vector<uint_pair> &assignments_by_default){
    for (unsigned int m = 0; m < assignments_by_default.size(); m++) {
        presets.deassign_family(assignments_by_default[m].first, false);
        for (unsigned int k = 0; k < K_MAX; k++)
            if (k != assignments_by_default[m].second)
                presets.forbid_assignment(assignments_by_default[m].first, k, false);
    }
    if (assignments_by_default.size() > 0)
        presets.compute_all_bounds(); presets.compute_feasibility();
}

unsigned int brute_force(Presets &presets,
                         unsigned int& nb_nodes,
                         clock_t &time_graph,
                         clock_t &time_bounds,
                         clock_t &time_anticipating,
                         std::vector<unsigned int> reference_costs,
                         const FamilyDistribution &distribution){
    unsigned int cost = presets.get_presets_costs() + presets.get_day_cost_lb();
    if (cost > BEST_SOLUTION)
        return cost;
    if (presets.get_nb_assignments() == NB_FAMILIES){
        presets.write_solution("flow_solution");
        return cost;
    }

    // do the default assignments (families which only have one choice left)
    std::vector<uint_pair> assignments_by_default = do_assignments_by_default(presets);
    if (!presets.is_feasible()) {
        undo_assignments_by_default(presets, assignments_by_default);
        return BEST_SOLUTION + 1;
    }

    nb_nodes++;
    clock_t t0 = clock();
    Graph G = Graph(presets);
    G.set_prior_paths(distribution);
    G.compute_max_flow_min_cost();
    time_graph += clock() - t0;
    std::cout << "flow is " << G.get_current_flow() << " with current cost " << G.get_flow_cost()
              << " and with true cost " << G.get_true_flow_cost() << std::endl;
    cost = G.get_true_flow_cost();
    if (!G.is_flow_maximal() || cost > BEST_SOLUTION) {
        // indicate that there is no acceptable solution from this preset
        undo_assignments_by_default(presets, assignments_by_default);
        return BEST_SOLUTION + 1;
    }
    if (G.get_current_flow() == 0) {
        G.check_day_costs_are_ok();
        presets.write_solution("flow_solution");
        undo_assignments_by_default(presets, assignments_by_default);
        return cost;
    }

    if (G.get_most_dispersed_family().second == 1 && G.day_costs_are_ok()) {
        G.check_day_costs_are_ok();
        G.get_solution().write_solution("flow_solution");
        undo_assignments_by_default(presets, assignments_by_default);
        return cost;
    }

    unsigned int current_family = (unsigned int)(presets.get_largest_unassigned_families()[0]);

    FamilyDistribution new_distribution = G.get_distribution();
    preset old_preset = presets.get_preset(current_family);
    unsigned int min_cost = BEST_SOLUTION + 1;
    for (unsigned int k_assign = 0; k_assign < K_MAX; k_assign++){
        if (old_preset[k_assign] == FORBIDDEN)
            continue;

        // we branch on current_family : here we assign it to k_assign.
        t0 = clock();
        presets.assign_family(current_family, k_assign);
        time_bounds += clock() - t0;

        // only do the next parts if it looks feasible
        if (presets.is_feasible()) {
            // anticipation stuff
            t0 = clock();
            std::vector<std::pair<unsigned int, preset>> saves_of_presets(0);
            std::vector<std::vector<unsigned int>> assignations(K_MAX, std::vector<unsigned int>(0));
            std::vector<std::vector<unsigned int>> counter_assignations(K_MAX, std::vector<unsigned int>(0));
            for (unsigned int k = 0; k < K_MAX; k++){
                if (old_preset[k] == FORBIDDEN)
                    continue;
                unsigned int current_day = presets.get_family_data(current_family, k);
                anticipate_on_day(presets, current_day, reference_costs[current_day], assignations[k], counter_assignations[k]);
                for (unsigned int& family_index: assignations[k]){
                    saves_of_presets.push_back(std::pair<unsigned int, preset>(family_index, presets.get_preset(family_index)));
                    presets.assign_family(family_index, current_day, false, false);
                }
                for (unsigned int& family_index: counter_assignations[k]){
                    presets.forbid_assignment(family_index, current_day, false, false);
                }
                if (assignations[k].size() > 0 || counter_assignations[k].size() > 0)
                    presets.compute_all_bounds(); presets.compute_feasibility();
            }
            time_anticipating += clock() - t0;

            // we have assigned current_family, and hopefully we've anticipated some assignations
            // and counter_assignations. Now we call the function to go down in the exploration tree.
            cost = brute_force(presets, nb_nodes, time_graph, time_bounds, time_anticipating, reference_costs, new_distribution);
            if (cost < min_cost) {
                min_cost = cost;
            }

            // now we will remove the anticipated stuff, it is no longer valid when we assign to another day.
            t0 = clock();
            for (unsigned int m = 0; m < saves_of_presets.size(); m++)
                presets.set_preset(saves_of_presets[m].first, saves_of_presets[m].second, false);
            for (unsigned int k = 0; k < K_MAX; k++)
                for (unsigned int s = 0; s < counter_assignations[k].size(); s++)
                    presets.enable_assignment(counter_assignations[k][s], presets.get_family_data(current_family, k), false, false);
            presets.compute_all_bounds(); presets.compute_feasibility();
            time_anticipating += clock() - t0;
        }

    }

    // now we put the preset for current_family back to what it was before
    t0 = clock();
    presets.set_preset(current_family, old_preset);
    time_bounds += clock() - t0;

    // finally we do the same for the assignments by default done at the beginning
    undo_assignments_by_default(presets, assignments_by_default);
    return min_cost;
}

unsigned int brute_force(Presets &presets,
                         unsigned int& nb_nodes,
                         clock_t &time_lagrangian_lb,
                         clock_t &time_bounds,
                         const std::vector<float> &lambda,
                         float &best_cost){

    // do the default assignments (families which only have one choice left)
    std::vector<uint_pair> assignments_by_default = do_assignments_by_default(presets);
    float cost = presets.get_presets_costs() + presets.get_day_cost_lb();
    if (!presets.is_feasible() || cost >= best_cost) {
        undo_assignments_by_default(presets, assignments_by_default);
        return best_cost + 1;
    }
    if (presets.get_nb_assignments() == NB_FAMILIES){
        presets.write_solution("lagrangian_bound_solution");
        undo_assignments_by_default(presets, assignments_by_default);
        return cost;
    }

    nb_nodes++;
    clock_t t0 = clock();
    std::vector<float> new_lambda(lambda);
    bool is_primal_feasible;
    cost = get_lagrangian_lb(presets, new_lambda, is_primal_feasible, true);
    // there is no proof that this is the local optimal solution but at least this solution is better than best.
    if (is_primal_feasible) {
        best_cost = cost;
        return cost;
    }
    if (cost >= best_cost) {
        // indicate that there is no acceptable solution from this preset
        undo_assignments_by_default(presets, assignments_by_default);
        return best_cost + 1;
    }

    time_lagrangian_lb += clock() - t0;

    unsigned int current_family = (unsigned int)(presets.get_largest_unassigned_families()[0]);
    new_lambda[current_family] = 0;

    preset old_preset = presets.get_preset(current_family);
    unsigned int min_cost = best_cost + 1;
    for (unsigned int k_assign = 0; k_assign < K_MAX; k_assign++){
        if (old_preset[k_assign] == FORBIDDEN)
            continue;

        // we branch on current_family : here we assign it to k_assign.
        t0 = clock();
        presets.assign_family(current_family, k_assign);
        time_bounds += clock() - t0;

        // only do the next parts if it looks feasible
        if (presets.is_feasible()) {
            cost = brute_force(presets, nb_nodes, time_lagrangian_lb, time_bounds, new_lambda, best_cost);
            if (cost < min_cost) {
                min_cost = cost;
            }
        }
    }

    // now we put the preset for current_family back to what it was before
    t0 = clock();
    presets.set_preset(current_family, old_preset);
    time_bounds += clock() - t0;

    // finally we do the same for the assignments by default done at the beginning
    undo_assignments_by_default(presets, assignments_by_default);
    return min_cost;
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
        unsigned int i_0 = rand() % NB_DAYS;
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


int main() {
    clock_t begin = clock();
    std::vector<unsigned int> initial_solution = read_solution("../../solutions/flow_solution_71441_59883784.csv");
    //std::vector<unsigned int> initial_solution = read_solution("../../solutions/local_search_solution_71480_.csv");
    //std::vector<unsigned int> initial_solution = read_solution("../../solutions/alocal_search_solution_74253_.csv");
    std::vector<std::vector<unsigned int>> family_data = read_instance(INSTANCE_PATH);

    Assignment A(family_data, initial_solution);A.stats();
//
//    std::vector<unsigned int> initial_solution2 = read_solution("../../solutions/flow_solution_71441_59883784.csv");
//    Assignment B(family_data, initial_solution2);
//    for (unsigned int i = 0; i < NB_DAYS; i++)
//        if (A.get_ith_day(i)->get_N() != B.get_ith_day(i)->get_N())
//            std::cout << "Day " << i << ": " << A.get_ith_day(i)->get_N() << " -> " << B.get_ith_day(i)->get_N()<<std::endl;
//    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
//        unsigned int k_A = A.get_ith_family(i)->get_k();
//        unsigned int k_B = B.get_ith_family(i)->get_k();
//        unsigned int day_A = A.get_ith_family(i)->get_assigned_day()->get_id();
//        unsigned int day_B = B.get_ith_family(i)->get_assigned_day()->get_id();
//        if (k_B != k_A) {
//            std::cout << "Family " << i << " (" << A.get_ith_family(i)->get_nb_people() << "): ";
//            std::cout << day_A << " -> " << day_B << " (" << k_A << " -> " << k_B << ")" << std::endl;
//        }
//    }
//    return 0;

//    LocalSearch LS(&A);
//    LS.run(10000000, 1000000);
//    //LS.run(1000, 100);
//    //LS.run_on_time_limit(begin + 60*60*CLOCKS_PER_SEC, 1000000);
//
//    std::cout << A.get_cost() << std::endl;
//    A.stats();
//    LS.stats();
//
//    A.write_solution("../../solutions/alocal_search_solution_" + std::to_string(A.get_cost()) + "_.csv");
//    return 0;

    Presets presets (family_data);
    presets.compute_all_bounds();
    std::vector<bool> useful_families(NB_FAMILIES, false);
    //get_differences_between_solutions(initial_solution, read_solution("../../solutions/flow_solution_71441_79041694.csv"), useful_families);
    //get_differences_between_solutions(initial_solution, read_solution("../../solutions/flow_solution_71490_47394673.csv"), useful_families);
    //get_differences_between_solutions(initial_solution, read_solution("../../solutions/flow_solution_71490_71088485.csv"), useful_families);
    //get_differences_between_solutions(initial_solution, read_solution("../../solutions/flow_solution_71490_81100696.csv"), useful_families);
    //farm(A, presets, useful_families, 5000); return 0;
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        unsigned int k = A.get_ith_family(i)->get_k();
        //if (i % 5 == 0)
        if (i < 4200)
        //if (i % 5 != 0)
        //if (i > 100 && k <= 2)
        //if (i % 50 != 49 && k <= 2)
        //if (k <= 1)
            presets.assign_family(i, k, false);
    }
    for (unsigned int i = 0; i < NB_DAYS; i++){
        if (A.get_ith_day(i)->is_Friday_like())
            presets.prescribe_occupancy_ub(i, MIN_NB_PEOPLE_PER_DAY);
    }
    for (unsigned int i = 0; i < NB_DAYS; i++){
        presets.prescribe_occupancy_lb(i, A.get_ith_day(i)->get_N() - 15);
        presets.prescribe_occupancy_ub(i, A.get_ith_day(i)->get_N() + 15);
    }
    presets.compute_all_bounds(); presets.compute_feasibility();
    //try_lagrangian_thing(presets); return 0;

//    for (int i = NB_DAYS - 1; i >= 1; i--) {
//        Day *day = A.get_ith_day(i);
//        unsigned int cost = day->get_assignments_costs();
//        std::cout << "We try assigning to day " << i << " of cost " << cost << std::endl;
//        std::vector<unsigned int> assignations(0), counter_assignations(0);
//        anticipate_on_day(presets, i, cost, assignations, counter_assignations);
//        for (unsigned int& family_index: assignations){
//            presets.assign_family(family_index, i, false, false);
//        }
//        for (unsigned int& family_index: counter_assignations){
//            presets.forbid_assignment(family_index, i, false, false);
//        }
//        presets.compute_all_bounds(); presets.compute_feasibility();
//    }
//    return 0;

    Graph G = Graph(presets);
    G.compute_max_flow_min_cost();
    std::cout << "flow is " << G.get_current_flow() << " with current cost " << G.get_flow_cost() <<" and with true cost " << G.get_true_flow_cost() << std::endl;
    G.check_flow();
    G.show_schedule();

    bool is_primal_feasible;
    std::vector<float> lambda(NB_FAMILIES, 0);
    std::cout << "The lagrangian bound ended with " << get_lagrangian_lb(presets, lambda, is_primal_feasible, true) << std::endl;
    //if (is_primal_feasible) return 0;

    std::vector<unsigned int> reference_costs(NB_DAYS, 0);
    for (unsigned int i = 0; i < NB_DAYS; i++)
        reference_costs[i] = A.get_ith_day(i)->get_assignments_costs();
    unsigned int nb_nodes = 0;
    clock_t t0 = clock();
    clock_t time_graph = 0;
    clock_t time_bounds = 0;
    clock_t time_anticipating = 0;
    float best_cost = BEST_SOLUTION;
    std::cout << "Best solution seen is " << brute_force(presets, nb_nodes, time_graph, time_bounds, //time_anticipating,
                                                         lambda, best_cost) << std::endl;
//                                                         reference_costs, G.get_distribution()) << std::endl;
    std::cout << "Number of nodes is " << nb_nodes << std::endl;
    std::cout << "Time spent computing lagrangian bound is  " << time_graph << std::endl;
    //std::cout << "Time spent computing flows is  " << time_graph << std::endl;
    std::cout << "Time spent computing bounds is " << time_bounds << std::endl;
    //std::cout << "Time spent anticipating is     " << time_anticipating << std::endl;
    std::cout << clock() - t0 << std::endl;

    return 0;
}

