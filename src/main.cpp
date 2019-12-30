#include <iostream>
#include<string>
#include <fstream>
#include <sstream>
#include <vector>

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
    if (presets.get_nb_assignments() == NB_FAMILIES){
        unsigned int cost = presets.get_presets_costs() + presets.get_day_cost_lb();
        if (cost <= BEST_SOLUTION)
            presets.write_solution();
        return cost;
    }

    // do the default assignements (families which only have one choice left)
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
    unsigned int cost = G.get_true_flow_cost();
    if (!G.is_flow_maximal() || cost > BEST_SOLUTION) {
        // indicate that there is no acceptable solution from this preset
        undo_assignments_by_default(presets, assignments_by_default);
        return BEST_SOLUTION + 1;
    }
    if (G.get_current_flow() == 0) {
        G.check_day_costs_are_ok();
        presets.write_solution();
        undo_assignments_by_default(presets, assignments_by_default);
        return cost;
    }

    if (G.get_most_dispersed_family().second == 1 && G.day_costs_are_ok()) {
        G.get_solution().write_solution();
        undo_assignments_by_default(presets, assignments_by_default);
        return cost;
    }
    //unsigned int current_family = G.get_overload_family();
    unsigned int current_family = (unsigned int)(presets.get_largest_unassigned_families()[0]);
    //unsigned int current_family = presets.get_family_max_size_min_choices();
    //unsigned int current_family = G.get_largest_least_dispersed_family();
    //unsigned int current_family = presets.get_largest_unassigned_strategic_family();
//    if (current_family == - 1) {
//        G.get_solution().write_solution("../../solutions/flow_solution_" + std::to_string(cost) + ".csv");
//        return cost;
//    }
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

//unsigned int branching_on_bounds(Presets &presets,
//                         unsigned int& nb_nodes,
//                         clock_t &time_graph,
//                         clock_t &time_bounds){
//    if (presets.get_nb_assignments() == NB_FAMILIES){
//        unsigned int cost = presets.get_presets_costs() + presets.get_day_cost_lb();
//        if (cost <= BEST_SOLUTION)
//            presets.write_solution("../../solutions/flow_solution_" + std::to_string(cost) + ".csv");
//        return cost;
//    }
//    // this needs to be changed to account for forbidden assignations
//    nb_nodes++;
//    clock_t t0 = clock();
//    Graph G = Graph(presets);
//    G.compute_max_flow_min_cost();
//    time_graph += clock() - t0;
//    std::cout << "flow is " << G.get_current_flow() << " with current cost " << G.get_flow_cost()
//              << " and with true cost " << G.get_true_flow_cost() << std::endl;
//    unsigned int cost = G.get_true_flow_cost();
//    if (!G.is_flow_maximal() || cost > BEST_SOLUTION)
//        return BEST_SOLUTION + 1; // indicate that there is no acceptable solution from this preset
//    if (G.get_current_flow() == 0) {
//        G.check_day_costs_are_ok();
//        presets.write_solution("../../solutions/flow_solution_" + std::to_string(cost) + ".csv");
//        return cost;
//    }
//
//    if (G.get_most_dispersed_family().second == 1 && G.day_costs_are_ok()) {
//        G.get_solution().write_solution("../../solutions/flow_solution_" + std::to_string(cost) + ".csv");
//        return cost;
//    }
//
//    uint_pair branching_pair = presets.get_bounds_to_branch();
//    //branching_pair.second = G.get_day_occupancy()[branching_pair.first];
//    if (branching_pair.second - presets.get_occupancy_lb(branching_pair.first) >= 10) {
//        unsigned int min_cost = BEST_SOLUTION + 1;
//
//        unsigned int previous_assignation_cost = G.get_true_flow_cost() - presets.get_day_cost_lb();
//        presets.prescribe_occupancy_ub(branching_pair.first, branching_pair.second);
//        t0 = clock();
//        presets.compute_all_bounds();
//        presets.compute_feasibility();
//        time_bounds += clock() - t0;
//        if (presets.is_feasible() && presets.get_day_cost_lb() + previous_assignation_cost <= BEST_SOLUTION) {
//            cost = branching_on_bounds(presets, nb_nodes, time_graph, time_bounds);
//            if (cost < min_cost) {
//                min_cost = cost;
//            }
//        }
//        presets.pop_last_occupancy_ub_prescription();
//        presets.prescribe_occupancy_lb(branching_pair.first, branching_pair.second + 1);
//        t0 = clock();
//        presets.compute_all_bounds();
//        presets.compute_feasibility();
//        time_bounds += clock() - t0;
//        if (presets.is_feasible() && presets.get_day_cost_lb() + previous_assignation_cost <= BEST_SOLUTION) {
//            cost = branching_on_bounds(presets, nb_nodes, time_graph, time_bounds);
//            if (cost < min_cost) {
//                min_cost = cost;
//            }
//        }
//        presets.pop_last_occupancy_lb_prescription();
//
//        return min_cost;
//    } else {
//        return brute_force(presets, nb_nodes, time_graph, time_bounds);
//    }
//}

unsigned int get_day_responsible_for_day_cost_inaccuracy(const Presets &presets, const Graph &G, const unsigned int &i) {
    if (i == NB_DAYS) return i;
    std::vector<unsigned int> V = G.get_day_occupancy();
    unsigned int occ_i = V[i];
    unsigned int occ_ip1 = V[i + 1];
    unsigned int occ_lb_i = presets.get_occupancy_lb(i);
    unsigned int occ_ub_i = presets.get_occupancy_ub(i);
    unsigned int occ_lb_ip1 = presets.get_occupancy_lb(i + 1);
    unsigned int occ_ub_ip1 = presets.get_occupancy_ub(i + 1);
    unsigned int margin_i = occ_i <= occ_ip1 ? occ_ub_i - occ_i : occ_i - occ_lb_i;
    unsigned int margin_ip1 = occ_i <= occ_ip1 ? occ_ip1 - occ_lb_ip1 : occ_ub_ip1 - occ_ip1;
    return margin_i >= margin_ip1 ? i : i + 1;
}

//unsigned int mix_branching(Presets &presets,
//                         unsigned int& nb_nodes,
//                         clock_t &time_graph,
//                         clock_t &time_bounds){
//    if (presets.get_nb_assignments() == NB_FAMILIES){
//        unsigned int cost = presets.get_presets_costs() + presets.get_day_cost_lb();
//        if (cost <= BEST_SOLUTION) {
//            Graph G = Graph(presets);
//            G.check_day_costs_are_ok();
//            presets.write_solution("../../solutions/flow_solution_" + std::to_string(cost) + ".csv");
//        }
//        return cost;
//    }
//    // this needs to be changed to account for forbidden assignations
//
//    nb_nodes++;
//    clock_t t0 = clock();
//    Graph G = Graph(presets);
//    G.compute_max_flow_min_cost();
//    time_graph += clock() - t0;
//    std::cout << "flow is " << G.get_current_flow() << " with current cost " << G.get_flow_cost()
//              << " and with true cost " << G.get_true_flow_cost() << std::endl;
//    unsigned int cost = G.get_true_flow_cost();
//    if (!G.is_flow_maximal() || cost > BEST_SOLUTION)
//        return BEST_SOLUTION + 1; // indicate that there is no acceptable solution from this preset
//    if (G.get_current_flow() == 0) {
//        G.check_day_costs_are_ok();
//        presets.write_solution("../../solutions/flow_solution_" + std::to_string(cost) + ".csv");
//        return cost;
//    }
//
//    if (G.get_most_dispersed_family().second == 1 && G.day_costs_are_ok()) {
//        G.get_solution().write_solution("../../solutions/flow_solution_" + std::to_string(cost) + ".csv");
//        return cost;
//    }
//
//    int day_to_branch = -1;
//    std::vector<float> real_day_costs = G.get_real_day_costs();
//    for (unsigned int i = 0; i < NB_DAYS && day_to_branch == -1; i++){
//        if (real_day_costs[i] > presets.get_day_cost_lb(i) + 150) {
//            day_to_branch = get_day_responsible_for_day_cost_inaccuracy(presets, G, i);
//        }
//    }
//
//
//    unsigned int min_cost = BEST_SOLUTION + 1;
//    if (day_to_branch >= 0){
//        unsigned int middle = floor(0.5*presets.get_occupancy_ub(day_to_branch) + 0.5*presets.get_occupancy_lb(day_to_branch));
//        unsigned int previous_assignation_cost = G.get_true_flow_cost() - presets.get_day_cost_lb();
//        presets.prescribe_occupancy_ub(day_to_branch, middle);
//        t0 = clock();
//        presets.compute_all_bounds();
//        presets.compute_feasibility();
//        time_bounds += clock() - t0;
//        if (presets.is_feasible() && presets.get_day_cost_lb() + previous_assignation_cost <= BEST_SOLUTION) {
//            cost = mix_branching(presets, nb_nodes, time_graph, time_bounds);
//            if (cost < min_cost) {
//                min_cost = cost;
//            }
//        }
//        presets.pop_last_occupancy_ub_prescription();
//        presets.prescribe_occupancy_lb(day_to_branch, middle + 1);
//        t0 = clock();
//        presets.compute_all_bounds();
//        presets.compute_feasibility();
//        time_bounds += clock() - t0;
//        if (presets.is_feasible() && presets.get_day_cost_lb() + previous_assignation_cost <= BEST_SOLUTION) {
//            cost = mix_branching(presets, nb_nodes, time_graph, time_bounds);
//            if (cost < min_cost) {
//                min_cost = cost;
//            }
//        }
//        presets.pop_last_occupancy_lb_prescription();
//
//    } else {
//
//        unsigned int current_family = (unsigned int) (presets.get_largest_unassigned_families()[0]);
//        for (unsigned int k = 0; k < K_MAX; k++) {
//            t0 = clock();
//            presets.assign_family(current_family, k);
//            time_bounds += clock() - t0;
//            if (presets.is_feasible()) {
//                cost = mix_branching(presets, nb_nodes, time_graph, time_bounds);
//                if (cost < min_cost) {
//                    min_cost = cost;
//                }
//            }
//        }
//
//        t0 = clock();
//        presets.deassign_family(current_family);
//        time_bounds += clock() - t0;
//    }
//
//    return min_cost;
//}


int main() {
    clock_t begin = clock();
    std::vector<unsigned int> initial_solution = read_solution("../../solutions/local_search_solution_72604_.csv");
    std::vector<std::vector<unsigned int>> family_data = read_instance(INSTANCE_PATH);

    Assignment A(family_data, initial_solution);A.stats();

//    std::vector<unsigned int> initial_solution2 = read_solution("../../solutions/flow_solution_72504_36891366.csv");
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
//    A.write_solution("../../solutions/local_search_solution_" + std::to_string(A.get_cost()) + "_.csv");


//    std::vector<unsigned int> schedule(NB_DAYS, 0);
//    for (unsigned int i = 0; i < NB_FAMILIES; i++){
//        schedule[greedy_solution[i]] += family_data[i][NB_CHOICES];
//    }
    Presets presets (family_data);
    presets.compute_all_bounds();
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        unsigned int k = A.get_ith_family(i)->get_k();
        //if (i % 5 == 0)
        //if (i % 2 != 0)
        //if (i % 5 != 0)
        if (i > 100 && k <= 2)
        //if (i % 50 != 9 && k <= 2)
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


//    presets.assign_family(2993, 92, true, false);
//    find_all_equivalent_solutions(presets, 92, 1300);
//    return 0;

    for (int i = NB_DAYS - 1; i >= 0; i--) {
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
    //return 0;

    Graph G = Graph(presets);
    G.compute_max_flow_min_cost();
    std::cout << "flow is " << G.get_current_flow() << " with current cost " << G.get_flow_cost() <<" and with true cost " << G.get_true_flow_cost() << std::endl;
    G.check_flow();
    G.show_schedule();
    //std::cout << presets.get_occupancy_lb(branching_pair.first) << " " << presets.get_occupancy_ub(branching_pair.first) << " " << branching_pair.second<<std::endl;


//    Graph G1 = Graph(presets);
//    G1.set_prior_paths(G.get_distribution());
//    G1.compute_max_flow_min_cost();
//    std::cout << "flow is " << G1.get_current_flow() << " with current cost " << G1.get_flow_cost() <<" and with true cost " << G1.get_true_flow_cost() << std::endl;
//    G1.check_flow();
    std::vector<unsigned int> reference_costs(NB_DAYS, 0);
    for (unsigned int i = 0; i < NB_DAYS; i++)
        reference_costs[i] = A.get_ith_day(i)->get_assignments_costs();
    unsigned int nb_nodes = 0;
    clock_t t0 = clock();
    clock_t time_graph = 0;
    clock_t time_bounds = 0;
    clock_t time_anticipating = 0;
    std::cout << "Best solution seen is " << brute_force(presets, nb_nodes, time_graph, time_bounds, time_anticipating,
                                                         reference_costs, G.get_distribution()) << std::endl;
    std::cout << "Number of nodes is " << nb_nodes << std::endl;
    std::cout << "Time spent computing flows is  " << time_graph << std::endl;
    std::cout << "Time spent computing bounds is " << time_bounds << std::endl;
    std::cout << "Time spent anticipating is     " << time_anticipating << std::endl;
    std::cout << clock() - t0 << std::endl;

    return 0;
}


//        Graph G = Graph(presets);
//        for (unsigned int j = 0; j < NB_FAMILIES; j++) {
//            G.add_flow_for_assigning_family(j, A.get_ith_family(j)->get_k());
//            G.check_flow(false);
//        }
//        G.check_flow();
//        cost = G.get_true_flow_cost() - presets.get_day_cost_lb();
//        if (cost != 67433)
//            throw std::logic_error("lhbhbj");
//        Graph G1 = Graph(presets);
//        G1.compute_max_flow_min_cost();
//        unsigned int cost_G1 = G1.get_true_flow_cost();
//        unsigned int cost_G = G.get_true_flow_cost();
//        if (cost_G < cost_G1)
//            throw std::logic_error("rfdjk");
//        for (unsigned int j = 0; j < NB_DAYS; j++){
//            unsigned int reference_day_cost = A.get_ith_day(j)->get_cost();
//            unsigned int day_cost_lb = presets.get_day_cost_lb(j);
//            if (day_cost_lb > reference_day_cost)
//                throw std::logic_error("pmqnx;");
//        }
//        for (unsigned int j = 0; j < NB_FAMILIES; j++) {
//            Family *f = A.get_ith_family(j);
//            unsigned int real_k = f->get_k();
//            if (presets[j][real_k] == FORBIDDEN)
//                throw std::logic_error("flrjfew");
//            for (unsigned int k = 0; k < K_MAX; k++)
//                if (presets[j][k] == COMPULSORY && k != real_k)
//                    throw std::logic_error("flrjfew");
//        }


//unsigned int tmp_n = presets.get_nb_assignments();
//if (presets.get_nb_assignments() != tmp_n)
//throw std::logic_error("iuhgrd#{]");