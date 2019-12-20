#include <iostream>
#include<string>
#include <fstream>
#include <sstream>
#include <vector>

#include "tools.h"
#include "Flow.h"
#include "LocalSearch.h"


int main() {
    clock_t begin = clock();
    std::vector<unsigned int> initial_solution = read_solution("../../solutions/local_search_solution_72604_.csv");
    //initial_solution = read_solution("../../solutions/flow_solution.csv");
    std::vector<std::vector<unsigned int>> family_data = read_instance(INSTANCE_PATH);

    Assignment A(family_data, initial_solution);A.stats();//return 0;
//    std::cout << A.get_cost() << std::endl;
//    A.check_solution_is_ok();
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


    std::vector<preset> presets = get_empty_presets();
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        unsigned int k = A.get_ith_family(i)->get_k();
        if (i % 100 != 0 && k < K_MAX)
        //if (k <= 1)
            presets[i] = get_assignation_preset(k);
    }

    Graph G = Graph(family_data, presets);
//    G1.clear_flow();
//    for (unsigned int i = 0; i < NB_FAMILIES; i++){
//        if (i % 2 == 0)
//            G1.add_flow_for_assigning_family(i, A.get_ith_family(i)->get_k(), family_data[i][NB_CHOICES]);
//    }
    G.compute_max_flow_min_cost();
    std::cout << "flow is " << G.get_current_flow() << " with current cost " << G.get_flow_cost() <<" and with true cost " << G.get_true_flow_cost() << std::endl;
    G.check_flow();
    //G.show_schedule();

    unsigned int current_family = 0;
    unsigned int k = 0;  // the k we are assigning
    unsigned int value_if_assign, value_if_counter_assign;
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        bool found = false;
        for (;k < K_MAX; k++){
            for (; current_family < NB_FAMILIES; current_family++){
                if (!is_an_assignation(presets[current_family]) &&
                    presets[current_family][k] != FORBIDDEN){
                    found = true;
                    break;
                }
            }
            if (found) break;
            current_family = 0;
        }

        if (k == K_MAX)
            break;

        presets[current_family][k] = FORBIDDEN;
        if (is_an_assignation(presets[current_family]))
            presets[current_family][k + 1] = COMPULSORY;
        if (are_presets_feasible(presets, family_data)) {
            Graph G1 = Graph(family_data, presets);
            G1.compute_max_flow_min_cost();
            std::cout << "flow is " << G1.get_current_flow() << " with current cost " << G1.get_flow_cost()
                      << " and with true cost " << G1.get_true_flow_cost() << std::endl;
            if (!G1.is_flow_maximal())
                value_if_counter_assign = 100*UPPER_BOUND;
            else
                value_if_counter_assign = G1.get_true_flow_cost();
        }
        else
            value_if_counter_assign = 100*UPPER_BOUND;

        preset p = presets[current_family];

        presets[current_family] = get_assignation_preset(k);
        if (are_presets_feasible(presets, family_data)) {
            Graph G = Graph(family_data, presets);
            G.compute_max_flow_min_cost();
            std::cout << "flow is " << G.get_current_flow() << " with current cost " << G.get_flow_cost()
                      << " and with true cost " << G.get_true_flow_cost() << std::endl;
            value_if_assign = G.get_true_flow_cost();
            if (!G.is_flow_maximal())
                value_if_assign = 100*UPPER_BOUND;
            else
                value_if_assign = G.get_true_flow_cost();
        }
        else
            value_if_assign = 100*UPPER_BOUND;

        if (value_if_counter_assign < value_if_assign)
            presets[current_family] = p;
    }

    write_solution_(family_data, presets, "../../solutions/flow_solution.csv");
    Assignment A2(family_data, get_solution(family_data, presets));A2.get_cost();A2.stats();




//    Graph G = Graph();
//    std::cout << "flow is " << G.get_current_flow() << std::endl;
//    G.find_and_apply_augmenting_path();
//    std::vector<Arc *> path = G.get_shortest_path();
//    std::cout << "shortest path" << std::endl;
//    for (unsigned int i = 0; i < path.size(); i++){
//        std::cout << path[i]->get_id() << std::endl;
//    }
//    G.show_distances();
//    std::cout << "flow is " << G.get_current_flow() << std::endl;
//    G.find_and_apply_augmenting_path();
//    std::cout << "flow is " << G.get_current_flow() << std::endl;
//    G.find_and_apply_augmenting_path();
//    std::cout << "flow is " << G.get_current_flow() << std::endl;
//
//    G = Graph();
//    G.compute_max_flow_min_cost();
//    std::cout << "flow is " << G.get_current_flow() << " with cost " << G.get_flow_cost() << std::endl;

        // std::vector<unsigned int> flow_solution = G.get_solution();
    //write_solution(flow_solution, "../../solutions/flow_solution.csv");

    return 0;
}