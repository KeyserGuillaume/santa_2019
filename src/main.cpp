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
    //std::vector<unsigned int> initial_solution = read_solution("../../solutions/local_search_solution_74565_.csv");
    std::vector<std::vector<unsigned int>> family_data = read_instance(INSTANCE_PATH);

//    Assignment A(family_data, initial_solution);A.stats();
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


    Graph G = Graph(family_data);
    G.compute_max_flow_min_cost();
    std::cout << "flow is " << G.get_current_flow() << " with true cost " << G.get_true_flow_cost() << std::endl;

//    Graph G = Graph();
//    std::cout << "flow is " << G.get_current_flow() << std::endl;
//    std::vector<Arc *> path = G.get_shortest_path();
//    std::cout << "shortest path" << std::endl;
//    for (unsigned int i = 0; i < path.size(); i++){
//        std::cout << path[i]->get_id() << std::endl;
//    }
//    G.find_and_apply_augmenting_path();
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