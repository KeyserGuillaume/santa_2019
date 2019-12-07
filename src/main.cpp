#include <iostream>
#include<string>
#include <fstream>
#include <sstream>
#include <vector>

#include "tools.h"
#include "Flow.h"
#include "LocalSearch.h"


int main() {
    std::vector<unsigned int> greedy_solution = read_solution("../../solutions/greedy_solution.csv");
    std::vector<std::vector<unsigned int>> family_data = read_instance(INSTANCE_PATH);

    Assignment A(family_data, greedy_solution);
    std::cout << A.get_cost() << std::endl;
    A.check_solution_is_ok();
    LocalSearch LS(&A);
    LS.run(10000000, 1000000);
    std::cout << A.get_cost() << std::endl;
    A.stats();

//    clock_t begin = clock();
//    LS.run_on_time_limit(begin + 10*60*CLOCKS_PER_SEC, 1000000);

    A.write_solution("../../solutions/local_search_solution_" + std::to_string(A.get_cost()) + "_.csv");


//    std::vector<unsigned int> schedule(NB_DAYS, 0);
//    for (unsigned int i = 0; i < NB_FAMILIES; i++){
//        schedule[greedy_solution[i]] += family_data[i][NB_CHOICES];
//    }
//
//    Graph G = Graph(schedule);
//    G.compute_max_flow_min_cost();
//    std::vector<unsigned int> flow_solution = G.get_solution();
//    write_solution(flow_solution, "../../solutions/flow_solution.csv");

    return 0;
}