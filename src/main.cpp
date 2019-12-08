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
    std::vector<unsigned int> greedy_solution = read_solution("../../solutions/local_search_solution_74565_.csv");
    std::vector<std::vector<unsigned int>> family_data = read_instance(INSTANCE_PATH);

    Assignment A(family_data, greedy_solution);
//    Family* f1 = A.get_ith_family(34);     // testing the move oriented cycles
//    Family* f2 = A.get_ith_family(35);
//    Day* d = f1->get_assigned_day();
//    f1->set_assigned_day(f2->get_assigned_day());
//    f2->set_assigned_day(d);
    std::cout << A.get_cost() << std::endl;
    A.check_solution_is_ok();
    LocalSearch LS(&A);
    //LS.run(10000000, 1000000);
    LS.run(10000, 1000);
    //LS.run_on_time_limit(begin + 60*60*CLOCKS_PER_SEC, 1000000);

    std::cout << A.get_cost() << std::endl;
    A.stats();
    LS.stats();

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