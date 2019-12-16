#pragma once

#include <iostream>
#include <vector>
#include <time.h>

#include "constants.h"
#include "Assignment.h"
#include "tools.h"

class LocalSearch {
private:
    unsigned int k;
    Assignment* A;
    std::vector<unsigned int> count_order_change = std::vector<unsigned int>(13, 0);
    std::vector<unsigned int> count_size_moved_by_twin_paths = std::vector<unsigned int>(9, 0);
    unsigned int nb_successful_augmenting_path_0 = 0;
    unsigned int nb_successful_augmenting_path_1 = 0;
    unsigned int nb_successful_augmenting_path_2 = 0;
    unsigned int nb_successful_twin_paths = 0;
    unsigned int nb_successful_oriented_cycle = 0;

    void augmenting_path_move();
    void twin_paths_move();
    void oriented_cycle_move();
    void remove_Friday_like();
    void jump();

public:
    LocalSearch(Assignment* A): A(A){srand((unsigned int)time(0)); k = 0;}

    void run_on_time_limit(const clock_t time_limit, const int &period_display);
    void run(const long &nb_iteration, const int &period_display);
    void display() const{
        std::cout << k << " " << A->get_cost() << std::endl;
    }
    void stats() const;
};
