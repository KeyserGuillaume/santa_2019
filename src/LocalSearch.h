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

    void augmenting_path_move();
    void jump();

public:
    LocalSearch(Assignment* A): A(A){srand((unsigned int)time(0)); k = 0;}

    void run_on_time_limit(const clock_t time_limit, const int &period_display);
    void run(const long &nb_iteration, const int &period_display);
    void display() const{
        std::cout << k << " " << A->get_cost() << std::endl;
    }
};
