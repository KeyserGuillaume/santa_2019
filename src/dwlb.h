#pragma once

#include "constants.h"
#include "tools.h"
#include "Presets.h"

#include <ilcplex/ilocplex.h>

class Column{
    std::vector<std::vector<unsigned int>> selected_choices;
    double cost;
public:
    Column(const std::vector<std::vector<unsigned int>> &selected_choices, const double& cost): selected_choices(selected_choices), cost(cost){}
    unsigned int get_family_use_nb(const unsigned int& i) const { return selected_choices[i].size(); }
    double get_cost() const { return cost; }
};

void test_stuff(Presets& presets, const std::vector<unsigned int> &initial_solution, const double& cost);