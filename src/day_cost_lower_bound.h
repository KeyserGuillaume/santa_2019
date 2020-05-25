#pragma once

#include "constants.h"
#include "tools.h"

double get_day_cost(const unsigned int& i, const unsigned int& ip1);
double get_day_cost(const std::vector<unsigned int> &occupancies);
double get_binary_day_cost_lb(const unsigned int &a_i, const unsigned int &b_i, const unsigned int &a_ip1, const unsigned int &b_ip1);
double get_ternary_day_cost_lb(const unsigned int &a_im1, const unsigned int &b_im1, const unsigned int &a_i,
                             const unsigned int &b_i, const unsigned int &a_ip1, const unsigned int &b_ip1);
double get_day_cost_DP_lb(const std::vector<std::vector<unsigned int>> &possible_quantities);
double get_day_cost_DP_ultimate_lb(const std::vector<std::vector<unsigned int>> &possible_quantities, const unsigned int& s);