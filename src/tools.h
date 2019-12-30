#pragma once

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "constants.h"
#include <math.h>

template<typename T1, typename T2>
bool myfunction(std::pair<T1, T2> i, std::pair<T1, T2> j) {
    return (i.second < j.second);
}

template<typename T1, typename T2>
void sort_by_second(std::vector<std::pair<T1, T2>> &list) {
    std::sort (list.begin(), list.end(), myfunction<T1, T2>);
}

unsigned int nb_chiffres(const unsigned int &i);
void print_nicely(const std::vector<unsigned int> &my_vector, const unsigned int &size_ub = 4);
std::vector<bool> get_binary_rep(unsigned int x, const unsigned int& n);
std::vector<unsigned int> get_unique_vector(const std::vector<unsigned int> &v);

class PriorPath{
public:
    int cost;
    unsigned int family_index;
    std::vector<unsigned int> turns;
    unsigned int flow_quantity;
    PriorPath(int cost, unsigned int family_index, std::vector<unsigned int> turns, unsigned int flow_quantity = 0):
            cost(cost), family_index(family_index), turns(turns), flow_quantity(flow_quantity){}
};

struct PriorPath_compare{
    bool operator()(PriorPath i, PriorPath j);
};
//void sort_prior_paths(std::vector<PriorPath> &list);

std::vector<unsigned int> read_solution(const std::string &filename);
std::vector<std::vector<unsigned int>> read_instance(const std::string &filename);
void write_solution_(std::vector<unsigned int> solution, const std::string &filename);
preset get_assignation_preset(const unsigned int &k);
preset get_counter_assignation_preset(const unsigned int &k);
bool is_an_assignation(const preset& p);
