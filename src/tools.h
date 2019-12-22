#pragma once

#include<string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

#include "constants.h"
#include <math.h>

typedef std::pair<unsigned int, unsigned int> uint_pair;

std::vector<unsigned int> read_solution(const std::string &filename);
std::vector<std::vector<unsigned int>> read_instance(const std::string &filename);
void write_solution_(std::vector<unsigned int> solution, const std::string &filename);
unsigned int nb_chiffres(unsigned int i);
preset get_assignation_preset(unsigned int k);
preset get_counter_assignation_preset(unsigned int k);
bool is_an_assignation(const preset& p);
bool myfunction (uint_pair i, uint_pair j);
void sort_by_second(std::vector<uint_pair> &list);
void print_nicely(const std::vector<unsigned int> &my_vector, const unsigned int &size_ub = 4);