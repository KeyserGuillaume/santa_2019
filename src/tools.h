#pragma once

#include<string>
#include <fstream>
#include <sstream>
#include <vector>

#include "constants.h"
#include <math.h>

typedef std::pair<unsigned int, unsigned int> uint_pair;

std::vector<unsigned int> read_solution(const std::string &filename);
std::vector<std::vector<unsigned int>> read_instance(const std::string &filename);
void write_solution_(std::vector<unsigned int> solution, const std::string &filename);
void write_solution_(const std::vector<std::vector<unsigned int>> &family_data, std::vector<preset> presets, const std::string &filename);
unsigned int nb_chiffres(unsigned int i);
std::vector<preset> get_empty_presets();
preset get_assignation_preset(unsigned int k);
preset get_counter_assignation_preset(unsigned int k);
bool is_an_assignation(const preset& p);
bool myfunction (uint_pair i, uint_pair j);
void sort_by_second(std::vector<uint_pair> &list);