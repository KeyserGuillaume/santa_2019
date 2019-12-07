#pragma once

#include<string>
#include <fstream>
#include <sstream>
#include <vector>

#include "constants.h"

std::vector<unsigned int> read_solution(const std::string &filename);
std::vector<std::vector<unsigned int>> read_instance(const std::string &filename);
void write_solution_(std::vector<unsigned int> solution, const std::string &filename);