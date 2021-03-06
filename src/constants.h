#pragma once
#include <string>
#include <vector>
#include <unordered_set>
#include <set>

enum status {ALLOWED, FORBIDDEN, COMPULSORY};

typedef std::vector<status> preset;
typedef std::pair<unsigned int, unsigned int> uint_pair;
typedef std::pair<float, float> float_pair;
typedef std::pair<double, double> double_pair;
typedef std::vector<std::pair<unsigned int, std::vector<unsigned int>>> FamilyDistribution;

const unsigned int NB_FAMILIES = 5000;
const unsigned int NB_DAYS = 100;
const unsigned int NB_CHOICES = 10;
const unsigned int MIN_NB_PEOPLE_PER_DAY = 125;
const unsigned int MAX_NB_PEOPLE_PER_DAY = 300;
const unsigned int UPPER_BOUND = 10000; // was 80000; I think I put 10000 here (not an upper bound per se) in order to change behaviour of B&B
const unsigned int LARGE_UPPER_BOUND = 100000; // we know that the instance we work on has optimal value just under 70000 so 100000 is a know UB
const std::string INSTANCE_PATH = "./instance/family_data.csv";

const unsigned int CONSTANT_COST [] = {0, 50, 50, 100, 200, 200, 300, 300, 400, 500, 500};
const unsigned int MARGINAL_COST [] = {0,  0,  9,   9,   9,  18,  18,  36,  36, 235, 434};

const unsigned int K_MAX = 5;
//const unsigned int BEST_SOLUTION = 68889;
const unsigned int BEST_SOLUTION = 1000000;