#pragma once
#include <string>

const unsigned int NB_FAMILIES = 5000;
const unsigned int NB_DAYS = 100;
const unsigned int NB_CHOICES = 10;
const unsigned int MIN_NB_PEOPLE_PER_DAY = 125;
const unsigned int MAX_NB_PEOPLE_PER_DAY = 300;
const unsigned int UPPER_BOUND = 10000; // was 80000
const std::string INSTANCE_PATH = "../../instance/family_data.csv";

const unsigned int CONSTANT_COST [] = {0, 50, 50, 100, 200, 200, 300, 300, 400, 500, 500};
const unsigned int MARGINAL_COST [] = {0,  0,  9,   9,   9,  18,  18,  36,  36, 235, 434};