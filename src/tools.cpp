#include "tools.h"
#include <algorithm>

std::vector<unsigned int> read_solution(const std::string &filename){
    std::ifstream targetFile (filename.c_str());
    if (!targetFile.is_open()) throw std::runtime_error("No targets file found");
    std::string input_line, entity;
    std::vector<unsigned int> res;
    // ignore header
    std::getline(targetFile, input_line);
    while (std::getline(targetFile, input_line)){
        std::stringstream line(input_line);
        std::getline(line, entity, ',');
        std::getline(line, entity, ',');
        res.push_back(std::stoi(entity) - 1);
    }
    return res;
}

std::vector<std::vector<unsigned int>> read_instance(const std::string &filename){
    std::ifstream targetFile (filename.c_str());
    if (!targetFile.is_open()) throw std::runtime_error("No targets file found");
    std::string input_line, entity;
    std::vector<std::vector<unsigned int>> res;
    // ignore header
    std::getline(targetFile, input_line);
    while (std::getline(targetFile, input_line)){
        std::vector<unsigned int> family_data;
        std::stringstream line(input_line);
        // ignore family_id
        std::getline(line, entity, ',');
        for (unsigned int i = 0; i < NB_CHOICES; i++){
            std::getline(line, entity, ',');
            family_data.push_back(std::stoi(entity) - 1); // we count days from 0 to 99, not from 1 to 100
        }
        std::getline(line, entity, ',');
        family_data.push_back(std::stoi(entity)); // number of people in the family
        res.push_back(family_data);
    }
    return res;
}

void write_solution_(std::vector<unsigned int> solution, const std::string &filename) {
    std::ofstream write(filename.c_str());
    write << "family_id,assigned_day" << std::endl;
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        write << i << "," << solution[i] + 1 << std::endl; // we restore here days from 1 to 100, never before
    }
}

unsigned int nb_chiffres(const unsigned int &i){
    if (i < 2)
        return 1;
    else
        return 1 + (unsigned int)(log(i)/log(10));
}

preset get_assignation_preset(const unsigned int &k) {
    preset result = std::vector<status>(K_MAX, FORBIDDEN);
    result[k] = COMPULSORY;
    return result;
}

preset get_counter_assignation_preset(const unsigned int &k) {
    preset result = std::vector<status>(K_MAX, ALLOWED);
    result[k] = FORBIDDEN;
    return result;
}

bool is_an_assignation(const preset &p) {
    unsigned int count = 0;
    for (unsigned int i = 0; i < K_MAX; i++)
        if (p[i] == FORBIDDEN)
            count++;
    return (count == K_MAX - 1);
}

bool myfunction (uint_pair i, uint_pair j) { return (i.second < j.second);}

void sort_by_second(std::vector<uint_pair> &myvector) {
    std::sort (myvector.begin(), myvector.end(), myfunction);
}

void print_nicely(const std::vector<unsigned int> &my_vector, const unsigned int &size_ub) {
    for (unsigned int i = 0; i < my_vector.size(); i++) {
        unsigned int n = nb_chiffres(my_vector[i]);
        if (n >= size_ub)
            throw std::logic_error("Why is this integer so big ?");
        std::cout << my_vector[i] << std::string(size_ub - n, ' ');
    }
    std::cout << std::endl;
}

