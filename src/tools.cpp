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


void write_solution_(const std::vector<std::vector<unsigned int>> &family_data, std::vector<preset> presets, const std::string &filename) {
    write_solution_(get_solution(family_data, presets), filename);
}

std::vector<unsigned int> get_solution(const std::vector<std::vector<unsigned int>> &family_data, std::vector<preset> presets) {
    std::vector<unsigned int> solution(0);
    unsigned int last_seen_allowed;
    unsigned int nb_allowed_seen;
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        nb_allowed_seen = 0;
        for (unsigned int k = 0; k < K_MAX; k++) {
            if (presets[i][k] == COMPULSORY)
                solution.push_back(family_data[i][k]);
            if (presets[i][k] == ALLOWED){
                nb_allowed_seen++;
                last_seen_allowed = family_data[i][k];
            }
        }
        if (nb_allowed_seen == 1 && solution.size() == i)
            solution.push_back(last_seen_allowed);
        else if (solution.size() == i)
            throw std::logic_error("Do these presets define a solutions ?");
    }
    return solution;
}


unsigned int nb_chiffres(unsigned int i){
    if (i < 2)
        return 1;
    else
        return 1 + (unsigned int)(log(i)/log(10));
}

preset get_assignation_preset(unsigned int k) {
    preset result = std::vector<status>(K_MAX, FORBIDDEN);
    result[k] = COMPULSORY;
    return result;
}

preset get_counter_assignation_preset(unsigned int k) {
    preset result = std::vector<status>(K_MAX, ALLOWED);
    result[k] = FORBIDDEN;
    return result;
}

std::vector<preset> get_empty_presets() {
    std::vector<preset> presets(0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        presets.push_back(std::vector<status>(K_MAX, ALLOWED));
    return presets;
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

bool are_presets_feasible(const std::vector<preset> &presets, const std::vector<std::vector<unsigned int>> &family_data) {
    std::vector<unsigned int> presets_schedule(NB_DAYS, 0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        for (unsigned int k = 0; k < K_MAX; k++) {
            if (presets[i][k] == COMPULSORY) {
                presets_schedule[family_data[i][k]] += family_data[i][NB_CHOICES];
                if (presets_schedule[family_data[i][k]] > MAX_NB_PEOPLE_PER_DAY)
                    return false;
            }
        }
    }
    return true;
}

