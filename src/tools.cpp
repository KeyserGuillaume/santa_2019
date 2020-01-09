#include "tools.h"

std::vector<bool>
cleverly_get_possible_quantities(const unsigned int &max_allowed, const std::vector<unsigned int> &max_per_size);

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

void print_nicely(const std::vector<unsigned int> &my_vector, const unsigned int &size_ub) {
    for (unsigned int i = 0; i < my_vector.size(); i++) {
        unsigned int n = nb_chiffres(my_vector[i]);
        if (n >= size_ub)
            throw std::logic_error("Why is this integer so big ?");
        std::cout << my_vector[i] << std::string(size_ub - n, ' ');
    }
    std::cout << std::endl;
}

bool PriorPath_compare::operator()(PriorPath i, PriorPath j) {
    return i.cost > j.cost;
}

std::vector<bool> get_binary_rep(unsigned int x, const unsigned int &n) {
    std::vector<bool> res(n, false);
    for (int i = n - 1; i >= 0; i--){
        if (x >= pow(2, i)){
            x -= pow(2, i);
            res[n - i - 1] = true;
        }
    }
    return res;
}

std::vector<unsigned int> get_unique_vector(const std::vector<unsigned int> &v) {
    std::unordered_set<unsigned int> my_set(0);
    for (unsigned int i = 0; i < v.size(); i++)
        my_set.emplace(v[i]);
    return std::vector<unsigned int>(my_set.begin(), my_set.end());
}

std::vector<unsigned int> get_possible_quantities(const std::vector<unsigned int> &sizes,
                                                  const unsigned int &min_allowed,
                                                  const unsigned int &max_allowed) {
    if (sizes.size() > 3) {
        std::vector<unsigned int> max_per_size(7, 0);
        for (unsigned int i = 0; i < sizes.size(); i++)
            max_per_size[sizes[i] - 2]++;
        std::vector<bool> quantity_is_possible = cleverly_get_possible_quantities(max_allowed, max_per_size);

        std::vector<unsigned int> res(0);
        // push them in this order in order to make the search for best solutions faster.
        for (int i = max_allowed; i >= int(min_allowed); i--)
            if (quantity_is_possible[i])
                res.push_back(i);
        return res;
    }
    unsigned long n = pow(2, sizes.size());
    std::vector<unsigned int> res (n, 0);
    for (unsigned int i = 0; i < n; i++){
        std::vector<bool> mask = get_binary_rep(i, sizes.size());
        for (unsigned int j = 0; j < sizes.size(); j++)
            if (mask[j])
                res[i] += sizes[j];
    }
    std::set<unsigned int, greater> my_set;
    for (unsigned int i = 0; i < res.size(); i++)
        if (res[i] <= max_allowed && res[i] >= min_allowed)
            my_set.emplace(res[i]);
    res = std::vector<unsigned int>(my_set.begin(), my_set.end());
    return res;
}

std::vector<bool>
cleverly_get_possible_quantities(const unsigned int &max_allowed, const std::vector<unsigned int> &max_per_size) {
    std::vector<bool> quantity_is_possible(max_allowed + 1, false);
    quantity_is_possible[0] = true;
    for (unsigned int i = 0; i < 7; i++) {
        if (max_per_size[i] == 0) continue;
        for (int m = max_allowed - 2; m >= 0; m--) {
            if (!quantity_is_possible[m]) continue;
            for (unsigned int s = 1; s <= max_per_size[i]; s++) {
                unsigned int next_pos = m + s*(2 + i);
                if (next_pos > max_allowed || quantity_is_possible[next_pos]) break;
                quantity_is_possible[next_pos] = true;
            }
        }
    }
    return quantity_is_possible;
}

std::vector<bool> random_selection(const unsigned int &p, const unsigned int &n) {
    std::vector<bool> res (n, false);
    for (unsigned int i = 0; i < p; i++)
        res[i] = true;
    std::random_shuffle(res.begin(), res.end());
    return res;
}

std::vector<bool> random_family_selection(const unsigned int &nb_families) {
    return random_selection(nb_families, NB_FAMILIES);
}

std::vector<bool> random_day_selection(const unsigned int &nb_days) {
    return random_selection(nb_days, NB_DAYS);
}


void get_differences_between_solutions(const std::vector<unsigned int> &sol1, const std::vector<unsigned int> &sol2, std::vector<bool> &is_different) {
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        if (sol1[i] != sol2[i])
            is_different[i] = true;
}
