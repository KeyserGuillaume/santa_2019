#pragma once

#include<iostream>
#include "constants.h"
#include "tools.h"


class Presets {
    std::vector<std::vector<unsigned int>> &family_data;
    std::vector<preset> presets;
    std::vector<bool> is_already_assigned;
    std::vector<unsigned int> preset_occupancy, preset_cardinal, occupancy_lower_bounds, occupancy_upper_bounds;
    std::vector<uint_pair> bottleneck_bounds;
    std::vector<std::vector<unsigned int>> k_fold_bottleneck_bounds;
    std::vector<float> day_costs_lower_bounds;
    unsigned int presets_costs = 0;
    float day_cost_lower_bound = 0;
    bool is_feasible_ = true;
    unsigned int nb_assignments = 0;
    unsigned int nb_forbidden_assignments = 0;
    std::vector<uint_pair> sorted_families;


//    void compute_occupancy_lb(const unsigned int &i);
//    void compute_occupancy_ub(const unsigned int &i);
    void compute_occupancy_bounds();
    void compute_bottleneck_bounds();
    void compute_k_fold_bottleneck_ub();
    void compute_day_costs_lb(const unsigned int &i);

public:
    // move up those two once they've been optimized....?
    void compute_all_bounds();
    void compute_feasibility();

    Presets(std::vector<std::vector<unsigned int>> &family_data);
    const preset& operator[](const unsigned int &i) const{return presets[i];}
    void assign_family(const unsigned int &i, const unsigned int &k, const bool &compute_stuff = true);
    void deassign_family(const unsigned int &i, const bool &compute_bounds = true);

    std::vector<unsigned int> get_solution() const;
    void write_solution(const std::string &filename) const;

    unsigned int get_family_data(const unsigned int &i, const unsigned int &k) const{return family_data[i][k];}
    unsigned int get_family_size(const unsigned int &i) const{return family_data[i][NB_CHOICES];}
    bool is_feasible() const{return is_feasible_;}
    bool is_family_alr_assigned(const unsigned int i) const{return is_already_assigned[i];}
    unsigned int get_nb_assignments() const{return nb_assignments;}
    unsigned int get_nb_forbidden_assignments() const{return nb_forbidden_assignments;}
    unsigned int get_presets_occupancy(const unsigned int &i) const{return preset_occupancy[i];}
    unsigned int get_presets_cardinal(const unsigned int &i) const{return preset_cardinal[i];}
    unsigned int get_occupancy_lb(const unsigned int &i) const{return occupancy_lower_bounds[i];}
    unsigned int get_occupancy_ub(const unsigned int &i) const{return occupancy_upper_bounds[i];}
    unsigned int get_bottleneck_lb(const unsigned int &i) const{return bottleneck_bounds[i].first;}
    unsigned int get_bottleneck_ub(const unsigned int &i) const{return bottleneck_bounds[i].second;}
    unsigned int get_bottleneck_ub(const unsigned int &i, const unsigned int &k) const{return k_fold_bottleneck_bounds[i][k];}
    unsigned int get_presets_costs() const{return presets_costs;}
    unsigned int get_day_cost_lower_bound() const{return std::max(4500, int(floor(day_cost_lower_bound)));}
    float get_day_cost_lower_bound(const unsigned int &i) const{return day_costs_lower_bounds[i];}
};
