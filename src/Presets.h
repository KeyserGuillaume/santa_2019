#pragma once

#include <iostream>
#include "constants.h"
#include "tools.h"
#include "day_cost_lower_bound.h"

class Presets {
    std::vector<std::vector<unsigned int>> *family_data;
    std::vector<preset> presets;
    std::vector<bool> is_already_assigned;
    std::vector<unsigned int> preset_occupancy, preset_cardinal, occupancy_lb, occupancy_ub, occupancy_explicit_ub;
    std::vector<uint_pair> prescribed_occupancy_ub, prescribed_occupancy_lb;
    std::vector<uint_pair> bottleneck_bounds;
    std::vector<std::vector<unsigned int>> k_fold_bottleneck_bounds;
    std::vector<double> day_costs_lb;
    double day_cost_DP_lb = 0;
    unsigned int presets_costs = 0;
    double day_cost_lower_bound = 0;
    bool is_feasible_ = true;
    unsigned int nb_assignments = 0;
    // too annoying to maintain, not used
    //unsigned int nb_forbidden_assignments = 0;
    std::vector<uint_pair> sorted_families;


//    void compute_occupancy_lb(const unsigned int &i);
//    void compute_occupancy_ub(const unsigned int &i);
    void compute_occupancy_bounds();
    void compute_bottleneck_bounds();
    void compute_k_fold_bottleneck_ub();
    void compute_day_cost_DP_lb();
    double get_day_cost_lb_(const unsigned int &i);
    double get_additional_day_cost_lb(const unsigned int &i);

public:

    // move up those two once they've been optimized....?
    void compute_all_bounds(const bool & including_costs = true);
    void compute_feasibility();

    Presets () {}
    Presets& operator=(const Presets& other);
    Presets (std::vector<std::vector<unsigned int>> *family_data);
    Presets (std::vector<std::vector<unsigned int>> *family_data, const std::string & filename);
    const preset& operator[](const unsigned int &i) const{return presets[i];}
    void assign_family(const unsigned int &i, const unsigned int &k, const bool &compute_stuff = true, const bool &k_is_relative = true);
    void deassign_family(const unsigned int &i, const bool &compute_bounds = true);
    void forbid_assignment(const unsigned int &i, const unsigned int &k, const bool &compute_stuff = true, const bool &k_is_relative = true);
    void enable_assignment(const unsigned int &i, const unsigned int &k, const bool &compute_stuff = true, const bool &k_is_relative = true);
    void prescribe_occupancy_ub(const unsigned int &i_day, const unsigned int &N){prescribed_occupancy_ub.push_back(uint_pair(i_day, N));}
    void prescribe_occupancy_lb(const unsigned int &i_day, const unsigned int &N){prescribed_occupancy_lb.push_back(uint_pair(i_day, N));}
    void pop_last_occupancy_ub_prescription(){prescribed_occupancy_ub.pop_back();}
    void pop_last_occupancy_lb_prescription(){prescribed_occupancy_lb.pop_back();}

    std::vector<unsigned int> get_solution() const;
    void write_solution(const std::string &filename) const;
    std::string hash() const;

    preset get_preset(const unsigned int& i) const { return presets[i]; }
    void set_preset(const unsigned int& i, const preset& p, const bool& compute_stuff = true);

    unsigned int get_family_data(const unsigned int &i, const unsigned int &k) const{return (*family_data)[i][k];}
    unsigned int get_family_size(const unsigned int &i) const{return (*family_data)[i][NB_CHOICES];}
    bool is_feasible() const{return is_feasible_;}
    bool is_family_alr_assigned(const unsigned int i) const{return is_already_assigned[i];}
    unsigned int get_nb_assignments() const{return nb_assignments;}
    unsigned int get_presets_occupancy(const unsigned int &i) const{return preset_occupancy[i];}
    unsigned int get_presets_cardinal(const unsigned int &i) const{return preset_cardinal[i];}
    unsigned int get_occupancy_lb(const unsigned int &i) const{return occupancy_lb[i];}
    unsigned int get_occupancy_ub(const unsigned int &i) const{return occupancy_ub[i];}
    unsigned int get_bottleneck_lb(const unsigned int &i) const{return bottleneck_bounds[i].first;}
    unsigned int get_bottleneck_ub(const unsigned int &i) const{return bottleneck_bounds[i].second;}
    unsigned int get_bottleneck_ub(const unsigned int &i, const unsigned int &k) const{return k_fold_bottleneck_bounds[i][k];}
    unsigned int get_presets_costs() const{return presets_costs;}
    unsigned int get_presets_costs(const unsigned int& day) const;
    unsigned int get_day_cost_lb() const{return int(floor(day_cost_lower_bound));}
    double get_float_day_cost_lb() const{return day_cost_lower_bound;}
    std::vector<unsigned int> get_largest_unassigned_families() const;
    uint_pair get_bounds_to_branch() const;
    std::vector<unsigned int> get_all_families_assigned_to_day(const unsigned int& i) const;
    unsigned int get_k(const unsigned int& family_index, const unsigned int& day_index) const;
    std::vector<uint_pair> get_assignments_by_default() const;
    unsigned int get_family_max_size_min_choices() const;

    bool should_we_compute_additional_day_costs() const;
    bool should_we_compute_day_costs_with_DP() const;

    std::vector<std::vector<unsigned int>> get_possible_quantities_per_day() const;
    bool are_there_possible_quantities_per_day() const;

    void write_presets(const std::string& filename);
};
