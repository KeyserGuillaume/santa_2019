#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "tools.h"
#include "Flow.h"
#include "LocalSearch.h"
#include "Presets.h"

std::vector<unsigned int>
get_possible_quantities(const std::vector<unsigned int> &sizes, const unsigned int &max_allowed);
void recursively_find_possible_quantities(const std::vector<unsigned int> &max_per_size,
                                          std::vector<unsigned int> &possible_quantities,
                                          const unsigned int &max_allowed, const unsigned int &current_quantity,
                                          const unsigned int &current_index);
void recursively_find_quantity_decompositions(const std::vector<unsigned int> &max_per_size,
                                              const unsigned int& quantity,
                                              std::vector<std::vector<unsigned int>> &decompositions,
                                              std::vector<unsigned int> &current_decomposition,
                                              const unsigned int &current_quantity,
                                              const unsigned int& current_index,
                                              const bool& one_is_enough = false);
void recursively_find_all_masks_fitting_decomposition(const std::vector<unsigned int> &sizes,
                                                      const std::vector<unsigned int> &decomposition,
                                                      std::vector<unsigned int> &current_decomposition,
                                                      std::vector<unsigned int> &nb_left_to_see_per_size,
                                                      std::vector<std::vector<bool>> &masks,
                                                      std::vector<bool> &current_mask,
                                                      unsigned int &count,
                                                      const unsigned int& current_index,
                                                      const bool& one_is_enough = false);
std::vector<std::vector<bool>> get_masks_fitting_quantity(const std::vector<unsigned int> &sizes, const unsigned int &q, const bool& one_is_enough = false);
std::vector<std::vector<unsigned int>> get_indexes_fitting_quantity(const std::vector<unsigned int> &sizes, const unsigned int &q, const bool& one_is_enough = false);
bool anticipate_on_day(const Presets &presets,
                       const unsigned int &day,
                       const unsigned int &reference_cost,
                       std::vector<unsigned int> &assignations,
                       std::vector<unsigned int> &counter_assignations);
void find_all_equivalent_solutions(const std::vector<std::vector<unsigned int>> &possible_quantities,
                                   const std::vector<float> &possible_costs,
                                   const unsigned int &cost_threshold,
                                   const unsigned int &min_occupancy,
                                   const unsigned int &max_occupancy,
                                   std::vector<std::vector<unsigned int>> &solutions,
                                   std::vector<unsigned int> &current_solution,
                                   const unsigned int &nb_assigned_people,
                                   const float &current_cost,
                                   const unsigned int &next_cost_index);
void find_best_solution(const std::vector<std::vector<unsigned int>> &possible_quantities,
                        const std::vector<float> &possible_costs,
                        const unsigned int& min_occupancy,
                        const unsigned int& max_occupancy,
                        float& best_solution,
                        const unsigned int &nb_assigned_people,
                        const float &current_cost,
                        const unsigned int &next_cost_index);
void find_one_optimal_solution(const std::vector<std::vector<unsigned int>> &possible_quantities,
                               const std::vector<double> &possible_costs,
                               const unsigned int& min_occupancy,
                               const unsigned int& max_occupancy,
                               double& best_solution_cost,
                               std::vector<unsigned int> &best_solution,
                               std::vector<unsigned int> &current_solution,
                               const unsigned int &nb_assigned_people,
                               const float &current_cost,
                               const unsigned int &next_cost_index);
void make_incomplete_counter_assignations(std::vector<std::pair<unsigned int, float>> families_to_consider,
                                          std::vector<std::vector<unsigned int>> corresponding_families,
                                          const std::vector<std::vector<unsigned int>> &solutions,
                                          const unsigned int &day,
                                          std::vector<unsigned int> &counter_assignations);
int get_day_lower_bound(const Presets &presets,
                        const unsigned int &day,
                        const unsigned int& min_allowed,
                        const unsigned int& max_allowed,
                        const std::vector<float> &lambda,
                        std::vector<unsigned int> &solution,
                        const bool& print);
std::vector<float> try_lagrangian_thing(const Presets& presets);
unsigned int get_lagrangian_lb(const Presets &presets, std::vector<float> &lambda, bool &is_primal_feasible,
                               const bool &stop_if_decrease = true);

void solve_day_subpb_with_dp(const Presets &presets,
                             const unsigned int &day,
                             const unsigned int& min_allowed,
                             const unsigned int& max_allowed,
                             const std::vector<float> &lambda,
                             std::vector<std::vector<unsigned int>>& solutions,
                             std::vector<float> &values);



class RLLB{
    std::vector<std::vector<unsigned int>> possible_Ns;
    std::vector<unsigned int> optimal_Ns_indexes;
    std::vector<float> lambda;
    bool _is_primal_feasible;
    std::vector<std::vector<std::vector<unsigned int>>> day_sub_pb_solutions;
    std::vector<std::vector<double>> day_subpb_sol_values;
    float lb;
    std::vector<unsigned int> family_uses_nb;

    void solve_day_subpbs(const Presets &presets);
    void solve_day_subpbs(const Presets &presets, const unsigned int &j, const bool& drop_prev_sol = false);
    void compute_lb(const Presets& presets);
    std::vector<std::vector<double>> compute_DP_lb(const Presets& presets);
    void compute_DP_lb_2(const Presets& presets);
    void compute_feasibility(const Presets& presets);
    bool be_silent = false;
public:
    bool solve_strong_dp = false; // whether we use DP_lb or DP_lb_2
    bool use_dp_for_sub_pbs = true; // how we solve the sub-problems. DP is much faster

    RLLB () {}
    RLLB(const Presets & presets);
    RLLB(const RLLB& rllb);
    RLLB(const Presets & presets, const RLLB& rllb);
    void initialize_properties(const Presets & presets);
    void notify_assignment(const Presets &presets, const unsigned int& i, const unsigned int& k_assign);
    void compute_family_uses_nb();
//    void make_step(const Presets& presets, const float& alpha=1);
    void optimize_lambda(const Presets &presets, const bool &stop_if_decrease = true, const unsigned int& goal = 1000000, const bool& once_only = false);

    float get_lb() const{return lb;}
    bool is_primal_feasible() const {return _is_primal_feasible;}
    std::vector<unsigned int> get_solution(const Presets& presets) const;
    unsigned int suggest_branching_family(const Presets & presets) const;
    unsigned int suggest_best_branching_family(Presets & presets) const;
    void carry_out_tests(Presets & presets);

    std::vector<int> get_partial_solution(const Presets& presets) const;
    void show_lagrangian_stats(const Presets& presets) const;
    void write_lambda(const std::string& filename) const;
    void read_lambda(const Presets &presets, const std::string& filename);
    void compute_true_day_occupancy_bounds(const Presets& presets, const float& upper_bound=68900);

};