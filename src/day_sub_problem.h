#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "tools.h"
#include "Flow.h"
#include "LocalSearch.h"
#include "Presets.h"

std::vector<unsigned int> get_possible_quantities(const std::vector<unsigned int> &sizes);
void recursively_find_possible_quantities(const std::vector<unsigned int> &max_per_size,
                                          std::vector<unsigned int> &possible_quantities,
                                          const unsigned int &current_quantity,
                                          const unsigned int& current_index);
void recursively_find_quantity_decompositions(const std::vector<unsigned int> &max_per_size,
                                              const unsigned int& quantity,
                                              std::vector<std::vector<unsigned int>> &decompositions,
                                              std::vector<unsigned int> &current_decomposition,
                                              const unsigned int &current_quantity,
                                              const unsigned int& current_index);
void recursively_find_all_masks_fitting_decomposition(const std::vector<unsigned int> &sizes,
                                                      const std::vector<unsigned int> &decomposition,
                                                      std::vector<unsigned int> &current_decomposition,
                                                      std::vector<unsigned int> &nb_left_to_see_per_size,
                                                      std::vector<std::vector<bool>> &masks,
                                                      std::vector<bool> &current_mask,
                                                      unsigned int &count,
                                                      const unsigned int& current_index);
std::vector<std::vector<bool>> get_masks_fitting_quantity(const std::vector<unsigned int> &sizes, const unsigned int &q);
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
void make_incomplete_counter_assignations(std::vector<std::pair<unsigned int, float>> families_to_consider,
                                          std::vector<std::vector<unsigned int>> corresponding_families,
                                          const std::vector<std::vector<unsigned int>> &solutions,
                                          const unsigned int &day,
                                          std::vector<unsigned int> &counter_assignations);
