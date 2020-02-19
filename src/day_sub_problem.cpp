#include "day_sub_problem.h"
#include <limits>

void lagrangian_stats(const Presets &presets, const std::vector<double> &lambda,
                      const std::vector<unsigned int> &family_uses_nb);

//void assignments_costs_only_lb(const Presets &presets, const std::vector<double> &lambda,
//                               const std::vector<std::vector<unsigned int>> &solutions, int &sum_cost_lbs,
//                               std::vector<unsigned int> &family_uses_nb);

void recursively_find_possible_quantities(const std::vector<unsigned int> &max_per_size,
                                          std::vector<unsigned int> &possible_quantities,
                                          const unsigned int &max_allowed, const unsigned int &current_quantity,
                                          const unsigned int &current_index) {
    // 2, 3, 4, 5, 6, 7, 8 make up 7 possible family sizes
    if (current_index == 7) {
        possible_quantities.push_back(current_quantity);
        return;
    }
    for (unsigned int i = 0; i <= max_per_size[current_index] && current_quantity + i * (current_index + 2) <= max_allowed; i++){
        recursively_find_possible_quantities(max_per_size, possible_quantities, max_allowed,
                                             current_quantity + i * (current_index + 2),
                                             current_index + 1);
    }
}

std::vector<unsigned int> get_possible_quantities(const std::vector<unsigned int> &sizes, const unsigned int &max_allowed) {
//    if (sizes.size() > 3) {
        std::vector<unsigned int> max_per_size(7, 0);
        for (unsigned int i = 0; i < sizes.size(); i++)
            max_per_size[sizes[i] - 2]++;
        std::vector<bool> quantity_is_possible = cleverly_get_possible_quantities(max_allowed, max_per_size);

        std::vector<unsigned int> res(0);
        // push then in this order in order to make the search for best solutions faster.
        for (int i = max_allowed; i >= 0; i--)
            if (quantity_is_possible[i])
                res.push_back(i);
        return res;
//    }
//    unsigned long n = pow(2, sizes.size());
//    std::vector<unsigned int> res (n, 0);
//    for (unsigned int i = 0; i < n; i++){
//        std::vector<bool> mask = get_binary_rep(i, sizes.size());
//        for (unsigned int j = 0; j < sizes.size(); j++)
//            if (mask[j])
//                res[i] += sizes[j];
//    }
//    std::set<unsigned int, greater> my_set;
//    for (unsigned int i = 0; i < res.size(); i++)
//        my_set.emplace(res[i]);
//    res = std::vector<unsigned int>(my_set.begin(), my_set.end());
//    return res;
}

void recursively_find_quantity_decompositions(const std::vector<unsigned int> &max_per_size,
                                              const unsigned int& quantity,
                                              std::vector<std::vector<unsigned int>> &decompositions,
                                              std::vector<unsigned int> &current_decomposition,
                                              const unsigned int &current_quantity,
                                              const unsigned int &current_index,
                                              const bool& one_is_enough) {
    if (one_is_enough && decompositions.size() > 0) return;
    if (current_quantity == quantity){
        decompositions.push_back(current_decomposition);
        return;
    }
    if (current_index == 7 || current_quantity > quantity) return;
    unsigned int quantity_upper_bound = current_quantity;
    for (unsigned int i = current_index; i < 7; i++)
        quantity_upper_bound += (i + 2) * max_per_size[i];
    if (quantity_upper_bound < quantity) return;
    for (unsigned int i = 0; i <= max_per_size[current_index]; i++){
        current_decomposition[current_index] = i;
        recursively_find_quantity_decompositions(max_per_size, quantity, decompositions, current_decomposition,
                                                 current_quantity + i * (current_index + 2), current_index + 1, one_is_enough);
        current_decomposition[current_index] = 0;
    }
}

void recursively_find_all_masks_fitting_decomposition(const std::vector<unsigned int> &sizes,
                                                      const std::vector<unsigned int> &decomposition,
                                                      std::vector<unsigned int> &current_decomposition,
                                                      std::vector<unsigned int> &nb_left_to_see_per_size,
                                                      std::vector<std::vector<bool>> &masks,
                                                      std::vector<bool> &current_mask,
                                                      unsigned int &count,
                                                      const unsigned int &current_index,
                                                      const bool& one_is_enough) {
    count++; if (count > 25000000 || (one_is_enough && masks.size() > 0)) return;
    bool decomposition_is_respected = true;
    for (unsigned int i = 0; i < 7; i++) {
        if (current_decomposition[i] > decomposition[i] ||
            current_decomposition[i] + nb_left_to_see_per_size[i] < decomposition[i])
            return;
        if (current_decomposition[i] != decomposition[i])
            decomposition_is_respected = false;
    }
    if (decomposition_is_respected){
        masks.push_back(current_mask);
        return;
    }
    if (current_index == sizes.size()) return;
    nb_left_to_see_per_size[sizes[current_index] - 2]--;
    recursively_find_all_masks_fitting_decomposition(sizes, decomposition, current_decomposition,
                                                     nb_left_to_see_per_size, masks, current_mask, count, current_index + 1, one_is_enough);
    current_mask[current_index] = true;
    current_decomposition[sizes[current_index] - 2]++;
    recursively_find_all_masks_fitting_decomposition(sizes, decomposition, current_decomposition,
                                                     nb_left_to_see_per_size, masks, current_mask, count, current_index + 1, one_is_enough);
    current_mask[current_index] = false;
    current_decomposition[sizes[current_index] - 2]--;
    nb_left_to_see_per_size[sizes[current_index] - 2]++;
}

void recursively_find_all_indexes_fitting_decomposition(const std::vector<unsigned int> &sizes,
                                                        const std::vector<unsigned int> &decomposition,
                                                        std::vector<unsigned int> &current_decomposition,
                                                        std::vector<unsigned int> &nb_left_to_see_per_size,
                                                        std::vector<std::vector<unsigned int>> &indexes,
                                                        std::vector<unsigned int> &current_solution,
                                                        unsigned int &count,
                                                        const unsigned int &current_index,
                                                        const bool& one_is_enough) {
    count++;
    if (count > 25000000 || (one_is_enough && indexes.size() > 0)) return;
    bool decomposition_is_respected = true;
    for (unsigned int i = 0; i < 7; i++) {
        if (current_decomposition[i] > decomposition[i] ||
            current_decomposition[i] + nb_left_to_see_per_size[i] < decomposition[i])
            return;
        if (current_decomposition[i] != decomposition[i])
            decomposition_is_respected = false;
    }
    if (decomposition_is_respected){
        indexes.push_back(current_solution);
        return;
    }
    if (current_index == sizes.size()) return;
    nb_left_to_see_per_size[sizes[current_index] - 2]--;
    recursively_find_all_indexes_fitting_decomposition(sizes, decomposition, current_decomposition,
                                                     nb_left_to_see_per_size, indexes, current_solution, count, current_index + 1, one_is_enough);
    current_solution.push_back(current_index);
    current_decomposition[sizes[current_index] - 2]++;
    recursively_find_all_indexes_fitting_decomposition(sizes, decomposition, current_decomposition,
                                                     nb_left_to_see_per_size, indexes, current_solution, count, current_index + 1, one_is_enough);
    current_solution.pop_back();
    current_decomposition[sizes[current_index] - 2]--;
    nb_left_to_see_per_size[sizes[current_index] - 2]++;
}

std::vector<std::vector<bool>> get_masks_fitting_quantity(const std::vector<unsigned int> &sizes, const unsigned int &q, const bool& one_is_enough) {
    if (sizes.size() > 0) {
        //throw std::invalid_argument("too many sizes");
        std::vector<unsigned int> max_per_size(7, 0);
        for (unsigned int i = 0; i < sizes.size(); i++)
            max_per_size[sizes[i] - 2]++;
        std::vector<std::vector<unsigned int>> possible_decompositions(0);
        std::vector<unsigned int> current_decomposition(7, 0);
        recursively_find_quantity_decompositions(max_per_size, q, possible_decompositions, current_decomposition, 0, 0, one_is_enough);
        std::vector<std::vector<bool>> res(0);
        std::vector<bool> current_mask(sizes.size(), false);
        unsigned int count = 3;
        for (std::vector<unsigned int>& decomposition: possible_decompositions)
            recursively_find_all_masks_fitting_decomposition(sizes, decomposition, current_decomposition, max_per_size, res, current_mask, count, 0, one_is_enough);
        if (count % 10000 == 1)
            return std::vector<std::vector<bool>>(0);
        if (res.size() == 0)
            return res;
        return res;
    }
    std::vector<std::vector<bool>> res(0);
    unsigned long n = pow(2, sizes.size());
    for (unsigned int i = 0; i < n; i++){
        unsigned int sum = 0;
        std::vector<bool> mask = get_binary_rep(i, sizes.size());
        for (unsigned int j = 0; j < sizes.size(); j++)
            if (mask[j])
                sum += sizes[j];
        if (sum == q) {
            res.push_back(mask);
            if (one_is_enough) return res;
        }
    }
    return res;
}

std::vector<std::vector<unsigned int>> get_indexes_fitting_quantity(const std::vector<unsigned int> &sizes, const unsigned int &q, const bool& one_is_enough) {
    if (sizes.size() > 0) {
        std::vector<unsigned int> max_per_size(7, 0);
        for (unsigned int i = 0; i < sizes.size(); i++)
            max_per_size[sizes[i] - 2]++;
        std::vector<std::vector<unsigned int>> possible_decompositions(0);
        std::vector<unsigned int> current_decomposition(7, 0);
        recursively_find_quantity_decompositions(max_per_size, q, possible_decompositions, current_decomposition, 0, 0, one_is_enough);
        std::vector<std::vector<unsigned int>> res(0);
        std::vector<unsigned int> current_solution(0);
        unsigned int count = 3;
        for (std::vector<unsigned int>& decomposition: possible_decompositions)
            recursively_find_all_indexes_fitting_decomposition(sizes, decomposition, current_decomposition, max_per_size, res, current_solution, count, 0, one_is_enough);
        if (count % 100000 == 1)
            return std::vector<std::vector<unsigned int>>(0);
        if (res.size() == 0)
            return res;
        return res;
    }
    std::vector<std::vector<unsigned int>> res(0);
    return res;
}

// for a given day and a certain day cost budget, find all solutions under budget and make assignations and
// counter_assignations based on that.
bool anticipate_on_day(const Presets &presets, const unsigned int &day, const unsigned int &reference_cost,
                       std::vector<unsigned int> &assignations, std::vector<unsigned int> &counter_assignations) {
    std::vector<std::vector<unsigned int>> solutions;
    std::vector<std::pair<unsigned int, double>> families_to_consider(0);
    // sort the families by increasing cost per person
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        for (unsigned int k = 0; k < K_MAX; k++)
            if (presets.get_family_data(i, k) == day && presets[i][k] == ALLOWED)
                families_to_consider.push_back(std::pair<unsigned int, double>(i, double(CONSTANT_COST[k])/double(presets.get_family_size(i)) + MARGINAL_COST[k]));
    sort_by_second(families_to_consider);

    if (families_to_consider.size() == 0)
        return true;

    // group the families by the ~20 levels of cost per person
    std::vector<double> possible_costs(1, families_to_consider[0].second);
    std::vector<std::vector<unsigned int>> possible_sizes(1, std::vector<unsigned int>(0));
    std::vector<std::vector<unsigned int>> corresponding_families(1, std::vector<unsigned int>(0));
    unsigned int current_index = 0;
    for(unsigned int i = 0; i < families_to_consider.size(); i++) {
        if (families_to_consider[i].second != possible_costs[current_index]){
            possible_costs.push_back(families_to_consider[i].second);
            possible_sizes.push_back(std::vector<unsigned int>(1, presets.get_family_size(families_to_consider[i].first)));
            corresponding_families.push_back(std::vector<unsigned int>(1, i));
            current_index++;
        } else {
            possible_sizes[current_index].push_back(presets.get_family_size(families_to_consider[i].first));
            corresponding_families[current_index].push_back(i);
        }
    }

    unsigned int min_allowed = (unsigned int)(std::max(int(presets.get_occupancy_lb(day)) - int(presets.get_presets_occupancy(day)), 0));
    unsigned int max_allowed = presets.get_occupancy_ub(day) - presets.get_presets_occupancy(day);

    // for every level of cost per person, get the possible number of people ( = quantities) that we can assign
    std::vector<std::vector<unsigned int>> unique_possible_quantities;
    for (unsigned int i = 0; i < possible_sizes.size(); i++) {
        unique_possible_quantities.push_back(get_possible_quantities(possible_sizes[i], max_allowed));
    }

    unsigned int presets_costs = presets.get_presets_costs(day);
    if (presets_costs > reference_cost)
        return false;
        //throw std::invalid_argument("This was not supposed to happen");

    // compute the best solution
    double best_solution = 2000;
    find_best_solution(unique_possible_quantities, possible_costs, min_allowed, max_allowed, best_solution, 0, 0, 0);
    std::cout << "Day " << day << ": reference_cost was " << reference_cost << " and best solution was " << best_solution + presets_costs << std::endl;
    unsigned int threshold = (reference_cost > 1.05 * best_solution) ? reference_cost - 15 : reference_cost; // not valid with presets_costs
    threshold = reference_cost - presets_costs; //std::max((unsigned int)(500), threshold);
    std::cout << "We set the threshold to " << threshold + presets_costs << std::endl;

    // compute all the vectors of quantities that give a solution under the reference_cost
    std::vector<unsigned int> current_solution(possible_costs.size(), 0);
    find_all_equivalent_solutions(unique_possible_quantities, possible_costs, threshold, min_allowed, max_allowed,
                                  solutions, current_solution, 0, 0, 0);

    if (solutions.size() > 10000){
        std::cout << "have to abort because too many solutions" << std::endl;
        return true;
    }
    if (solutions.size() == 0){
        std::cout << "have to abort because no solution" << std::endl;
        return false;
    }

    // detect families always used and families never used
    std::vector<unsigned int> never_appears(families_to_consider.size(), true);
    std::vector<unsigned int> always_appears(families_to_consider.size(), true);


    for (unsigned int j = 0; j < possible_costs.size(); j++){
        std::vector<unsigned int> size_never_appears(7, true);
        std::vector<unsigned int> size_always_appears(7, true);

        std::vector<unsigned int> max_per_size(7, 0);
        for (unsigned int s = 0; s < possible_sizes[j].size(); s++)
            max_per_size[possible_sizes[j][s] - 2]++;
        std::vector<unsigned int> unique_solution_values(0);
        
        for (unsigned int i = 0; i < solutions.size(); i++) 
           unique_solution_values.push_back(solutions[i][j]);
        unique_solution_values = get_unique_vector(unique_solution_values);
        
        for (unsigned int m = 0; m < 7; m++){
            if (max_per_size[m] == 0) continue;
            max_per_size[m]--;
            
            // if there is one solution value which we can realize despite this, then size_always_appears[m] = false.
            for (unsigned int l = 0; l < unique_solution_values.size() && size_always_appears[m]; l++){
                std::vector<bool> quantity_is_possible = cleverly_get_possible_quantities(unique_solution_values[l], max_per_size);
                if (quantity_is_possible[unique_solution_values[l]])
                    size_always_appears[m] = false;
            }
            // is there one solution value x such that we can realize the value x - (m + 2) ? if yes, then size_never_appears[m] = false.
            for (unsigned int l = 0; l < unique_solution_values.size() && size_never_appears[m]; l++){
                int target = unique_solution_values[l] - (m + 2);
                if (target < 0) continue;
                std::vector<bool> quantity_is_possible = cleverly_get_possible_quantities(target, max_per_size);
                if (quantity_is_possible[target])
                    size_never_appears[m] = false;
            }

            max_per_size[m]++;
        }
        for (unsigned int s = 0; s < possible_sizes[j].size(); s++) {
            never_appears[corresponding_families[j][s]] = size_never_appears[possible_sizes[j][s] - 2];
            always_appears[corresponding_families[j][s]] = size_always_appears[possible_sizes[j][s] - 2];
        }

    }

    for (unsigned int i = 0; i < families_to_consider.size(); i++) {
        if (always_appears[i]) {
            assignations.push_back(families_to_consider[i].first);
            std::cout << "The family " << families_to_consider[i].first << " is always at day " << day << std::endl;
        }
        else if (never_appears[i]) {
            counter_assignations.push_back(families_to_consider[i].first);
            std::cout << "The family " << families_to_consider[i].first << " is never at day " << day << std::endl;
        }
        else {
            std::cout << "We can't say for family " << families_to_consider[i].first << std::endl;
        }
    }
    return true;
}

void find_all_equivalent_solutions(const std::vector<std::vector<unsigned int>> &possible_quantities,
                                   const std::vector<double> &possible_costs,
                                   const unsigned int &cost_threshold,
                                   const unsigned int &min_occupancy,
                                   const unsigned int &max_occupancy,
                                   std::vector<std::vector<unsigned int>> &solutions,
                                   std::vector<unsigned int> &current_solution,
                                   const unsigned int &nb_assigned_people,
                                   const double &current_cost,
                                   const unsigned int &next_cost_index) {
    if (solutions.size() > 10000) return;
    double nb_needed_people = (nb_assigned_people > min_occupancy) ? 0 : min_occupancy - nb_assigned_people;
    // this must be changed if there are any negative costs. It would likely go faster with the new version anyway.
    if (current_cost + nb_needed_people * possible_costs[next_cost_index] > cost_threshold)
        return;
    if (next_cost_index == possible_costs.size()){
        if (nb_assigned_people < min_occupancy || nb_assigned_people > max_occupancy)
            return;
        //std::cout << current_cost << std::endl;
        solutions.push_back(current_solution);
        return;
    }
    // It's worth continuing so we do the following recursion
    for (unsigned int i = 0; i < possible_quantities[next_cost_index].size(); i++){
        unsigned int q = possible_quantities[next_cost_index][i];
        current_solution[next_cost_index] = q;
        find_all_equivalent_solutions(possible_quantities, possible_costs, cost_threshold, min_occupancy,
                                      max_occupancy, solutions, current_solution, nb_assigned_people + q,
                                      current_cost + q * possible_costs[next_cost_index], next_cost_index + 1);
    }
    current_solution[next_cost_index] = 0;
}

void find_best_solution(const std::vector<std::vector<unsigned int>> &possible_quantities,
                        const std::vector<double> &possible_costs,
                        const unsigned int &min_occupancy,
                        const unsigned int &max_occupancy,
                        double &best_solution,
                        const unsigned int &nb_assigned_people,
                        const double &current_cost,
                        const unsigned int &next_cost_index) {
    double nb_needed_people = (nb_assigned_people > min_occupancy) ? 0 : min_occupancy - nb_assigned_people;
    // this must be changed if there are any negative costs. It would likely go faster with the new version anyway.
    if (current_cost + nb_needed_people * possible_costs[next_cost_index] > best_solution)
        return;
    if (nb_assigned_people > max_occupancy)
        return;
    if (next_cost_index == possible_costs.size()){
        if (nb_assigned_people < min_occupancy)
            return;
        //std::cout << current_cost << std::endl;
        best_solution = std::min(best_solution, current_cost);
        return;
    }
    // It's worth continuing so we do the following recursion
    for (unsigned int i = 0; i < possible_quantities[next_cost_index].size(); i++){
        unsigned int q = possible_quantities[next_cost_index][i];
        find_best_solution(possible_quantities, possible_costs, min_occupancy, max_occupancy, best_solution,
                           nb_assigned_people + q, current_cost + q * possible_costs[next_cost_index], next_cost_index + 1);
    }
}

void make_incomplete_counter_assignations(std::vector<std::pair<unsigned int, double>> families_to_consider,
                                          std::vector<std::vector<unsigned int>> corresponding_families,
                                          const std::vector<std::vector<unsigned int>> &solutions,
                                          const unsigned int& day,
                                          std::vector<unsigned int> &counter_assignations) {
    for (unsigned int i = 0; i < corresponding_families.size(); i++){
        bool used = false;
        for (unsigned int j = 0; j < solutions.size(); j++)
            if (solutions[j][i] > 0)
                used = true;
        if (!used){
            for (unsigned int m = 0; m < corresponding_families[i].size(); m++){
                counter_assignations.push_back(families_to_consider[corresponding_families[i][m]].first);
                std::cout << "We could at least find that the family " << families_to_consider[corresponding_families[i][m]].first << " is never at day " << day << std::endl;
            }
        }
    }
    return;
}

void find_one_optimal_solution(const std::vector<std::vector<unsigned int>> &possible_quantities,
                               const std::vector<double> &possible_costs, const unsigned int &min_occupancy,
                               const unsigned int &max_occupancy, double &best_solution_cost,
                               std::vector<unsigned int> &best_solution, std::vector<unsigned int> &current_solution,
                               const unsigned int &nb_assigned_people, const double &current_cost,
                               const unsigned int &next_cost_index) {
    //unsigned int nb_needed_people = (nb_assigned_people > min_occupancy) ? 0 : min_occupancy - nb_assigned_people;
    unsigned int N = nb_assigned_people;
    double lb = current_cost;
    for (unsigned int i = next_cost_index; i < possible_costs.size(); i++){
        if (N == max_occupancy) break;
        if (possible_costs[i] > 0 && N >= min_occupancy) break;
        unsigned int q;
        if (possible_costs[i] > 0)
            q = std::min(min_occupancy - N, possible_quantities[i][0]);
        else
            q = std::min(max_occupancy - N, possible_quantities[i][0]);
        N += q;
        lb += q * possible_costs[i];
    }
    if (lb > best_solution_cost + 0.000001 || N < min_occupancy || N > max_occupancy)
        return;
    if (nb_assigned_people > max_occupancy)
        return;
    if (next_cost_index == possible_costs.size() || nb_assigned_people == max_occupancy){
        if (nb_assigned_people < min_occupancy)
            return;
        if (current_cost <= best_solution_cost + 0.000001){
            best_solution_cost = current_cost;
            best_solution = std::vector<unsigned int>(current_solution.begin(), current_solution.end());
        }
        return;
    }
    // It's worth continuing so we do the following recursion
    for (unsigned int i = 0; i < possible_quantities[next_cost_index].size(); i++){
        unsigned int q = possible_quantities[next_cost_index][i];
        current_solution[next_cost_index] = q;
        find_one_optimal_solution(possible_quantities, possible_costs, min_occupancy, max_occupancy, best_solution_cost,
                           best_solution, current_solution, nb_assigned_people + q,
                           current_cost + q * possible_costs[next_cost_index], next_cost_index + 1);
    }
    current_solution[next_cost_index] = 0;
}

// solves knapsack-like day subproblem for one occupancy value with a b&b algorithm. It's very efficient when the lambda are
// zero because we factorize the families of identical ratio. For some reason it's still quite good when lambda is integer but when
// lambda is not integer it becomes really slow
int get_day_lower_bound(const Presets &presets, const unsigned int &day, const unsigned int &min_allowed,
                        const unsigned int &max_allowed, const std::vector<double> &lambda,
                        std::vector<unsigned int> &solution, const bool &print) {
    clock_t t0 = clock();
    std::vector<std::pair<unsigned int, double>> families_to_consider(0);

    // if the passed solution already contains some stuff, use it to get a lower bound on the sub-problem optimal solution
    double best_solution_cost = 15000;
    if (solution.size() > 0){
        best_solution_cost = 0;
        unsigned int N = 0;
        for (unsigned int s = 0; s < solution.size(); s++){
            unsigned int i = solution[s];
            unsigned int k = 0;
            for (; k < K_MAX && presets.get_family_data(i, k) != day; k++) {}
            if (k == K_MAX) throw std::logic_error("What's going on here ?");
            if (presets.is_family_alr_assigned(i))  throw std::logic_error("What's going on here ?");
            best_solution_cost += presets.get_family_size(i) * ((double(CONSTANT_COST[k]) + lambda[i]) / double(presets.get_family_size(i)) + MARGINAL_COST[k]);
            N += presets.get_family_size(i);
        }
        if (N > max_allowed || N < min_allowed) throw std::logic_error("What's going on here ?");
    }

    // sort the families by increasing cost per person
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        for (unsigned int k = 0; k < K_MAX; k++)
            if (presets.get_family_data(i, k) == day && presets[i][k] == ALLOWED) {
                double cost = (double(CONSTANT_COST[k]) + lambda[i]) / double(presets.get_family_size(i)) + MARGINAL_COST[k];
                families_to_consider.push_back(std::pair<unsigned int, double>(i, cost));
            }

    sort_by_second(families_to_consider);
    if (print)
        std::cout << "time spent on identifying and sorting the appropriate families: " << clock() - t0 << std::endl;
    t0 = clock();

    if (families_to_consider.size() == 0)
        return presets.get_presets_costs(day);

    // group the families by the ~20 levels of cost per person (~100 levels if many lambdas are non-null)
    std::vector<double> possible_costs(1, families_to_consider[0].second);
    std::vector<std::vector<unsigned int>> possible_sizes(1, std::vector<unsigned int>(0));
    std::vector<std::vector<unsigned int>> corresponding_families(1, std::vector<unsigned int>(0));
    unsigned int current_index = 0;
    for(unsigned int i = 0; i < families_to_consider.size(); i++) {
        if (families_to_consider[i].second != possible_costs[current_index]){
            possible_costs.push_back(families_to_consider[i].second);
            possible_sizes.push_back(std::vector<unsigned int>(1, presets.get_family_size(families_to_consider[i].first)));
            corresponding_families.push_back(std::vector<unsigned int>(1, i));
            current_index++;
        } else {
            possible_sizes[current_index].push_back(presets.get_family_size(families_to_consider[i].first));
            corresponding_families[current_index].push_back(i);
        }
    }

    // for every level of cost per person, get the possible number of people ( = quantities) that we can assign
    std::vector<std::vector<unsigned int>> unique_possible_quantities;
    for (unsigned int i = 0; i < possible_sizes.size(); i++) {
        unique_possible_quantities.push_back(get_possible_quantities(possible_sizes[i], max_allowed));
    }

    if (print) std::cout << "time spent on getting the possible quantities: " << clock() - t0 << std::endl;
    t0 = clock();

        // compute a best solution
    std::vector<unsigned int> best_solution(0);
    std::vector<unsigned int> current_solution(possible_costs.size(), 0);
    find_one_optimal_solution(unique_possible_quantities, possible_costs, min_allowed, max_allowed,
                              best_solution_cost,
                              best_solution, current_solution, 0, 0, 0);

    if (best_solution_cost == 15000) return 15000; // temporary, change it back asap
//        throw std::logic_error("mdk");

    t0 = clock();
    if (print) {
        std::cout << "time spent on getting an optimal solution: " << clock() - t0 << std::endl;
        std::cout << "The day " << day << " has lower bound " << best_solution_cost << std::endl;
    }
    solution.clear();
    for (unsigned int j = 0; j < possible_costs.size(); j++) {
        std::vector<std::vector<bool>> masks = get_masks_fitting_quantity(possible_sizes[j], best_solution[j], true);
        for (unsigned int l = 0; l < masks[0].size(); l++)
            if (masks[0][l])
                solution.push_back(families_to_consider[corresponding_families[j][l]].first);
    }

    if (print) std::cout << "time spent on the rest: " << clock() - t0 << std::endl;
    return best_solution_cost + presets.get_presets_costs(day);
}

// this function uses dynamic programing to solve the knapsack-like day sub-problems in much better time than the b&b way
// because the b&b way starts over for every possible day occupancy value whereas this method does them all at once.
// it is likely that with very close optimal ocupancy bounds it becomes slower than the b&b way but not sure with double lambda
void solve_day_subpb_with_dp(const Presets &presets, const unsigned int &day, const unsigned int &min_allowed,
                             const unsigned int &max_allowed,
                             const std::vector<double> &lambda, std::vector<std::vector<unsigned int>> &solutions,
                             std::vector<double> &values) {
    clock_t t0 = clock();
    std::vector<std::pair<unsigned int, double>> families_to_consider(0);

    // if the passed solution already contains some stuff, use it to get an upper bound on the sub-problem optimal solution
    // that's the idea, but never used it because difficult implementation and not fully convinced it would improve computation time a lot
//    double best_solution_cost = 15000;
//    if (solution.size() > 0){
//        best_solution_cost = 0;
//        unsigned int N = 0;
//        for (unsigned int s = 0; s < solution.size(); s++){
//            unsigned int i = solution[s];
//            unsigned int k = 0;
//            for (; k < K_MAX && presets.get_family_data(i, k) != day; k++) {}
//            if (k == K_MAX) throw std::logic_error("What's going on here ?");
//            if (presets.is_family_alr_assigned(i))  throw std::logic_error("What's going on here ?");
//            best_solution_cost += presets.get_family_size(i) * ((double(CONSTANT_COST[k]) + lambda[i]) / double(presets.get_family_size(i)) + MARGINAL_COST[k]);
//            N += presets.get_family_size(i);
//        }
//        if (N > max_allowed || N < min_allowed) throw std::logic_error("What's going on here ?");
//    }

    // sort the families by increasing cost per person
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        for (unsigned int k = 0; k < K_MAX; k++)
            if (presets.get_family_data(i, k) == day && presets[i][k] == ALLOWED) {
                double cost = (double(CONSTANT_COST[k]) + lambda[i]) / double(presets.get_family_size(i)) + MARGINAL_COST[k];
                families_to_consider.push_back(std::pair<unsigned int, double>(i, cost));
            }

    sort_by_second(families_to_consider);
//    if (print)
//        std::cout << "time spent on identifying and sorting the appropriate families: " << clock() - t0 << std::endl;
//    t0 = clock();

    if (families_to_consider.size() == 0) {
        values[0] = presets.get_presets_costs(day);
        solutions.clear();
        solutions.push_back(std::vector<unsigned int>(0));
        return;
    }

    // group the families by the ~20 levels of cost per person (between 30 and 200 levels if many lambdas are non-null)
    // (day 0 can easily reach 300 levels)
    std::vector<double> possible_costs(1, families_to_consider[0].second);
    std::vector<std::vector<unsigned int>> possible_sizes(1, std::vector<unsigned int>(0));
    std::vector<std::vector<unsigned int>> corresponding_families(1, std::vector<unsigned int>(0));
    unsigned int current_index = 0;
    for(unsigned int i = 0; i < families_to_consider.size(); i++) {
        if (families_to_consider[i].second != possible_costs[current_index]){
            possible_costs.push_back(families_to_consider[i].second);
            possible_sizes.push_back(std::vector<unsigned int>(1, presets.get_family_size(families_to_consider[i].first)));
            corresponding_families.push_back(std::vector<unsigned int>(1, i));
            current_index++;
        } else {
            possible_sizes[current_index].push_back(presets.get_family_size(families_to_consider[i].first));
            corresponding_families[current_index].push_back(i);
        }
    }

//    std::cout << possible_costs.size() << " " << max_allowed - min_allowed + 1 << " ";
//    std::cout << clock() - t0 << " ";
//    t0 = clock();
    // for every level of cost per person, get the possible number of people ( = quantities) that we can assign
    std::vector<std::vector<unsigned int>> unique_possible_quantities;
    for (unsigned int i = 0; i < possible_sizes.size(); i++) {
        unique_possible_quantities.push_back(get_possible_quantities(possible_sizes[i], max_allowed));
    }

    // value_functions[s][n] will be the minimal cost to fill the day with n people
    // using only the s first cost levels
    std::vector<std::vector<double>> value_functions(possible_costs.size(), std::vector<double>(MAX_NB_PEOPLE_PER_DAY + 1, 15000));
    // optimal_quanttities[s][n] will be the quantity of people that should be set
    // at cost level s for n people
    std::vector<std::vector<unsigned int>> optimal_quantities(possible_costs.size(), std::vector<unsigned int>(MAX_NB_PEOPLE_PER_DAY + 1));
    for (unsigned int m = 0; m < unique_possible_quantities[0].size(); m++) {
        unsigned int q = unique_possible_quantities[0][m];
        value_functions[0][q] = q * possible_costs[0];
        optimal_quantities[0][q] = q;
    }
    // above was initialization, below is computation of the value functions
    for (unsigned int s = 1; s < unique_possible_quantities.size(); s++){
        for (unsigned int N = 0; N <= max_allowed; N++){
            double cost, min_cost = 15000;
            for (int m = unique_possible_quantities[s].size() - 1; m >= 0; m--) {
                unsigned int q = unique_possible_quantities[s][m];
                if (q > N) break;
                cost = possible_costs[s] * q + value_functions[s - 1][N - q];
                if (cost < min_cost){
                    min_cost = cost;
                    optimal_quantities[s][N] = q;
                }
            }
            value_functions[s][N] = min_cost;
        }
    }

//    std::cout << clock() - t0 << " ";
//    t0 = clock();
    // here we build the solutions
    std::vector<std::vector<std::vector<unsigned int>>> memoization_list(possible_costs.size(), std::vector<std::vector<unsigned int>> (MAX_NB_PEOPLE_PER_DAY+1, std::vector<unsigned int>(0)));
    for (unsigned int N = min_allowed; N <= max_allowed; N++){
        solutions[N].clear();
        unsigned int n = N;
        for (int s = possible_costs.size() - 1; s >= 0; s--) {
            unsigned int &q = optimal_quantities[s][n];
            if (q != 0) {
                //the following test determines whether a family set matching q on cost level s was previously computed
                // in this block of code, q > 0 so we need at least 1 family thus this test is correct
                if (memoization_list[s][q].size() == 0)
                    memoization_list[s][q] = get_indexes_fitting_quantity(possible_sizes[s], q, true)[0];
                for (unsigned int l = 0; l < memoization_list[s][q].size(); l++)
                        solutions[N].push_back(families_to_consider[corresponding_families[s][memoization_list[s][q][l]]].first);
            }
            if (s >= 1)
                n -= q;
        }
    }

//    std::cout << clock() - t0 << "; ";
    // this operation's complexity is NB_FAMILIES * K_MAX !!
    unsigned int preset_cost = presets.get_presets_costs(day);
    for (unsigned int N = min_allowed; N <= max_allowed; N++)
        values[N] = value_functions[possible_costs.size() - 1][N] + preset_cost;

//    verify the results with the b&b, single-N version
//    for (unsigned int N = min_allowed; N <= max_allowed; N++) {
//        std::vector<unsigned int> solution(0);
//        double ref = get_day_lower_bound(presets, day, N, N, lambda, solution, false);
//        double comp = values[N];
//        if (comp != ref && ref != 15000)
//            throw std::logic_error("lsmjfr");
//    }
    return;
}

void lagrangian_stats(const Presets &presets, const std::vector<double> &lambda,
                      const std::vector<unsigned int> &family_uses_nb) {
    std::vector<unsigned int> histo_family_uses_nb(5, 0);
    std::vector<unsigned int> histo_family_never_used_sizes(7, 0);
    std::vector<unsigned int> histo_family_used_twice_sizes(7, 0);
    std::vector<unsigned int> histo_family_never_used_lambdas(10, 0);
    std::vector<unsigned int> histo_family_used_twice_lambdas(10, 0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        if (presets.is_family_alr_assigned(i)) continue;
        histo_family_uses_nb[family_uses_nb[i]]++;
        if (family_uses_nb[i] == 0) {
            histo_family_never_used_sizes[presets.get_family_size(i) - 2]++;
            histo_family_never_used_lambdas[(unsigned int)((lambda[i] + 200) / 40)]++;
        }
        if (family_uses_nb[i] == 2) {
            histo_family_used_twice_sizes[presets.get_family_size(i) - 2]++;
            histo_family_used_twice_lambdas[(unsigned int) ((lambda[i] + 200) / 40)]++;
        }
    }

    std::cout << "stats on family uses nb" << std::endl;
    print_nicely(std::vector<unsigned int>({0, 1, 2, 3, 4}), 5);
    print_nicely(histo_family_uses_nb, 5);
    std::cout << "stats on families never used : what sizes ?" << std::endl;
    print_nicely(std::vector<unsigned int>({2, 3, 4, 5, 6, 7, 8}), 4);
    print_nicely(histo_family_never_used_sizes, 4);
//    std::cout << "stats on families never used : what lambdas ?" << std::endl;
//    print_nicely(histo_family_never_used_lambdas, 5);
//    std::cout << "stats on families used twice : what sizes ?" << std::endl;
//    print_nicely(std::vector<unsigned int>({2, 3, 4, 5, 6, 7, 8}), 4);
//    print_nicely(histo_family_used_twice_sizes, 4);
//    std::cout << "stats on families never twice : what lambdas ?" << std::endl;
//    print_nicely(histo_family_used_twice_lambdas, 5);
}



RLLB::RLLB(const Presets & presets) {
    initialize_properties(presets);
    solve_day_subpbs(presets);
    compute_lb(presets);
//    lagrangian_stats(presets, lambda, family_uses_nb);
}

RLLB::RLLB(const Presets &presets, const std::vector<double> &lambda) {
    initialize_properties(presets);
    this->lambda = lambda;
    solve_day_subpbs(presets);
    compute_lb(presets);
}

RLLB::RLLB(const RLLB &rllb) {
    possible_Ns = rllb.possible_Ns;
    optimal_Ns_indexes = rllb.optimal_Ns_indexes;
    lambda = rllb.lambda;
    _is_primal_feasible = rllb._is_primal_feasible;
    day_sub_pb_solutions = rllb.day_sub_pb_solutions;
    day_subpb_sol_values = rllb.day_subpb_sol_values;
    lb = rllb.lb;
    family_uses_nb = rllb.family_uses_nb;
}

RLLB::RLLB(const Presets &presets, const RLLB &rllb) {
    initialize_properties(presets);
    lambda = rllb.lambda;
    solve_day_subpbs(presets);
    compute_lb(presets);
}

void RLLB::initialize_properties(const Presets &presets) {
    lambda = std::vector<double>(NB_FAMILIES, 0);
    possible_Ns = presets.get_possible_quantities_per_day();
    _is_primal_feasible = false;
    day_sub_pb_solutions = std::vector<std::vector<std::vector<unsigned int>>> (NB_DAYS, std::vector<std::vector<unsigned int>>(0));
    day_subpb_sol_values = std::vector<std::vector<double>> (NB_DAYS, std::vector<double>(0));
    for (unsigned int j = 0; j < NB_DAYS; j++) {
        day_sub_pb_solutions[j] = std::vector<std::vector<unsigned int>>(possible_Ns[j].size(), std::vector<unsigned int>(0));
        day_subpb_sol_values[j] = std::vector<double>(possible_Ns[j].size(), 0);
    }
}

void RLLB::solve_day_subpbs(const Presets &presets){// solve all the day subproblems
    for (unsigned int j = 0; j < NB_DAYS; j++) {
        solve_day_subpbs(presets, j);
    }
}

void RLLB::solve_day_subpbs(const Presets &presets, const unsigned int &j, const bool& drop_prev_sol){
    if (!use_dp_for_sub_pbs) {
        for (unsigned int n = 0; n < possible_Ns[j].size(); n++) {
            unsigned int N = possible_Ns[j][n] - presets.get_presets_occupancy(j);
            if (drop_prev_sol) day_sub_pb_solutions[j][n] = std::vector<unsigned int>(0);
            day_subpb_sol_values[j][n] = get_day_lower_bound(presets, j, N, N, lambda, day_sub_pb_solutions[j][n], false);
        }
    } else {
        unsigned int min = 1000, max = 0;
        for (unsigned int n = 0; n < possible_Ns[j].size(); n++) {
            min = std::min(min, possible_Ns[j][n]);
            max = std::max(max, possible_Ns[j][n]);
        }
        min -=  presets.get_presets_occupancy(j);
        max -=  presets.get_presets_occupancy(j);
        std::vector<std::vector<unsigned int>> solutions(MAX_NB_PEOPLE_PER_DAY + 1, std::vector<unsigned int>(0));
//        if (!drop_prev_sol){
//            code something to pass the previous solutions as argument
//        }
        std::vector<double> values(MAX_NB_PEOPLE_PER_DAY + 1, 15000);
        solve_day_subpb_with_dp(presets, j, min, max, lambda, solutions, values);
        for (unsigned int n = 0; n < possible_Ns[j].size(); n++) {
            unsigned int N = possible_Ns[j][n] - presets.get_presets_occupancy(j);
            day_subpb_sol_values[j][n] = values[N];
            day_sub_pb_solutions[j][n] = solutions[N];
        }
    }
}

// this function must be called whenever the presets are modified because otherwise the sub-problem values are
// false and so is the lower bound.
void RLLB::notify_assignment(const Presets &presets, const unsigned int &i, const unsigned int& k_assign) {
    lambda[i] = 0;
    std::vector<std::vector<unsigned int>> old_possible_Ns = possible_Ns;
    possible_Ns = presets.get_possible_quantities_per_day();
    // we recompute for the days on which this family could have been assigned.
    // in fact for the k such that k != k_assign and the N such that the solution
    // does not use family i, we do not need to redo the subpb optimization.
    for (unsigned int k = 0; k < K_MAX; k++){
        unsigned int j = presets.get_family_data(i, k);
        day_sub_pb_solutions[j] = std::vector<std::vector<unsigned int>>(possible_Ns[j].size(), std::vector<unsigned int>(0));
        day_subpb_sol_values[j] = std::vector<double>(possible_Ns[j].size(), 0);
        solve_day_subpbs(presets, j, true);
    }
    // because of the quirks of presets, other days could have had their possible quantities changed (via occupancy_lb
    // and occupancy_ub) but there is no need to redo those computations. We just need to do the reindexing.
    for (unsigned int j = 0; j < NB_DAYS; j++){
        if (possible_Ns[j].size() == old_possible_Ns[j].size()) continue;
        bool was_recomputed = false;
        unsigned int k = 0;
        while (k < K_MAX && !was_recomputed) {
            if (presets.get_family_data(i, k) == j)
                was_recomputed = true;
            k++;
        }
        if (was_recomputed) continue;
        // the day j has fewer possible quantities than before, the solutions are still optimal but must be reindexed.
        std::vector<double> old_sol_values = day_subpb_sol_values[j];
        std::vector<std::vector<unsigned int>> old_solutions = day_sub_pb_solutions[j];
        unsigned int m_old = 0, m_new = 0;
        while (m_new < possible_Ns[j].size()){
            // probably most frequent case: the quantities correspond so copy the values
            if (old_possible_Ns[j][m_old] == possible_Ns[j][m_new]) {
                day_subpb_sol_values[j][m_new] = old_sol_values[m_old];
                day_sub_pb_solutions[j][m_new] = old_solutions [m_old];
                m_new++;
            }
            m_old++;
        }
        for (; m_new < m_old; m_new++){
            day_sub_pb_solutions[j].pop_back();
            day_subpb_sol_values[j].pop_back();
        }
    }
    // that is not always necessary but complicated enough already
    compute_lb(presets);
}

void RLLB::compute_feasibility(const Presets &presets) {
    _is_primal_feasible = true;
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        if (!presets.is_family_alr_assigned(i) && family_uses_nb[i] != 1)
            _is_primal_feasible = false;
}

void RLLB::compute_lb(const Presets &presets) {
    // we assume the day subproblems have already been solved
    if (solve_strong_dp)
        compute_DP_lb_2(presets);
    else
        compute_DP_lb(presets);
    return;
}

// this method uses dynamic programing to get a lower bound from the subproblem optimal values
// the value functions are returned.
std::vector<std::vector<double>> RLLB::compute_DP_lb(const Presets& presets) {
    std::vector<std::vector<double>> value_functions(NB_DAYS, std::vector<double>(0));
    std::vector<std::vector<unsigned int>> argmin_value_functions(NB_DAYS, std::vector<unsigned int>(0));
    std::vector<unsigned int> &last_poss_qties = possible_Ns[NB_DAYS - 1];
    for (unsigned int m = 0; m < last_poss_qties.size(); m++) {
        value_functions[NB_DAYS - 1].push_back(
                get_day_cost(last_poss_qties[m], last_poss_qties[m]) + day_subpb_sol_values[NB_DAYS - 1][m]);
    }
    for (int i = NB_DAYS - 2; i >= 0; i--) {
        for (unsigned int n = 0; n < possible_Ns[i].size(); n++) {
            double mini = std::numeric_limits<double>::max() , current;
            unsigned int argmin;
            for (unsigned int m = 0; m < possible_Ns[i + 1].size(); m++) {
                current = value_functions[i + 1][m] +
                          get_day_cost(possible_Ns[i][n], possible_Ns[i + 1][m]) +
                          day_subpb_sol_values[i][n];
                if (current < mini) {
                    mini = current;
                    argmin = m;
                }
            }
            if (mini == std::numeric_limits<double>::max() ) throw std::logic_error("dmlke ?");
            value_functions[i].push_back(mini);
            argmin_value_functions[i + 1].push_back(argmin);
        }
    }

    double lb_psum_lambda = 100000000;
    unsigned int argmin;
    for (unsigned int m = 0; m < possible_Ns[0].size(); m++) {
        if (value_functions[0][m] < lb_psum_lambda) {
            lb_psum_lambda = value_functions[0][m];
            argmin = m;
        }
    }
    if (lb_psum_lambda == 10000000) throw std::logic_error("dmlke ?");

    double sum_lambda = 0;
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        if (!presets.is_family_alr_assigned(i))
            sum_lambda += lambda[i];
    lb = lb_psum_lambda - sum_lambda;

    std::vector<unsigned int> optimal_Ns(NB_DAYS, 0);
    optimal_Ns_indexes = std::vector<unsigned int>(NB_DAYS, 0);
    optimal_Ns[0] = possible_Ns[0][argmin];
    optimal_Ns_indexes[0] = argmin;
    for (unsigned int i = 1; i < NB_DAYS; i++) {
        argmin = argmin_value_functions[i][argmin];
        optimal_Ns[i] = possible_Ns[i][argmin];
        optimal_Ns_indexes[i] = argmin;
    }

    unsigned int c = 0;
    for (unsigned int j = 0; j < NB_DAYS; j++)
        c += optimal_Ns[j];
    if (!be_silent) std::cout << "current lb: " << lb << "  " << c  << std::endl;
//    print_nicely(optimal_Ns);

    compute_family_uses_nb();
    compute_feasibility(presets);

    return value_functions;
}

// this function improves the lower bound but the computation time is too great
// maybe there is a way to improve it but I don't know it
// We do DP to find an optimal day occupancy trajectory but with the added constraint that
// the sum of the occupancies is set to the sum of the family sizes
void RLLB::compute_DP_lb_2(const Presets& presets) {
    // version where we use the contraint \sum_j N_j = \sum_i n_i.
    // it works but it really takes too much time
    std::vector<unsigned int> Nmax(NB_DAYS, MIN_NB_PEOPLE_PER_DAY);
    std::vector<unsigned int> Nmin(NB_DAYS, MAX_NB_PEOPLE_PER_DAY);
    for (unsigned int i = 0; i < NB_DAYS; i++){
        for (unsigned int j = 0; j < possible_Ns[i].size(); j++){
            Nmax[i] = std::max(Nmax[i], possible_Ns[i][j]);
            Nmin[i] = std::min(Nmin[i], possible_Ns[i][j]);
        }
    }

    // those funny lines of code provide lower and upper bounds on the M for which value_function[t][n][M - Mmin]
    // should be computed
    unsigned int sum_ni = 0;
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        sum_ni += presets.get_family_size(i);
    std::vector<unsigned int> Mmin(NB_DAYS);
    std::vector<unsigned int> Mmax(NB_DAYS);
    for (unsigned int i = 0; i < NB_DAYS; i++){
        unsigned int s1 = 0;
        unsigned int s2 = 0;
        unsigned int s3 = 0;
        unsigned int s4 = 0;
        for (unsigned int j = 0; j < i; j++) {
            s1 += Nmin[j];
            s2 += Nmax[j];
        }
        for (unsigned int j = i; j < NB_DAYS; j++) {
            s3 += Nmin[j];
            s4 += Nmax[j];
        }
        Mmin[i] = std::max(int(s3), int(sum_ni - s2));
        Mmax[i] = std::min(int(s4), int(sum_ni - s1));
    }

    // value_function[t][n][M - Mmin],
    // 0       <= t <= 99,
    // corresponds to nth element in possible_Ns[t],
    // Mmin    <= M <= Mmax.
    std::vector<std::vector<std::vector<double>>> value_functions(NB_DAYS, std::vector<std::vector<double>>(0));
    std::vector<std::vector<std::vector<unsigned int>>> argmin_value_functions(NB_DAYS, std::vector<std::vector<unsigned int>>(0));
    for (unsigned int m = 0; m < possible_Ns[NB_DAYS - 1].size(); m++){
        unsigned int q = possible_Ns[NB_DAYS - 1][m];
        double cost = get_day_cost(q, q);
        value_functions[NB_DAYS - 1].push_back(std::vector<double>(Mmax[NB_DAYS - 1] - Mmin[NB_DAYS - 1] + 1, 100000));
        value_functions[NB_DAYS - 1][m][q - Mmin[NB_DAYS - 1]] = cost + day_subpb_sol_values[NB_DAYS - 1][m];
    }
    for (int i = NB_DAYS - 2; i >= 0; i--){
        for (unsigned int n = 0; n < possible_Ns[i].size(); n++) {
            value_functions[i].push_back(std::vector<double>(Mmax[i] - Mmin[i] + 1, 100000));
            argmin_value_functions[i + 1].push_back(std::vector<unsigned int>(Mmax[i] - Mmin[i] + 1));
            for (unsigned int M = Mmin[i]; M <= Mmax[i]; M++) {
                unsigned int q = possible_Ns[i][n];
                if (M - q < Mmin[i + 1] || M - q > Mmax[i + 1])
                    continue;
                double mini = 1000000, current;
                unsigned int argmin;
                for (unsigned int m = 0; m < possible_Ns[i + 1].size(); m++) {
                    current = value_functions[i + 1][m][M - q - Mmin[i + 1]] +
                              get_day_cost(q, possible_Ns[i + 1][m]) +
                              day_subpb_sol_values[i][n];
                    if (current < mini) {
                        mini = current;
                        argmin = m;
                    }
                }
                value_functions[i][n][M - Mmin[i]] = mini;
                argmin_value_functions[i + 1][n][M - Mmin[i]] = argmin;
            }
        }
    }

    double lb_psum_lambda = 1000000;
    unsigned int argmin;
    for (unsigned int m = 0; m < possible_Ns[0].size(); m++) {
        if (value_functions[0][m][0] < lb_psum_lambda){
            lb_psum_lambda = value_functions[0][m][0];
            argmin = m;
        }
    }

    double sum_lambda = 0;
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        if (!presets.is_family_alr_assigned(i))
            sum_lambda += lambda[i];
    lb = lb_psum_lambda - sum_lambda;
    if (!be_silent) std::cout << lb << std::endl;

    std::vector<unsigned int> optimal_Ns(NB_DAYS, 0);
    optimal_Ns_indexes = std::vector<unsigned int>(NB_DAYS, 0);
    optimal_Ns[0] = possible_Ns[0][argmin];
    optimal_Ns_indexes[0] = argmin;
    unsigned int M = sum_ni;
    for (unsigned int i = 1; i < NB_DAYS; i++){
        unsigned int prev_argmin = argmin;
        argmin = argmin_value_functions[i][argmin][M - Mmin[i - 1]];
        M -= possible_Ns[i - 1][prev_argmin];
        optimal_Ns[i] = possible_Ns[i][argmin];
        optimal_Ns_indexes[i] = argmin;
    }

    unsigned int c = 0;
    for (unsigned int j = 0; j < NB_DAYS; j++)
        c += optimal_Ns[j];
//    std::cout << lb << " " << c << std::endl;

    compute_family_uses_nb();
    compute_feasibility(presets);
}

void RLLB::compute_family_uses_nb() {
    family_uses_nb = std::vector<unsigned int>(NB_FAMILIES, 0);
    for (unsigned int j = 0; j < NB_DAYS; j++){
        std::vector<unsigned int> &day_optimal_solution = day_sub_pb_solutions[j][optimal_Ns_indexes[j]];
        for (unsigned int l = 0; l < day_optimal_solution.size(); l++)
            family_uses_nb[day_optimal_solution[l]]++;
    }
}

// Implements a subgradient algorithm to optimize the lagrangian relaxation
// The step alpha is important, we can manually set it (toggle code comments) to constant ( = 1 ), or
// decreasing in a sum_square converging and sum diverging alpha_t = (1/(a + bt))
// the second option gives the best results but it was quicker to use the first for
// the greedy oriented search heuristic
void RLLB::optimize_lambda(const Presets &presets, const bool &stop_if_decrease, const unsigned int& goal, const bool& once_only) {
    // never do this if we are primal-feasible
    if (_is_primal_feasible) throw std::logic_error("why do this ?");
    unsigned int time_since_last_improvement = 0;
    double best_lb = lb;
    std::vector<double> best_lambda = lambda;
    std::vector<std::vector<double>> best_sol_values = day_subpb_sol_values;
    std::vector<std::vector<std::vector<unsigned int>>> best_solutions = day_sub_pb_solutions;
    for (unsigned int t = 0; t < 600; t++) {
        // modify the lambdas with subgradient method (constant step alpha)
        // small detail: the weird commented code corresponded to
        // making sure the lambda of family 0 is always equal to 0
        // to avoid having all the lambdas rush off to -inf or +inf.
//        double alpha = 1/(1. + 20.*t/600.);
        double alpha = 1/(5 + 15.*t/60.);
//        double alpha = 1;
        if (t > 0 && once_only) break;
        std::vector<bool> need_recomputing(NB_DAYS, false);
        for (unsigned int i = 0; i < NB_FAMILIES; i++) {
            if (presets.is_family_alr_assigned(i)) continue;
            double delta;
            if (family_uses_nb[i] > 1) {
                delta = alpha * (family_uses_nb[i] - 1);
            } else if (family_uses_nb[i] < 1) {
                delta = -alpha;
            }
            if (family_uses_nb[i] != 1) {
                for (unsigned int k = 0; k < K_MAX; k++)
                    need_recomputing[presets.get_family_data(i, k)] = true;
//                if (i != 0) {
                    lambda[i] += delta;
//                } else {
//                    for (unsigned int j = 1; j < NB_FAMILIES; j++)
//                        lambda[j] -= delta;
//                }
            }
        }

        // compute the new costs
        for (unsigned int j = 0; j < NB_DAYS; j++)
            if (need_recomputing[j])
                solve_day_subpbs(presets, j);
        compute_lb(presets);
        if (_is_primal_feasible) return;

        if (lb <= best_lb){
            time_since_last_improvement++;
            if (time_since_last_improvement > 5) break;
        } else {
            best_lb = lb;
            best_lambda = lambda;
            best_sol_values = day_subpb_sol_values;
            best_solutions = day_sub_pb_solutions;
            time_since_last_improvement = 0;
        }
        if (lb < best_lb && stop_if_decrease) {
            break;
        }
        if (lb >= goal)
            return;

        if (t % 100 == 0 && !be_silent)
            lagrangian_stats(presets, lambda, family_uses_nb);
    }
    lb = best_lb;
    lambda = best_lambda;
    day_subpb_sol_values = best_sol_values;
    day_sub_pb_solutions = best_solutions;
    return;
}


std::vector<unsigned int> RLLB::get_solution(const Presets& presets) const {
    if (!_is_primal_feasible) throw std::logic_error("not feasible");
    std::vector<unsigned int> res (NB_FAMILIES, 0);
    for (unsigned int j = 0; j < NB_DAYS; j++){
        unsigned int m = optimal_Ns_indexes[j];
        for (unsigned int s = 0; s < day_sub_pb_solutions[j][m].size(); s++)
            res[day_sub_pb_solutions[j][m][s]] = j;
    }
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        if (presets.is_family_alr_assigned(i)) {
            unsigned int k = 0;
            while (k < K_MAX && presets[i][k] != COMPULSORY)
                k++;
            if (k == K_MAX) throw std::logic_error("kjvhbn");
            res[i] = presets.get_family_data(i, k);
        }
    }
    return res;
}

unsigned int RLLB::suggest_branching_family(const Presets & presets) const {
    // the following families are families I know are important from the way
    // they improve the lower bound for the first three and the way
    // they often make me close nodes in the b&b for the others.
    if (!presets.is_family_alr_assigned(4257)) return 4257;
    if (!presets.is_family_alr_assigned(202)) return 202;
    if (!presets.is_family_alr_assigned(3625)) return 3625;
    if (!presets.is_family_alr_assigned(1615)) return 1615;
    if (!presets.is_family_alr_assigned(4471)) return 4471;
    if (!presets.is_family_alr_assigned(874)) return 874;
    if (!presets.is_family_alr_assigned(1723)) return 1723;
    if (!presets.is_family_alr_assigned(557)) return 557;
    std::vector<uint_pair> not_1(0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        if (!presets.is_family_alr_assigned(i) && family_uses_nb[i] != 1)
            not_1.push_back(uint_pair(i, presets.get_family_size(i)));
    sort_by_second(not_1);
    return not_1[not_1.size() - 1].first;
}

// look for the best branching family
// may take a lot of time...
// presets is not const but it should not be modified permanently
unsigned int RLLB::suggest_best_branching_family(Presets &presets) const {
    std::cout << "Looking for the best branching family" << std::endl;
    std::vector<std::pair<unsigned int, double>> quick_lb_per_family(0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        if (presets.is_family_alr_assigned(i) || family_uses_nb[i] == 1)
            continue;
        preset old_preset = presets.get_preset(i);
        RLLB specific_rllb;
        double min_value = 1000000;
        for (unsigned int k_assign = 0; k_assign < K_MAX; k_assign++) {
            specific_rllb = *this;
            specific_rllb.be_silent = true;
            if (old_preset[k_assign] == FORBIDDEN)
                continue;
            presets.assign_family(i, k_assign);
            if (!presets.is_feasible()) {
                std::cout << "unfeasibility" << std::endl;
                continue;
            }
            specific_rllb.notify_assignment(presets, i, k_assign);
            specific_rllb.optimize_lambda(presets, false);
            if (specific_rllb.get_lb() < min_value)
                min_value = specific_rllb.get_lb();
            if (min_value < lb + 1)
                break;
        }
        presets.set_preset(i, old_preset);
        std::cout << "i: " << i << ", min value: " << min_value << std::endl;
        quick_lb_per_family.push_back(std::pair<unsigned int, double>(i, min_value));

        if (min_value > lb + 3)  // win some time, give up on best one
            return i;
    }

    // no need to sort, just iterate ?
    sort_by_second(quick_lb_per_family);

    std::pair<unsigned int, double> best = quick_lb_per_family[quick_lb_per_family.size() - 1];
    std::cout << "Best family found is " << best.first << " with cost of " << best.second << std::endl;
    return best.first;
}


std::vector<std::vector<unsigned int>> RLLB::get_selected_choices(const Presets &presets) const {
    std::vector<std::vector<unsigned int>> selected_choices(NB_FAMILIES, std::vector<unsigned int>(0));

    for (unsigned int j = 0; j < NB_DAYS; j++){
        std::vector<unsigned int> families_of_the_day = day_sub_pb_solutions[j][optimal_Ns_indexes[j]];
        for (unsigned int i: families_of_the_day){
            for (unsigned int k = 0; k < K_MAX; k++){
                if (presets.get_family_data(i, k) == j) {
                    selected_choices[i].push_back(k);
                    break;
                }
            }
        }
    }

    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        if (presets.is_family_alr_assigned(i)) {
            for (unsigned int k = 0; k < K_MAX; k++){
                if (presets[i][k] == COMPULSORY) {
                    selected_choices[i].push_back(k);
                    break;
                }
            }
        }
    }

    return selected_choices;

}

// this one could be refactored using the method get_selected_choices
std::vector<int> RLLB::get_partial_solution(const Presets &presets) const {
    std::vector<int> partial_solution(NB_FAMILIES, -1);
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        if (presets.is_family_alr_assigned(i)) {
            for (unsigned int k = 0; k < K_MAX; k++)
                if (presets[i][k] == COMPULSORY)
                    partial_solution[i] = k;
            continue;
        }
        if (family_uses_nb[i] != 1)
            continue;
        bool found_day = false;
        for (unsigned int k = 0; k < K_MAX && !found_day; k++){
            unsigned int j = presets.get_family_data(i, k);
            std::vector<unsigned int> families_of_the_day = day_sub_pb_solutions[j][optimal_Ns_indexes[j]];
            bool found_family = false;
            for (unsigned int m = 0; m < families_of_the_day.size() && !found_family; m++)
                if (families_of_the_day[m] == i)
                    found_family = true;
            if (found_family){
                partial_solution[i] = k;
                found_day = true;
            }
        }
    }
    return partial_solution;
}

void RLLB::show_lagrangian_stats(const Presets &presets) const {
    lagrangian_stats(presets, lambda, family_uses_nb);
}

void RLLB::write_lambda(const std::string &filename) const {
    std::ofstream write(filename.c_str());
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        write << lambda[i] << std::endl;
    }
}

// a method that delivered me from always waiting for 5 minutes when I run anything
// loads pre-computed lambdas (the best I found)
void RLLB::read_lambda(const Presets &presets, const std::string &filename) {
    initialize_properties(presets);
    lambda = std::vector<double>(0);
    std::ifstream targetFile (filename.c_str());
    if (!targetFile.is_open()) throw std::runtime_error("No targets file found");
    std::string input_line, entity;
    while (std::getline(targetFile, input_line)){
        std::stringstream line(input_line);
        std::getline(line, entity, ',');
        lambda.push_back(std::stof(entity));
    }
    solve_day_subpbs(presets);
    compute_lb(presets);
}

// this identified 3 families which, given the lambda I had, the ~2000 optimal assignations and the optimal occupancy
// bounds, greatly improved my lower bound (from 68848 to 68856 for each of them)
void RLLB::carry_out_tests(Presets &presets) {
    std::vector<double_pair> res(0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        if (presets.is_family_alr_assigned(i) || family_uses_nb[i] == 1)
            continue;
        preset old_preset = presets.get_preset(i);
        RLLB specific_rllb;
        double min_value = 1000000;
        for (unsigned int k_assign = 0; k_assign < K_MAX; k_assign++) {
            specific_rllb = *this;
            if (old_preset[k_assign] == FORBIDDEN)
                continue;
            presets.assign_family(i, k_assign);
            if (!presets.is_feasible()) {
                std::cout << "unfeasibility" << std::endl;
                continue;
            }
            specific_rllb.notify_assignment(presets, i, k_assign);
            specific_rllb.optimize_lambda(presets, false);
            if (specific_rllb.get_lb() < min_value)
                min_value = specific_rllb.get_lb();
        }
        presets.set_preset(i, old_preset);
        std::cout << "lambda is " << lambda[i] << " and min value is " << min_value << std::endl;
        res.push_back(double_pair(lambda[i], min_value));
    }
    sort_by_second(res);
    for (unsigned int m = 0; m < res.size(); m++)
        std::cout << "lambda is " << res[m].first << " and min value is " << res[m].second << std::endl;
}

void RLLB::compute_true_day_occupancy_bounds(const Presets &presets, const double &upper_bound) {
    std::vector<std::vector<double>> lower_bounds = compute_DP_lb(presets);

    // max in the sense that if a solution's days from j to 100 cost more than max_values(N, j) (all costs included),
    // then this solution cannot be optimal. coupled with the lower bounds we get the true occupancy bounds in
    // polynomial time (in fact we get the exact values)
    std::vector<std::vector<double>> max_values(NB_DAYS, std::vector<double>(0));

    double sum_lambda = 0;
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        if (!presets.is_family_alr_assigned(i))
            sum_lambda += lambda[i];

    std::vector<unsigned int> &last_poss_qties = possible_Ns[NB_DAYS - 1];
    for (unsigned int m = 0; m < last_poss_qties.size(); m++)
        max_values[0].push_back(upper_bound + sum_lambda);
    for (unsigned int j = 1; j < NB_DAYS; j++) {
        for (unsigned int n = 0; n < possible_Ns[j].size(); n++) {
            double maxi = -std::numeric_limits<double>::max() , current;
            for (unsigned int m = 0; m < possible_Ns[j - 1].size(); m++) {
                double tmp = get_day_cost(possible_Ns[j - 1][m], possible_Ns[j][n]);
                double tmp2 = day_subpb_sol_values[j - 1][m];
                double tmp3 = max_values[j - 1][m];
                current = max_values[j - 1][m] -
                          get_day_cost(possible_Ns[j - 1][m], possible_Ns[j][n]) -
                          day_subpb_sol_values[j - 1][m];
                if (current > maxi)
                    maxi = current;
            }
            if (maxi == -std::numeric_limits<double>::max() ) throw std::logic_error("dmlke ?");
            max_values[j].push_back(maxi);
        }
    }

    std::vector<unsigned int> occ_lb(0);
    std::vector<unsigned int> occ_ub(0);

    for (unsigned int j = 0; j < NB_DAYS; j++){
        for (unsigned int n = 0; n < possible_Ns[j].size(); n++){
            if (lower_bounds[j][n] <= max_values[j][n]){
                std::cout << "day " << j << ": N <=  " << possible_Ns[j][n] << std::endl;
                occ_ub.push_back(possible_Ns[j][n]);
                break;
            }
        }
        for (int n = possible_Ns[j].size() - 1; n >= 0; n--){
            if (lower_bounds[j][n] <= max_values[j][n]){
                std::cout << "day " << j << ":  " << possible_Ns[j][n] << " <= N" << std::endl;
                occ_lb.push_back(possible_Ns[j][n]);
                break;
            }
        }
    }

    unsigned int c_discontinuous = 0;
    for (unsigned int j = 0; j < NB_DAYS; j++)
        for (unsigned int n = 0; n < possible_Ns[j].size(); n++)
            if (lower_bounds[j][n] <= max_values[j][n])
                c_discontinuous++;

    unsigned int c_continuous = 0;
    for (unsigned int j = 0; j < NB_DAYS; j++)
        c_continuous += occ_ub[j] - occ_lb[j];

    std::cout << c_continuous << " " << c_discontinuous << std::endl;

}
