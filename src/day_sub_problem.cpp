#include "day_sub_problem.h"

void lagrangian_stats(const Presets &presets, const std::vector<float> &lambda,
                      const std::vector<unsigned int> &family_uses_nb);

void assignments_costs_only_lb(const Presets &presets, const std::vector<float> &lambda,
                               const std::vector<std::vector<unsigned int>> &solutions, int &sum_cost_lbs,
                               std::vector<unsigned int> &family_uses_nb);

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
    if (sizes.size() > 3) {
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
        my_set.emplace(res[i]);
    res = std::vector<unsigned int>(my_set.begin(), my_set.end());
    return res;
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


bool anticipate_on_day(const Presets &presets, const unsigned int &day, const unsigned int &reference_cost,
                       std::vector<unsigned int> &assignations, std::vector<unsigned int> &counter_assignations) {
    std::vector<std::vector<unsigned int>> solutions;
    std::vector<std::pair<unsigned int, float>> families_to_consider(0);
    // sort the families by increasing cost per person
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        for (unsigned int k = 0; k < K_MAX; k++)
            if (presets.get_family_data(i, k) == day && presets[i][k] == ALLOWED)
                families_to_consider.push_back(std::pair<unsigned int, float>(i, float(CONSTANT_COST[k])/float(presets.get_family_size(i)) + MARGINAL_COST[k]));
    sort_by_second(families_to_consider);

    if (families_to_consider.size() == 0)
        return true;

    // group the families by the ~20 levels of cost per person
    std::vector<float> possible_costs(1, families_to_consider[0].second);
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
    float best_solution = 2000;
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
//    unsigned int total = 0;
//    for (unsigned int i = 0; i < solutions.size(); i++){
//        std::vector<unsigned int> count_appearances_per_family (families_to_consider.size(), 0);
//        unsigned int c = 1;
//        for (unsigned int j = 0; j < possible_costs.size(); j++){
//            std::vector<std::vector<bool>> masks = get_masks_fitting_quantity(possible_sizes[j], solutions[i][j]);
//            if (masks.size() == 0) {
//                std::cout << "We had to stop due to excessive time needed to determine possible assignations or counter-assignations." << std::endl;
//                make_incomplete_counter_assignations(families_to_consider, corresponding_families, solutions, day, counter_assignations);
//                return true;
//            }
//            c = c * masks.size();
//            for (unsigned int m = 0; m < masks.size(); m++)
//                for (unsigned int l = 0; l < masks[m].size(); l++)
//                    if (masks[m][l])
//                        count_appearances_per_family[corresponding_families[j][l]]++;
//            for (unsigned int s = 0; s < corresponding_families[j].size(); s++) {
//                if (count_appearances_per_family[corresponding_families[j][s]] != masks.size())
//                    always_appears[corresponding_families[j][s]] = false;
//                if (count_appearances_per_family[corresponding_families[j][s]] > 0)
//                    never_appears[corresponding_families[j][s]] = false;
//
//            }
//        }
//        total += c;
//    }
//
//    std::cout << "Found " << solutions.size() << " equivalent solutions, yielding " << total << " distinct solutions." << std::endl;

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
                                   const std::vector<float> &possible_costs,
                                   const unsigned int &cost_threshold,
                                   const unsigned int &min_occupancy,
                                   const unsigned int &max_occupancy,
                                   std::vector<std::vector<unsigned int>> &solutions,
                                   std::vector<unsigned int> &current_solution,
                                   const unsigned int &nb_assigned_people,
                                   const float &current_cost,
                                   const unsigned int &next_cost_index) {
    if (solutions.size() > 10000) return;
    float nb_needed_people = (nb_assigned_people > min_occupancy) ? 0 : min_occupancy - nb_assigned_people;
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
                        const std::vector<float> &possible_costs,
                        const unsigned int &min_occupancy,
                        const unsigned int &max_occupancy,
                        float &best_solution,
                        const unsigned int &nb_assigned_people,
                        const float &current_cost,
                        const unsigned int &next_cost_index) {
    float nb_needed_people = (nb_assigned_people > min_occupancy) ? 0 : min_occupancy - nb_assigned_people;
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

void make_incomplete_counter_assignations(std::vector<std::pair<unsigned int, float>> families_to_consider,
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
                               const std::vector<float> &possible_costs, const unsigned int &min_occupancy,
                               const unsigned int &max_occupancy, float &best_solution_cost,
                               std::vector<unsigned int> &best_solution, std::vector<unsigned int> &current_solution,
                               const unsigned int &nb_assigned_people, const float &current_cost,
                               const unsigned int &next_cost_index) {
    //unsigned int nb_needed_people = (nb_assigned_people > min_occupancy) ? 0 : min_occupancy - nb_assigned_people;
    unsigned int N = nb_assigned_people;
    float lb = current_cost;
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
//    if (current_cost + nb_needed_people * possible_costs[next_cost_index] > best_solution_cost)
    if (lb > best_solution_cost || N < min_occupancy || N > max_occupancy)
        return;
    if (nb_assigned_people > max_occupancy)
        return;
    if (next_cost_index == possible_costs.size() || nb_assigned_people == max_occupancy){
        if (nb_assigned_people < min_occupancy)
            return;
        //std::cout << current_cost << std::endl;
        if (current_cost < best_solution_cost) {
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

int get_day_lower_bound(const Presets &presets, const unsigned int &day, const unsigned int& min_allowed, const unsigned int& max_allowed, const std::vector<float> &lambda, std::vector<unsigned int> &solution, const bool& print) {
    clock_t t0 = clock();
    std::vector<std::pair<unsigned int, float>> families_to_consider(0);

    // if the passed solution already contains some stuff, use it to get a lower bound on the sub-problem optimal solution
    float best_solution_cost = 10000;
    if (solution.size() > 0){
        best_solution_cost = 0;
        unsigned int N = 0;
        for (unsigned int s = 0; s < solution.size(); s++){
            unsigned int i = solution[s];
            unsigned int k = 0;
            for (; k < K_MAX && presets.get_family_data(i, k) != day; k++) {}
            if (k == K_MAX) throw std::logic_error("What's going on here ?");
            best_solution_cost += float(presets.get_family_size(i)) * ((float(CONSTANT_COST[k]) + lambda[i]) / float(presets.get_family_size(i)) + MARGINAL_COST[k]);
            N += presets.get_family_size(i);
        }
        if (N > max_allowed || N < min_allowed) throw std::logic_error("What's going on here ?");
    }

    // sort the families by increasing cost per person
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        for (unsigned int k = 0; k < K_MAX; k++)
            if (presets.get_family_data(i, k) == day && presets[i][k] == ALLOWED) {
                float cost = (float(CONSTANT_COST[k]) + lambda[i]) / float(presets.get_family_size(i)) + MARGINAL_COST[k];
                families_to_consider.push_back(std::pair<unsigned int, float>(i, cost));
            }

    if (print) std::cout << "time spent on identifying the appropriate families: " << clock() - t0 << std::endl;
    t0 = clock();
    sort_by_second(families_to_consider);
    if (print) std::cout << "time spent on sorting the appropriate families: " << clock() - t0 << std::endl;
    t0 = clock();

    if (families_to_consider.size() == 0)
        return presets.get_presets_costs(day);

    // group the families by the ~20 levels of cost per person
    std::vector<float> possible_costs(1, families_to_consider[0].second);
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
    std::vector<unsigned int> best_solution(possible_costs.size(), 0);
    std::vector<unsigned int> current_solution(possible_costs.size());
    find_one_optimal_solution(unique_possible_quantities, possible_costs, min_allowed, max_allowed,
                              best_solution_cost,
                              best_solution, current_solution, 0, 0, 0);

//    // check with the best solution computed less quickly
//    if (day == 80){
//        float best_solution_cost2 = 2000;
//        find_best_solution(unique_possible_quantities, possible_costs, min_allowed, max_allowed, best_solution_cost2, 0, 0, 0);
//        if (best_solution_cost != best_solution_cost2) throw std::logic_error("lkbhk");
//    }

    if (best_solution_cost == 10000)
        throw std::logic_error("mdk");

    if (print) std::cout << "time spent on getting an optimal solution: " << clock() - t0 << std::endl;
    t0 = clock();

    if (print) std::cout << "The day " << day << " has lower bound " << best_solution_cost << std::endl;

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

void assignments_costs_only_lb(const Presets &presets, const std::vector<float> &lambda, double &sum_cost_lbs,
                               std::vector<unsigned int> &family_uses_nb) {
    sum_cost_lbs = 0;
    family_uses_nb = std::vector<unsigned int>(NB_FAMILIES, 0);
    std::vector<std::vector<unsigned int>> solutions(NB_DAYS, std::vector<unsigned int> (0));
    for (unsigned int j = 0; j < NB_DAYS; j++) {
        unsigned int min_allowed = (unsigned int)(std::max(int(presets.get_occupancy_lb(j)) - int(presets.get_presets_occupancy(j)), 0));
        unsigned int max_allowed = presets.get_occupancy_ub(j) - presets.get_presets_occupancy(j);
        sum_cost_lbs += get_day_lower_bound(presets, j, min_allowed, max_allowed, lambda, solutions[j], false);
        for (unsigned int l = 0; l < solutions[j].size(); l++)
            family_uses_nb[solutions[j][l]]++;
    }
}

void joint_costs_DP_lb(const Presets &presets, const std::vector<float> &lambda, double &costs_lb,
                       std::vector<unsigned int> &family_uses_nb) {
    family_uses_nb = std::vector<unsigned int>(NB_FAMILIES, 0);
    std::vector<std::vector<std::vector<unsigned int>>> solutions(NB_DAYS, std::vector<std::vector<unsigned int>>(0));

    // solve all the day subproblems
    std::vector<std::vector<unsigned int>> possible_quantities = presets.get_possible_quantities_per_day();
    std::vector<std::vector<double>> day_subpb_sol_value(NB_DAYS, std::vector<double>(0));
    for (unsigned int j = 0; j < NB_DAYS; j++) {
        solutions[j] = std::vector<std::vector<unsigned int>>(possible_quantities[j].size(),
                                                              std::vector<unsigned int>(0));
        day_subpb_sol_value[j] = std::vector<double>(possible_quantities[j].size(), 0);
        for (unsigned int n = 0; n < possible_quantities[j].size(); n++) {
            unsigned int N = possible_quantities[j][n] - presets.get_presets_occupancy(j);
            day_subpb_sol_value[j][n] = get_day_lower_bound(
                    presets, j, N, N, lambda, solutions[j][n], false);
        }
    }

    // dynamic programming
    std::vector<std::vector<double>> value_functions(NB_DAYS, std::vector<double>(0));
    std::vector<std::vector<unsigned int>> argmin_value_functions(NB_DAYS, std::vector<unsigned int>(0));
    std::vector<unsigned int> &last_poss_qties = possible_quantities[NB_DAYS - 1];
    for (unsigned int m = 0; m < last_poss_qties.size(); m++) {
        value_functions[NB_DAYS - 1].push_back(
                get_day_cost(last_poss_qties[m], last_poss_qties[m]) + day_subpb_sol_value[NB_DAYS - 1][m]);
    }
    for (int i = NB_DAYS - 2; i >= 0; i--) {
        for (unsigned int n = 0; n < possible_quantities[i].size(); n++) {
            double mini = 10000000, current;
            unsigned int argmin;
            for (unsigned int m = 0; m < possible_quantities[i + 1].size(); m++) {
                current = value_functions[i + 1][m] +
                          get_day_cost(possible_quantities[i][n], possible_quantities[i + 1][m]) +
                          day_subpb_sol_value[i][n];
                if (current < mini) {
                    mini = current;
                    argmin = m;
                }
            }
            if (mini == 10000000) throw std::logic_error("dmlke ?");
            value_functions[i].push_back(mini);
            argmin_value_functions[i + 1].push_back(argmin);
        }
    }

    costs_lb = 100000000;
    unsigned int argmin;
    for (unsigned int m = 0; m < possible_quantities[0].size(); m++) {
        if (value_functions[0][m] < costs_lb) {
            costs_lb = value_functions[0][m];
            argmin = m;
        }
    }
    if (costs_lb == 10000000) throw std::logic_error("dmlke ?");

    std::vector<unsigned int> optimal_Ns(NB_DAYS, 0);
    std::vector<unsigned int> optimal_Ns_indexes(NB_DAYS, 0);
    optimal_Ns[0] = possible_quantities[0][argmin];
    for (unsigned int i = 1; i < NB_DAYS; i++) {
        argmin = argmin_value_functions[i][argmin];
        optimal_Ns[i] = possible_quantities[i][argmin];
        optimal_Ns_indexes[i] = argmin;
    }
//    print_nicely(optimal_Ns);

    for (unsigned int j = 0; j < NB_DAYS; j++){
        std::vector<unsigned int> &day_optimal_solution = solutions[j][optimal_Ns_indexes[j]];
        for (unsigned int l = 0; l < day_optimal_solution.size(); l++) {
            family_uses_nb[day_optimal_solution[l]]++;
        }
    }

    bool is_primal_feasible = true;
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        if (!presets.is_family_alr_assigned(i) && family_uses_nb[i] != 1)
            is_primal_feasible = false;
    }
    if (is_primal_feasible && costs_lb < BEST_SOLUTION) {
        presets.write_solution("lagrangian_bound_solution");
    }
}

unsigned int get_lagrangian_lb(const Presets &presets, std::vector<float> &lambda, bool &is_primal_feasible,
                               const bool &stop_if_decrease) {
    double costs_lb = 0;
    is_primal_feasible = false;
    float sum_lambda = 0;
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        sum_lambda += lambda[i];
    float res = 0;
    std::vector<unsigned int> family_uses_nb;
    std::vector<unsigned int> history(0);
    std::vector<float> old_lambda;
    for (unsigned int t = 0; t < 1000; t++) {
    //for (unsigned int t = 0; t < 10 || t < min_iter || history[t - 1] > history[t - 10]; t++) {
        float alpha = 1;//floor(30./float(sqrt(t) + 1)) + 1;

        //assignments_costs_only_lb(presets, lambda, costs_lb, family_uses_nb);
        joint_costs_DP_lb(presets, lambda, costs_lb, family_uses_nb);
        double new_res = costs_lb - sum_lambda;
        std::cout << new_res << std::endl;
        if (new_res < res && stop_if_decrease) {
            lambda = old_lambda;
            return new_res;
        }
        else
            res = new_res;
        if (res > BEST_SOLUTION)
            return new_res;
        history.push_back(new_res);
        bool found = false;
        old_lambda = lambda;
        for (unsigned int i = 0; i < NB_FAMILIES; i++) {
            if (presets.is_family_alr_assigned(i)) continue;
            float delta;
            if (family_uses_nb[i] > 1) {
                delta = alpha * (family_uses_nb[i] - 1);
            } else if (family_uses_nb[i] < 1) {
                delta = -alpha;
            }
            if (family_uses_nb[i] != 1) {
                found = true;
                if (i != 0) {
                    lambda[i] += delta;
                    sum_lambda += delta;
                } else {
                    for (unsigned int j = 1; j < NB_FAMILIES; j++) {
                        lambda[j] -= delta;
                    }
                    sum_lambda -= (NB_FAMILIES - 1) * delta;
                }
            }
        }
        if (!found) {
            is_primal_feasible = true;
            return res;
        }

        if (t % 100 == 0)
            lagrangian_stats(presets, lambda, family_uses_nb);
    }
    return res;
}

void lagrangian_stats(const Presets &presets, const std::vector<float> &lambda,
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
    std::cout << "stats on families never used : what lambdas ?" << std::endl;
    print_nicely(histo_family_never_used_lambdas, 5);
    std::cout << "stats on families used twice : what sizes ?" << std::endl;
    print_nicely(std::vector<unsigned int>({2, 3, 4, 5, 6, 7, 8}), 4);
    print_nicely(histo_family_used_twice_sizes, 4);
    std::cout << "stats on families never twice : what lambdas ?" << std::endl;
    print_nicely(histo_family_used_twice_lambdas, 5);
}


std::vector<float> try_lagrangian_thing(const Presets &presets) {
    std::vector<float> lambda(NB_FAMILIES, 0);
    bool b;
    std::cout << "ended with " << get_lagrangian_lb(presets, lambda, b, false) << std::endl;
    return lambda;
}


