#include "day_sub_problem.h"

void recursively_find_possible_quantities(const std::vector<unsigned int> &max_per_size,
                                          std::vector<unsigned int> &possible_quantities,
                                          const unsigned int &current_quantity, const unsigned int &current_index) {
    // 2, 3, 4, 5, 6, 7, 8 make up 7 possible family sizes
    if (current_index == 7) {
        possible_quantities.push_back(current_quantity);
        return;
    }
    for (unsigned int i = 0; i <= max_per_size[current_index]; i++){
        recursively_find_possible_quantities(max_per_size, possible_quantities, current_quantity + i * (current_index + 2), current_index + 1);
    }
}

std::vector<unsigned int> get_possible_quantities(const std::vector<unsigned int> &sizes) {
    if (sizes.size() > 12) {
        //throw std::invalid_argument("too many sizes");
        std::vector<unsigned int> max_per_size(7, 0);
        for (unsigned int i = 0; i < sizes.size(); i++)
            max_per_size[sizes[i] - 2]++;
        std::vector<unsigned int> res(0);
        recursively_find_possible_quantities(max_per_size, res, 0, 0);
        return get_unique_vector(res);
    }
    unsigned long n = pow(2, sizes.size());
    std::vector<unsigned int> res (n, 0);
    for (unsigned int i = 0; i < n; i++){
        std::vector<bool> mask = get_binary_rep(i, sizes.size());
        for (unsigned int j = 0; j < sizes.size(); j++)
            if (mask[j])
                res[i] += sizes[j];
    }
    return get_unique_vector(res);
}

void recursively_find_quantity_decompositions(const std::vector<unsigned int> &max_per_size,
                                              const unsigned int& quantity,
                                              std::vector<std::vector<unsigned int>> &decompositions,
                                              std::vector<unsigned int> &current_decomposition,
                                              const unsigned int &current_quantity, const unsigned int &current_index) {
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
                                                 current_quantity + i * (current_index + 2), current_index + 1);
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
                                                      const unsigned int &current_index) {
    count++; if (count > 25000000) return;
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
                                                     nb_left_to_see_per_size, masks, current_mask, count, current_index + 1);
    current_mask[current_index] = true;
    current_decomposition[sizes[current_index] - 2]++;
    recursively_find_all_masks_fitting_decomposition(sizes, decomposition, current_decomposition,
                                                     nb_left_to_see_per_size, masks, current_mask, count, current_index + 1);
    current_mask[current_index] = false;
    current_decomposition[sizes[current_index] - 2]--;
    nb_left_to_see_per_size[sizes[current_index] - 2]++;
}

std::vector<std::vector<bool>> get_masks_fitting_quantity(const std::vector<unsigned int> &sizes, const unsigned int &q) {
    if (sizes.size() > 12) {
        //throw std::invalid_argument("too many sizes");
        std::vector<unsigned int> max_per_size(7, 0);
        for (unsigned int i = 0; i < sizes.size(); i++)
            max_per_size[sizes[i] - 2]++;
        std::vector<std::vector<unsigned int>> possible_decompositions(0);
        std::vector<unsigned int> current_decomposition(7, 0);
        recursively_find_quantity_decompositions(max_per_size, q, possible_decompositions, current_decomposition, 0, 0);
        std::vector<std::vector<bool>> res(0);
        std::vector<bool> current_mask(sizes.size(), false);
        unsigned int count = 0;
        for (std::vector<unsigned int>& decomposition: possible_decompositions)
            recursively_find_all_masks_fitting_decomposition(sizes, decomposition, current_decomposition, max_per_size, res, current_mask, count, 0);
        if (count % 10000 == 1)
            return std::vector<std::vector<bool>>(0);
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
        if (sum == q)
            res.push_back(mask);
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
    // for every level of cost per person, get the possible number of people ( = quantities) that we can assign
    std::vector<std::vector<unsigned int>> unique_possible_quantities;
    for (unsigned int i = 0; i < possible_sizes.size(); i++) {
        unique_possible_quantities.push_back(get_possible_quantities(possible_sizes[i]));
    }

    unsigned int min_allowed = (unsigned int)(std::max(int(presets.get_occupancy_lb(day)) - int(presets.get_presets_occupancy(day)), 0));
    unsigned int max_allowed = presets.get_occupancy_ub(day) - presets.get_presets_occupancy(day);


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
    unsigned int total = 0;
    for (unsigned int i = 0; i < solutions.size(); i++){
        std::vector<unsigned int> count_appearances_per_family (families_to_consider.size(), 0);
        unsigned int c = 1;
        for (unsigned int j = 0; j < possible_costs.size(); j++){
            std::vector<std::vector<bool>> masks = get_masks_fitting_quantity(possible_sizes[j], solutions[i][j]);
            if (masks.size() == 0) {
                std::cout << "We had to stop due to excessive time needed to determine possible assignations or counter-assignations." << std::endl;
                make_incomplete_counter_assignations(families_to_consider, corresponding_families, solutions, day, counter_assignations);
                return true;
            }
            c = c * masks.size();
            for (unsigned int m = 0; m < masks.size(); m++)
                for (unsigned int l = 0; l < masks[m].size(); l++)
                    if (masks[m][l])
                        count_appearances_per_family[corresponding_families[j][l]]++;
            for (unsigned int s = 0; s < corresponding_families[j].size(); s++) {
                if (count_appearances_per_family[corresponding_families[j][s]] != masks.size())
                    always_appears[corresponding_families[j][s]] = false;
                if (count_appearances_per_family[corresponding_families[j][s]] > 0)
                    never_appears[corresponding_families[j][s]] = false;

            }
        }
        total += c;
    }

    std::cout << "Found " << solutions.size() << " equivalent solutions, yielding " << total << " distinct solutions." << std::endl;

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



