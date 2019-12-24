#include "Presets.h"
#include <unordered_set>

Presets::Presets(std::vector<std::vector<unsigned int>> &family_data): family_data(family_data) {
    presets = std::vector<preset> (0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        presets.push_back(std::vector<status>(K_MAX, ALLOWED));
    preset_occupancy = std::vector<unsigned int>(NB_DAYS, 0);
    preset_cardinal = std::vector<unsigned int>(NB_DAYS, 0);
    occupancy_lb = std::vector<unsigned int>(NB_DAYS, MIN_NB_PEOPLE_PER_DAY);
    occupancy_ub = std::vector<unsigned int>(NB_DAYS, MAX_NB_PEOPLE_PER_DAY);
    occupancy_explicit_ub = std::vector<unsigned int>(NB_DAYS, MAX_NB_PEOPLE_PER_DAY);
    day_costs_lb = std::vector<float>(NB_DAYS, 0);
    is_already_assigned = std::vector<bool>(NB_FAMILIES, false);
    sorted_families = std::vector<uint_pair>(0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        sorted_families.push_back(uint_pair(i, family_data[i][NB_CHOICES]));
    }
    sort_by_second(sorted_families);
}

void Presets::assign_family(const unsigned int &i, const unsigned int &k, const bool &compute_stuff) {
    deassign_family(i, false);
    presets[i] = get_assignation_preset(k);
    presets_costs += CONSTANT_COST[k] + family_data[i][NB_CHOICES] * MARGINAL_COST[k];
    preset_occupancy[family_data[i][k]] += family_data[i][NB_CHOICES];
    preset_cardinal[family_data[i][k]] += 1;
    is_already_assigned[i] = true;
    nb_assignments++;
    if (compute_stuff) {
        compute_all_bounds();
        compute_feasibility();
    }
}

void Presets::deassign_family(const unsigned int &i, const bool &compute_bounds) {
    for (unsigned int k = 0; k < K_MAX && is_already_assigned[i]; k++)
        if (presets[i][k] == COMPULSORY) {
            presets_costs -= CONSTANT_COST[k] + family_data[i][NB_CHOICES] * MARGINAL_COST[k];
            preset_occupancy[family_data[i][k]] -= family_data[i][NB_CHOICES];
            preset_cardinal[family_data[i][k]] -= 1;
            is_already_assigned[i] = false;
            nb_assignments--;
        }
    presets[i] = std::vector<status>(K_MAX, ALLOWED);
    if (compute_bounds){
        compute_all_bounds();
        compute_feasibility();
    }
}

void Presets::write_solution(const std::string &filename) const {
    write_solution_(get_solution(), filename);
}

std::vector<unsigned int> Presets::get_solution() const {
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

void Presets::compute_feasibility() {
    is_feasible_ = true;
    for (unsigned int i = 0; i < NB_DAYS && is_feasible_; i++)
        if (preset_occupancy[i] > occupancy_ub[i] || occupancy_ub[i] < occupancy_lb[i])
            is_feasible_ = false;
}

void Presets::compute_occupancy_bounds() {
    unsigned int max_day_cost = 500;
    occupancy_explicit_ub = std::vector<unsigned int> (NB_DAYS, 0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        for (unsigned int k = 0; k < K_MAX; k++)
            if (presets[i][k] != FORBIDDEN)
               occupancy_explicit_ub[family_data[i][k]] += family_data[i][NB_CHOICES];

    for (unsigned int i = 0; i < NB_DAYS; i++){
        occupancy_ub[i] = std::min(MAX_NB_PEOPLE_PER_DAY, occupancy_explicit_ub[i]);
        occupancy_lb[i] = std::max(MIN_NB_PEOPLE_PER_DAY, preset_occupancy[i]);
    }

    for (int i = NB_DAYS - 1; i > 0; i--){
        if (i < NB_DAYS - 1 && occupancy_ub[i] > occupancy_ub[i + 1] &&
            (occupancy_ub[i] - MIN_NB_PEOPLE_PER_DAY)/400.*pow(occupancy_ub[i], 0.5+(occupancy_ub[i] - occupancy_ub[i + 1])/50.) > 500){
            unsigned int a;
            for (unsigned int j = 0; j < 10; j++) {
                a = occupancy_ub[i];
                occupancy_ub[i] = ceil(occupancy_ub[i + 1] + 50*(log(400.*max_day_cost/float(a - MIN_NB_PEOPLE_PER_DAY))/log(a) - 0.5));
            }
        }
    }

    for (unsigned int i = 0; i < NB_DAYS; i++){
        if (i < NB_DAYS - 1 && occupancy_lb[i] > MIN_NB_PEOPLE_PER_DAY && occupancy_ub[i] <= occupancy_lb[i + 1]){// &&
            //(occupancy_lb[i] - MIN_NB_PEOPLE_PER_DAY)/400.*pow(occupancy_lb[i], 0.5+(occupancy_lb[i + 1] - occupancy_lb[i])/50.) > 500){
            //std::cout << i << " " << occupancy_lb[i] << " " << occupancy_lb[i + 1] << " lower bound became ";
            unsigned int a;
            for (unsigned int j = 0; j < 10; j++) {
                a = occupancy_lb[i];
                occupancy_lb[i] = std::max(a, (unsigned int)(std::max(0., ceil(occupancy_lb[i + 1] - 50.*(log(400.*max_day_cost/float(a - MIN_NB_PEOPLE_PER_DAY))/log(a) - 0.5)))));
            }
            //std::cout << occupancy_lb[i] << std::endl;
            a = occupancy_lb[i];
            occupancy_ub[i + 1] = std::min(occupancy_ub[i + 1], (unsigned int)(floor(occupancy_ub[i] + 50.*(log(400.*max_day_cost/float(a - MIN_NB_PEOPLE_PER_DAY))/log(a) - 0.5))));
            //std::cout << preset_lower_bound[i] << " " << occupancy_lb[i] << " " << occupancy_ub[i] << " " << occupancy_lb[i + 1] << " " << occupancy_ub[i + 1] << " " << preset_upper_bound[i + 1] << std::endl;
        }
    }

    std::cout << "Upper bounds:" << std::endl;
    print_nicely(occupancy_ub, 5);
}

float Presets::get_day_cost_lb_(const unsigned int &i) {
    unsigned int gap;
    if (i == NB_DAYS - 1)
        gap = 0;
    else if (occupancy_ub[i] <= occupancy_lb[i + 1])
        gap = occupancy_lb[i + 1] - occupancy_ub[i];
    else if (occupancy_ub[i + 1] <= occupancy_lb[i])
        gap = occupancy_lb[i] - occupancy_ub[i + 1];
    else
        gap = 0;
    return (occupancy_lb[i] - MIN_NB_PEOPLE_PER_DAY)/400.*pow(occupancy_lb[i], 0.5 + gap/50.);
}

unsigned int Presets::get_additional_day_cost_lb(const unsigned int &i) {
    // i is a family index. we compute a lower bound on the day cost of assigning the family.
    if (is_already_assigned[i]) return 0;
    std::unordered_set<unsigned int> impacted_days(0);
    for (unsigned int k = 0; k < K_MAX; k++){
        impacted_days.emplace(family_data[i][k]);
        if (int(family_data[i][k]) - 1 >= 0)
            impacted_days.emplace(family_data[i][k] - 1);
    }
    unsigned int min_possible_increase = 500;
    for (unsigned int k_assign = 0; k_assign < K_MAX; k_assign++){
        if (preset_occupancy[family_data[i][k_assign]] + family_data[i][NB_CHOICES] > occupancy_ub[family_data[i][k_assign]])
            continue;
        // simulate assigning to k
        std::vector<unsigned int> occupancy_ub_save(K_MAX);
        std::vector<unsigned int> occupancy_lb_save(K_MAX);
        for (unsigned int k = 0; k < K_MAX; k++){
            unsigned int &current_day = family_data[i][k];
            occupancy_lb_save[k] = occupancy_lb[current_day];
            occupancy_ub_save[k] = occupancy_ub[current_day];
            if (k == k_assign)
                occupancy_lb[current_day] = std::max(occupancy_lb[current_day], preset_occupancy[current_day] + family_data[i][NB_CHOICES]);
            else
                occupancy_ub[current_day] = std::min(occupancy_ub[current_day], occupancy_explicit_ub[current_day] - family_data[i][NB_CHOICES]);
        }
        // compute the increase
        int possible_increase = 0;
        for (const unsigned int& current_day: impacted_days)
            possible_increase += get_day_cost_lb_(current_day) - day_costs_lb[current_day];
        if (possible_increase < min_possible_increase)
            min_possible_increase = possible_increase;
        // restore the previous bounds
        for (int k = K_MAX - 1; k >= 0; k--){
            unsigned int &current_day = family_data[i][k];
            occupancy_lb[current_day] = occupancy_lb_save[k];
            occupancy_ub[current_day] = occupancy_ub_save[k];
        }
    }
    return min_possible_increase;
}

void Presets::compute_all_bounds() {
    clock_t t0 = clock();
    compute_occupancy_bounds();
    std::cout << "time spent on occupancy_bounds: " << clock() - t0 << std::endl;

    t0 = clock();
    std::cout << "Nb of families preset and number of people" << std::endl;
    print_nicely(preset_occupancy, 5);
    print_nicely(preset_cardinal, 5);
    std::cout << "time spent on printing: " << clock() - t0 << std::endl;

    t0 = clock();
    for (unsigned int i = 0; i < NB_DAYS; i++)
        day_costs_lb[i] = get_day_cost_lb_(i);
    compute_bottleneck_bounds();
    std::cout << "time spent on bottleneck bounds: " << clock() - t0 << std::endl;
    t0 = clock();
    compute_k_fold_bottleneck_ub();
    std::cout << "time spent on k-fold bottleneck bounds: " << clock() - t0 << std::endl;
    t0 = clock();
    day_cost_lower_bound = 0;
    for (unsigned int i = 0; i < NB_DAYS; i++)
        day_cost_lower_bound += get_day_cost_lb_(i);
    std::cout << "time spent on computing the day costs lower bound: " << clock() - t0 << std::endl;
    t0 = clock();
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        day_cost_lower_bound += get_additional_day_cost_lb(i);
    std::cout << "time spent on computing the per-family additional day costs: " << clock() - t0 << std::endl;
    std::cout << "cost from days is at least " << day_cost_lower_bound;
    if (day_cost_lower_bound < 4500)
        std::cout << " before being set to 4500";
    std::cout << std::endl;
}

void Presets::compute_bottleneck_bounds() {
    // this gives the bottleneck bound. To know how many families can be accommodated during the day, you need to add the
    // nb of families already put on that day by the preset.
    // first is the min number of *additional* families that can be assigned to the day, second is the max number.
    unsigned int current_family, current_day;
    bottleneck_bounds = std::vector<uint_pair> (NB_DAYS, uint_pair(0, 0));
    std::vector<unsigned int> possible_family_cum_sizes(NB_DAYS, 0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        for (unsigned int k = 0; k < K_MAX; k++) {
            current_family = sorted_families[i].first;
            current_day = family_data[current_family][k];
            if (presets[current_family][k] == ALLOWED &&
                possible_family_cum_sizes[current_day] + family_data[current_family][NB_CHOICES] <= occupancy_ub[current_day] - preset_occupancy[current_day]){
                possible_family_cum_sizes[current_day] += family_data[current_family][NB_CHOICES];
                bottleneck_bounds[current_day].second++;
            }
        }
    }
    possible_family_cum_sizes = std::vector<unsigned int>(NB_DAYS, 0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        for (unsigned int k = 0; k < K_MAX; k++) {
            current_family = sorted_families[NB_FAMILIES - i - 1].first;
            current_day = family_data[current_family][k];
            if (presets[current_family][k] == ALLOWED &&
                possible_family_cum_sizes[current_day] < occupancy_lb[current_day] - preset_occupancy[current_day]){
                possible_family_cum_sizes[current_day] += family_data[current_family][NB_CHOICES];
                bottleneck_bounds[current_day].first++;
            }
        }
    }
    std::cout << "bottleneck bounds (min and then max):" << std::endl;
    for (unsigned int i = 0; i < NB_DAYS; i++) {
        unsigned int n = nb_chiffres(bottleneck_bounds[i].first);
        if (n >= 4)
            throw std::logic_error("Why is this integer so big ?");
        std::cout << bottleneck_bounds[i].first << std::string(5 - n, ' ');
    }
    std::cout << std::endl;
    for (unsigned int i = 0; i < NB_DAYS; i++) {
        unsigned int n = nb_chiffres(bottleneck_bounds[i].second);
        if (n >= 4)
            throw std::logic_error("Why is this integer so big ?");
        std::cout << bottleneck_bounds[i].second << std::string(5 - n, ' ');
    }
    std::cout << std::endl;
}

void Presets::compute_k_fold_bottleneck_ub() {
    unsigned int current_fam, current_day;
    k_fold_bottleneck_bounds = std::vector<std::vector<unsigned int>>(NB_DAYS, std::vector<unsigned int>(K_MAX, 0));
    std::vector<std::vector<unsigned int>> pos_fam_cum_sizes(NB_DAYS, std::vector<unsigned int>(K_MAX, 0));
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        for (unsigned int k = 0; k < K_MAX; k++){
            current_fam = sorted_families[i].first;
            current_day = family_data[current_fam][k];
            if (presets[current_fam][k] == ALLOWED &&
                pos_fam_cum_sizes[current_day][k] + family_data[current_fam][NB_CHOICES] <= occupancy_ub[current_day] - preset_occupancy[current_day]){
                pos_fam_cum_sizes[current_day][k] += family_data[current_fam][NB_CHOICES];
                k_fold_bottleneck_bounds[current_day][k]++;
            }
        }
    }
}

std::vector<unsigned int> Presets::get_largest_unassigned_families() const {
    std::vector<unsigned int> largest(0);
    unsigned int maxi = 1;
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        if (!is_already_assigned[i] && family_data[i][NB_CHOICES] > maxi){
            maxi = family_data[i][NB_CHOICES];
            largest.clear();
            largest.push_back(i);
        }
        if (!is_already_assigned[i] && family_data[i][NB_CHOICES] == maxi)
            largest.push_back(i);
    }
    return largest;
}

int Presets::get_largest_unassigned_strategic_family() const {
    std::vector<unsigned int> largest = get_largest_unassigned_families();
    if (largest.size() == 0)
        return -1;
    std::vector<unsigned int> best_to_second_differences(largest.size(), 0);
    for (unsigned int m = 0; m < largest.size(); m++) {
        unsigned int i = largest[m];
        unsigned int best_k = K_MAX, second_best_k = K_MAX;
        for (unsigned int k = 0; k < K_MAX; k++){
            if (preset_occupancy[family_data[i][k]] + family_data[i][NB_CHOICES] > occupancy_ub[family_data[i][k]])
                continue;
            if (k > best_k && k < second_best_k)
                second_best_k = k;
            if (k < best_k)
                best_k = k;
        }
        best_to_second_differences[m] = second_best_k - best_k;
    }
    unsigned int mini = K_MAX;
    unsigned int i_mini = 0;
    for (unsigned int m = 0; m < largest.size(); m++) {
        if (best_to_second_differences[m] < mini){
            mini = best_to_second_differences[m];
            i_mini = largest[m];
        }
    }
    return i_mini;
}
