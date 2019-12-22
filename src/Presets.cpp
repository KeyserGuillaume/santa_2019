#include "Presets.h"

Presets::Presets(std::vector<std::vector<unsigned int>> &family_data): family_data(family_data) {
    presets = std::vector<preset> (0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        presets.push_back(std::vector<status>(K_MAX, ALLOWED));
    preset_occupancy = std::vector<unsigned int>(NB_DAYS, 0);
    preset_cardinal = std::vector<unsigned int>(NB_DAYS, 0);
    occupancy_lower_bounds = std::vector<unsigned int>(NB_DAYS, MIN_NB_PEOPLE_PER_DAY);
    occupancy_upper_bounds = std::vector<unsigned int>(NB_DAYS, MAX_NB_PEOPLE_PER_DAY);
    day_costs_lower_bounds = std::vector<float>(NB_DAYS, 0);
    is_already_assigned = std::vector<bool>(NB_FAMILIES, false);
    sorted_families = std::vector<uint_pair>(0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        sorted_families.push_back(uint_pair(i, family_data[i][NB_CHOICES]));
    }
    sort_by_second(sorted_families);
    compute_all_bounds();
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
        if (preset_occupancy[i] > occupancy_upper_bounds[i] || occupancy_upper_bounds[i] < occupancy_lower_bounds[i])
            is_feasible_ = false;
}

void Presets::compute_occupancy_bounds() {
    std::vector<unsigned int> preset_lower_bound(NB_DAYS, 0);
    std::vector<unsigned int> preset_upper_bound(NB_DAYS, 0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        for (unsigned int k = 0; k < K_MAX; k++) {
            if (presets[i][k] == COMPULSORY) {
                preset_lower_bound[family_data[i][k]] += family_data[i][NB_CHOICES];
            }
            if (presets[i][k] != FORBIDDEN) {
                preset_upper_bound[family_data[i][k]] += family_data[i][NB_CHOICES];
            }
        }
    }

    for (unsigned int i = 0; i < NB_DAYS; i++){
        occupancy_upper_bounds[i] = std::min(MAX_NB_PEOPLE_PER_DAY, preset_upper_bound[i]);
        occupancy_lower_bounds[i] = std::max(MIN_NB_PEOPLE_PER_DAY, preset_lower_bound[i]);
    }

    std::cout << "Upper bounds:" << std::endl;
    print_nicely(occupancy_upper_bounds);

    for (int i = NB_DAYS - 1; i > 0; i--){
        if (i < NB_DAYS - 1 && occupancy_upper_bounds[i] > occupancy_upper_bounds[i + 1] &&
            (occupancy_upper_bounds[i] - MIN_NB_PEOPLE_PER_DAY)/400.*pow(occupancy_upper_bounds[i], 0.5+(occupancy_upper_bounds[i] - occupancy_upper_bounds[i + 1])/50.) > 500){
            unsigned int a;
            for (unsigned int j = 0; j < 10; j++) {
                a = occupancy_upper_bounds[i];
                occupancy_upper_bounds[i] = ceil(occupancy_upper_bounds[i + 1] + 50*(log(400.*500./float(a - MIN_NB_PEOPLE_PER_DAY))/log(a) - 0.5));
            }
        }
    }

    std::cout << "Upper bounds:" << std::endl;
    print_nicely(occupancy_upper_bounds);

    for (unsigned int i = 0; i < NB_DAYS; i++){
        if (i < NB_DAYS - 1 && occupancy_lower_bounds[i] > MIN_NB_PEOPLE_PER_DAY && occupancy_upper_bounds[i] <= occupancy_lower_bounds[i + 1]){// &&
            //(occupancy_lower_bounds[i] - MIN_NB_PEOPLE_PER_DAY)/400.*pow(occupancy_lower_bounds[i], 0.5+(occupancy_lower_bounds[i + 1] - occupancy_lower_bounds[i])/50.) > 500){
            //std::cout << i << " " << occupancy_lower_bounds[i] << " " << occupancy_lower_bounds[i + 1] << " lower bound became ";
            unsigned int a;
            for (unsigned int j = 0; j < 10; j++) {
                a = occupancy_lower_bounds[i];
                occupancy_lower_bounds[i] = std::max(a, (unsigned int)(std::max(0., ceil(occupancy_lower_bounds[i + 1] - 50.*(log(400.*500./float(a - MIN_NB_PEOPLE_PER_DAY))/log(a) - 0.5)))));
            }
            //std::cout << occupancy_lower_bounds[i] << std::endl;
            a = occupancy_lower_bounds[i];
            occupancy_upper_bounds[i + 1] = std::min(occupancy_upper_bounds[i + 1], (unsigned int)(floor(occupancy_upper_bounds[i] + 50.*(log(400.*500./float(a - MIN_NB_PEOPLE_PER_DAY))/log(a) - 0.5))));
            //std::cout << preset_lower_bound[i] << " " << occupancy_lower_bounds[i] << " " << occupancy_upper_bounds[i] << " " << occupancy_lower_bounds[i + 1] << " " << occupancy_upper_bounds[i + 1] << " " << preset_upper_bound[i + 1] << std::endl;
        }
    }

    std::cout << "Upper bounds:" << std::endl;
    print_nicely(occupancy_upper_bounds);
}

void Presets::compute_day_costs_lb(const unsigned int &i) {
    unsigned int gap;
    if (i == NB_DAYS - 1)
        gap = 0;
    else if (occupancy_upper_bounds[i] <= occupancy_lower_bounds[i + 1])
        gap = occupancy_lower_bounds[i + 1] - occupancy_upper_bounds[i];
    else if (occupancy_upper_bounds[i + 1] <= occupancy_lower_bounds[i])
        gap = occupancy_lower_bounds[i] - occupancy_upper_bounds[i + 1];
    else
        gap = 0;
    day_costs_lower_bounds[i] = (occupancy_lower_bounds[i] - MIN_NB_PEOPLE_PER_DAY)/400.*pow(occupancy_lower_bounds[i], 0.5 + gap/50.);
}

void Presets::compute_all_bounds() {
    compute_occupancy_bounds();

    std::cout << "Nb of families preset and number of people" << std::endl;
    print_nicely(preset_occupancy);
    print_nicely(preset_cardinal);

    for (unsigned int i = 0; i < NB_DAYS; i++)
        compute_day_costs_lb(i);
    compute_bottleneck_bounds();
    compute_k_fold_bottleneck_ub();
    day_cost_lower_bound = 0;
    for (unsigned int i = 0; i < NB_DAYS; i++)
        day_cost_lower_bound += get_day_cost_lower_bound(i);
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
                possible_family_cum_sizes[current_day] + family_data[current_family][NB_CHOICES] <= occupancy_upper_bounds[current_day] - preset_occupancy[current_day]){
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
                possible_family_cum_sizes[current_day] < occupancy_lower_bounds[current_day] - preset_occupancy[current_day]){
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
        std::cout << bottleneck_bounds[i].first << std::string(4 - n, ' ');
    }
    std::cout << std::endl;
    for (unsigned int i = 0; i < NB_DAYS; i++) {
        unsigned int n = nb_chiffres(bottleneck_bounds[i].second);
        if (n >= 4)
            throw std::logic_error("Why is this integer so big ?");
        std::cout << bottleneck_bounds[i].second << std::string(4 - n, ' ');
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
                pos_fam_cum_sizes[current_day][k] + family_data[current_fam][NB_CHOICES] <= occupancy_upper_bounds[current_day] - preset_occupancy[current_day]){
                pos_fam_cum_sizes[current_day][k] += family_data[current_fam][NB_CHOICES];
                k_fold_bottleneck_bounds[current_day][k]++;
            }
        }
    }
}