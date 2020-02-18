#include "Presets.h"
#include <functional>

Presets::Presets(std::vector<std::vector<unsigned int>> *family_data): family_data(family_data) {
    presets = std::vector<preset> (0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        presets.push_back(std::vector<status>(K_MAX, ALLOWED));
    preset_occupancy = std::vector<unsigned int>(NB_DAYS, 0);
    preset_cardinal = std::vector<unsigned int>(NB_DAYS, 0);
    occupancy_lb = std::vector<unsigned int>(NB_DAYS, MIN_NB_PEOPLE_PER_DAY);
    occupancy_ub = std::vector<unsigned int>(NB_DAYS, MAX_NB_PEOPLE_PER_DAY);
    occupancy_explicit_ub = std::vector<unsigned int>(NB_DAYS, MAX_NB_PEOPLE_PER_DAY);
    day_costs_lb = std::vector<double>(NB_DAYS, 0);
    is_already_assigned = std::vector<bool>(NB_FAMILIES, false);
    sorted_families = std::vector<uint_pair>(0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        sorted_families.push_back(uint_pair(i, (*family_data)[i][NB_CHOICES]));
    }
    sort_by_second(sorted_families);
    prescribed_occupancy_ub = std::vector<uint_pair>(0);
    prescribed_occupancy_lb = std::vector<uint_pair>(0);
}

Presets::Presets(std::vector<std::vector<unsigned int>> *family_data, const std::string &filename): family_data(family_data) {
    // copied from method just above, why not do some refacto
    preset_occupancy = std::vector<unsigned int>(NB_DAYS, 0);
    preset_cardinal = std::vector<unsigned int>(NB_DAYS, 0);
    occupancy_lb = std::vector<unsigned int>(NB_DAYS, MIN_NB_PEOPLE_PER_DAY);
    occupancy_ub = std::vector<unsigned int>(NB_DAYS, MAX_NB_PEOPLE_PER_DAY);
    occupancy_explicit_ub = std::vector<unsigned int>(NB_DAYS, MAX_NB_PEOPLE_PER_DAY);
    day_costs_lb = std::vector<double>(NB_DAYS, 0);
    is_already_assigned = std::vector<bool>(NB_FAMILIES, false);
    sorted_families = std::vector<uint_pair>(0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        sorted_families.push_back(uint_pair(i, (*family_data)[i][NB_CHOICES]));
    }
    sort_by_second(sorted_families);
    prescribed_occupancy_ub = std::vector<uint_pair>(0);
    prescribed_occupancy_lb = std::vector<uint_pair>(0);

    presets = std::vector<preset> (0);
    std::ifstream targetFile (filename.c_str());
    if (!targetFile.is_open()) throw std::runtime_error("No targets file found");
    std::string input_line, entity;
    while (std::getline(targetFile, input_line)){
        preset p;
        std::vector<unsigned int> family_data;
        std::stringstream line(input_line);
        std::getline(line, entity, ',');
        unsigned int i = std::stoi(entity);
        for (unsigned int k = 0; k < K_MAX; k++){
            std::getline(line, entity, ',');
            p.push_back(static_cast<status>(std::stoi(entity)));
            if (p[k] == COMPULSORY){
                presets_costs += CONSTANT_COST[k] + (*(this->family_data))[i][NB_CHOICES] * MARGINAL_COST[k];
                preset_occupancy[(*(this->family_data))[i][k]] += (*(this->family_data))[i][NB_CHOICES];
                preset_cardinal[(*(this->family_data))[i][k]] += 1;
                is_already_assigned[i] = true;
            }
        }
        std::getline(line, entity, ',');
        presets.push_back(p);
    }
}


unsigned int Presets::get_k(const unsigned int &family_index, const unsigned int &day_index) const {
    unsigned int k = 0;
    while ((*family_data)[family_index][k] != day_index && k < K_MAX)
        k++;
    if (k == K_MAX)
        throw std::invalid_argument("This day and this family have no link");
    return k;
}

void Presets::forbid_assignment(const unsigned int &i, const unsigned int &k, const bool &compute_stuff, const bool &k_is_relative) {
    if (!k_is_relative){
        unsigned int real_k = get_k(i, k);
        forbid_assignment(i, real_k, compute_stuff);
        return;
    }
    if (is_already_assigned[i])
        throw std::invalid_argument("We don't forbid assignments regarding previously assigned families");
    presets[i][k] = FORBIDDEN;
    if (compute_stuff) {
        compute_all_bounds();
        compute_feasibility();
    }
}

void Presets::enable_assignment(const unsigned int &i, const unsigned int &k, const bool &compute_stuff, const bool &k_is_relative) {
    if (!k_is_relative){
        unsigned int real_k = get_k(i, k);
        enable_assignment(i, real_k, compute_stuff);
        return;
    }
    if (presets[i][k] != FORBIDDEN)
        throw std::invalid_argument("This was not forbidden to begin with");
    if (is_already_assigned[i])
        throw std::invalid_argument("This is already assigned so we cannot enable other choices");
    presets[i][k] = ALLOWED;
    if (compute_stuff) {
        compute_all_bounds();
        compute_feasibility();
    }
}

void Presets::assign_family(const unsigned int &i, const unsigned int &k, const bool &compute_stuff, const bool &k_is_relative) {
    if (!k_is_relative){
        unsigned int real_k = get_k(i, k);
        assign_family(i, real_k, compute_stuff);
        return;
    }
    deassign_family(i, false);
    presets[i] = get_assignation_preset(k);
    presets_costs += CONSTANT_COST[k] + (*family_data)[i][NB_CHOICES] * MARGINAL_COST[k];
    preset_occupancy[(*family_data)[i][k]] += (*family_data)[i][NB_CHOICES];
    preset_cardinal[(*family_data)[i][k]] += 1;
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
            presets_costs -= CONSTANT_COST[k] + (*family_data)[i][NB_CHOICES] * MARGINAL_COST[k];
            preset_occupancy[(*family_data)[i][k]] -= (*family_data)[i][NB_CHOICES];
            preset_cardinal[(*family_data)[i][k]] -= 1;
            is_already_assigned[i] = false;
            nb_assignments--;
        }
    presets[i] = std::vector<status>(K_MAX, ALLOWED);
    if (compute_bounds){
        compute_all_bounds();
        compute_feasibility();
    }
}

void Presets::set_preset(const unsigned int &i, const preset &p, const bool &compute_stuff) {
    for (unsigned int k = 0; k < K_MAX; k++){
        if (p[k] == COMPULSORY) {
            assign_family(i, k);
            return;
        }
    }
    if (is_already_assigned[i])
        deassign_family(i);
    presets[i] = p;
    if (compute_stuff) {
        compute_all_bounds();
        compute_feasibility();
    }
}

std::string Presets::hash() const {
    std::vector<bool> S (NB_CHOICES * NB_FAMILIES, false);
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        for (unsigned int k = 0; k < K_MAX; k++)
            if (presets[i][k] == COMPULSORY)
                S[NB_CHOICES * i + k] = true;
    std::hash<std::vector<bool>> vector_hash;

    return std::to_string(vector_hash(S) % 100000000);
}

void Presets::write_solution(const std::string &filename) const {
    write_solution_(get_solution(), "../../solutions/" + filename + "_" + std::to_string(presets_costs + int(day_cost_lower_bound)) + "_" + hash() + ".csv");
}

std::vector<unsigned int> Presets::get_solution() const {
    std::vector<unsigned int> solution(0);
    unsigned int last_seen_allowed;
    unsigned int nb_allowed_seen;
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        nb_allowed_seen = 0;
        for (unsigned int k = 0; k < K_MAX; k++) {
            if (presets[i][k] == COMPULSORY)
                solution.push_back((*family_data)[i][k]);
            if (presets[i][k] == ALLOWED){
                nb_allowed_seen++;
                last_seen_allowed = (*family_data)[i][k];
            }
        }
        if (nb_allowed_seen == 1 && solution.size() == i)
            solution.push_back(last_seen_allowed);
        else if (solution.size() == i)
            throw std::logic_error("Do these presets define a day_sub_pb_solutions ?");
    }
    return solution;
}

void Presets::compute_feasibility() {
    is_feasible_ = true;
    for (unsigned int i = 0; i < NB_DAYS && is_feasible_; i++)
        if (occupancy_ub[i] < occupancy_lb[i])
            is_feasible_ = false;
    for (unsigned int i = 0; i < NB_DAYS && is_feasible_; i++)
        if (preset_occupancy[i] > occupancy_ub[i] || occupancy_ub[i] < occupancy_lb[i])
            is_feasible_ = false;
    for (unsigned int i = 0; i < NB_FAMILIES && is_feasible_; i++){
        if (is_already_assigned[i]) continue;
        unsigned int count = 0;
        for (unsigned int k = 0; k < K_MAX && !count; k++)
            if (presets[i][k] == ALLOWED)
                count++;
        if (!count)
            is_feasible_ = false;
    }
    if (!are_there_possible_quantities_per_day())
        is_feasible_ = false;
//    for (unsigned int i = 0; i < NB_FAMILIES && is_feasible_; i++){
//        bool b1 = false, b2 = false;
//        for (unsigned int k = 0; k < K_MAX; k++){
//            if (presets[i][k] == ALLOWED)
//                b1 = true;
//            if (presets[i][k] == COMPULSORY)
//                b2 = true;
//        }
//        if (b1 && b2)
//            throw std::logic_error("klvjd!#");
//    }
}

void Presets::compute_occupancy_bounds() {
    unsigned int max_day_cost = 500;
    occupancy_explicit_ub = std::vector<unsigned int> (NB_DAYS, 0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        for (unsigned int k = 0; k < K_MAX; k++)
            if (presets[i][k] != FORBIDDEN)
               occupancy_explicit_ub[(*family_data)[i][k]] += (*family_data)[i][NB_CHOICES];

    for (unsigned int i = 0; i < NB_DAYS; i++){
        occupancy_ub[i] = std::min(MAX_NB_PEOPLE_PER_DAY, occupancy_explicit_ub[i]);
        occupancy_lb[i] = std::max(MIN_NB_PEOPLE_PER_DAY, preset_occupancy[i]);
    }

    for (uint_pair &pair: prescribed_occupancy_ub){
        occupancy_ub[pair.first] = std::min(occupancy_ub[pair.first], pair.second);
    }

    for (uint_pair &pair: prescribed_occupancy_lb){
        occupancy_lb[pair.first] = std::max(occupancy_lb[pair.first], pair.second);
    }
//
//    for (int i = NB_DAYS - 1; i > 0; i--){
//        if (i < NB_DAYS - 1 && occupancy_ub[i] > occupancy_ub[i + 1] &&
//            (occupancy_ub[i] - MIN_NB_PEOPLE_PER_DAY)/400.*pow(occupancy_ub[i], 0.5+(occupancy_ub[i] - occupancy_ub[i + 1])/50.) > 500){
//            unsigned int a;
//            for (unsigned int j = 0; j < 10; j++) {
//                a = occupancy_ub[i];
//                occupancy_ub[i] = ceil(occupancy_ub[i + 1] + 50*(log(400.*max_day_cost/double(a - MIN_NB_PEOPLE_PER_DAY))/log(a) - 0.5));
//            }
//        }
//    }
//
//    for (unsigned int i = 0; i < NB_DAYS; i++)
//    if (i < NB_DAYS - 1 && occupancy_lb[i] == MIN_NB_PEOPLE_PER_DAY && occupancy_ub[i] <= occupancy_lb[i + 1] &&
//        1./400.*pow(occupancy_lb[i], 0.5+(occupancy_lb[i + 1] - occupancy_ub[i])/50.) > max_day_cost)
//        occupancy_ub[i] = MIN_NB_PEOPLE_PER_DAY;
//
//    for (unsigned int i = 0; i < NB_DAYS; i++){
//        if (i < NB_DAYS - 1 && occupancy_lb[i] > MIN_NB_PEOPLE_PER_DAY && occupancy_ub[i] <= occupancy_lb[i + 1]){// &&
//            //(occupancy_lb[i] - MIN_NB_PEOPLE_PER_DAY)/400.*pow(occupancy_lb[i], 0.5+(occupancy_lb[i + 1] - occupancy_lb[i])/50.) > 1000){
//            //std::cout << i << " " << occupancy_lb[i] << " " << occupancy_lb[i + 1] << " lower bound became ";
//            unsigned int a;
//            for (unsigned int j = 0; j < 10; j++) {
//                a = occupancy_lb[i];
//                occupancy_lb[i] = std::max(a, (unsigned int)(std::max(0., ceil(occupancy_lb[i + 1] - 50.*(log(400.*max_day_cost/double(a - MIN_NB_PEOPLE_PER_DAY))/log(a) - 0.5)))));
//            }
//            //std::cout << occupancy_lb[i] << std::endl;
//            a = occupancy_lb[i];
//            occupancy_ub[i + 1] = std::min(occupancy_ub[i + 1], (unsigned int)(floor(occupancy_ub[i] + 50.*(log(400.*max_day_cost/double(a - MIN_NB_PEOPLE_PER_DAY))/log(a) - 0.5))));
//            //std::cout << preset_lower_bound[i] << " " << occupancy_lb[i] << " " << occupancy_ub[i] << " " << occupancy_lb[i + 1] << " " << occupancy_ub[i + 1] << " " << preset_upper_bound[i + 1] << std::endl;
//        }
//    }
//
//    for (unsigned int i = 0; i < NB_DAYS - 1; i++){
//        // it's worth checking whether in this circumstance, minimizing the gap with day i + 1 always gives the minimum day cost...
//        unsigned int a = occupancy_ub[i];
//        unsigned int b = occupancy_lb[i + 1];
//        if (occupancy_lb[i] == MIN_NB_PEOPLE_PER_DAY && b > a){
//            float tmp = pow(a, 0.5 + (b - a) / 50.);
//            float tmp2 = (a - 125) / 400.;
//            float c = (a - 125) / 400. * pow(a, 0.5 + (b - a) / 50.);
//            if (c > 500)
//                occupancy_ub[i] = 125;
//        }
//    }

//    std::cout << "lower bounds:" << std::endl;
//    print_nicely(occupancy_lb, 5);

//    std::cout << "Upper bounds:" << std::endl;
//    print_nicely(occupancy_ub, 5);
}

double Presets::get_day_cost_lb_(const unsigned int &i) {
//    if (i == NB_DAYS - 1)
//        return get_binary_day_cost_lb(occupancy_lb[i], occupancy_ub[i], occupancy_lb[i], occupancy_ub[i]);
//    else
//        return get_binary_day_cost_lb(occupancy_lb[i], occupancy_ub[i], occupancy_lb[i + 1], occupancy_ub[i + 1]);
    if (i == 0)
        return 0.5 * get_binary_day_cost_lb(occupancy_lb[i], occupancy_ub[i], occupancy_lb[i + 1], occupancy_ub[i + 1]);
    else if (i == NB_DAYS - 1)
        return 0.5 * get_binary_day_cost_lb(occupancy_lb[i - 1], occupancy_ub[i - 1], occupancy_lb[i], occupancy_ub[i]) +
               get_binary_day_cost_lb(occupancy_lb[i], occupancy_ub[i], occupancy_lb[i], occupancy_ub[i]);
    else
        return get_ternary_day_cost_lb(occupancy_lb[i - 1], occupancy_ub[i - 1], occupancy_lb[i],
                                       occupancy_ub[i], occupancy_lb[i + 1], occupancy_ub[i + 1]);
}

double Presets::get_additional_day_cost_lb(const unsigned int &i) {//return 0;
    // i is a family index. we compute a lower bound on the day cost of assigning the family.
    if (is_already_assigned[i]) return 0;
    std::unordered_set<unsigned int> impacted_days(0);
    for (unsigned int k = 0; k < K_MAX; k++){
        impacted_days.emplace((*family_data)[i][k]);
        if (int((*family_data)[i][k]) - 1 >= 0)
            impacted_days.emplace((*family_data)[i][k] - 1);
        if ((*family_data)[i][k] < NB_DAYS - 1)
            impacted_days.emplace((*family_data)[i][k] + 1);
    }
    double min_possible_increase = 10000;
    for (unsigned int k_assign = 0; k_assign < K_MAX; k_assign++){
        if (preset_occupancy[(*family_data)[i][k_assign]] + (*family_data)[i][NB_CHOICES] > occupancy_ub[(*family_data)[i][k_assign]])
            continue;
        if (presets[i][k_assign] == FORBIDDEN) continue;
        // simulate assigning to k
        std::vector<unsigned int> occupancy_ub_save(K_MAX);
        std::vector<unsigned int> occupancy_lb_save(K_MAX);
        for (unsigned int k = 0; k < K_MAX; k++){
            const unsigned int &current_day = (*family_data)[i][k];
            occupancy_lb_save[k] = occupancy_lb[current_day];
            occupancy_ub_save[k] = occupancy_ub[current_day];
            if (k == k_assign)
                occupancy_lb[current_day] = std::max(occupancy_lb[current_day], preset_occupancy[current_day] + (*family_data)[i][NB_CHOICES]);
            else if (presets[i][k] != FORBIDDEN)
                occupancy_ub[current_day] = std::min(occupancy_ub[current_day], occupancy_explicit_ub[current_day] - (*family_data)[i][NB_CHOICES]);
        }
        // compute the increase
        double possible_increase = 0;
        for (const unsigned int& current_day: impacted_days) {
            possible_increase += get_day_cost_lb_(current_day) - day_costs_lb[current_day];
//            if (get_day_cost_lb_(current_day) < day_costs_lb[current_day])
//                throw std::logic_error("impossible " + std::to_string(day_costs_lb[current_day]) + " " + std::to_string(get_day_cost_lb_(current_day)));
        }
        if (possible_increase < min_possible_increase)
            min_possible_increase = possible_increase;
        // restore the previous bounds
        for (int k = K_MAX - 1; k >= 0; k--){
            const unsigned int &current_day = (*family_data)[i][k];
            occupancy_lb[current_day] = occupancy_lb_save[k];
            occupancy_ub[current_day] = occupancy_ub_save[k];
        }
    }
    return min_possible_increase;
}

void Presets::compute_all_bounds(const bool & including_costs) {
    clock_t t0 = clock();
    compute_occupancy_bounds();
//    std::cout << "time spent on occupancy_bounds: " << clock() - t0 << std::endl;

    t0 = clock();
//    std::cout << "Nb of families preset and number of people" << std::endl;
//    print_nicely(preset_occupancy, 5);
    //print_nicely(preset_cardinal, 5);
    compute_bottleneck_bounds();
//    std::cout << "time spent on bottleneck bounds: " << clock() - t0 << std::endl;
//    t0 = clock();
    compute_k_fold_bottleneck_ub();
//    std::cout << "time spent on k-fold bottleneck bounds: " << clock() - t0 << std::endl;
    if (including_costs) {
        t0 = clock();
        for (unsigned int i = 0; i < NB_DAYS; i++)
            day_costs_lb[i] = get_day_cost_lb_(i);
        day_cost_lower_bound = 0;
        for (unsigned int i = 0; i < NB_DAYS; i++)
            day_cost_lower_bound += day_costs_lb[i];
//        std::cout << "time spent on computing the day costs lower bound: " << clock() - t0 << std::endl;
        t0 = clock();
        if (should_we_compute_additional_day_costs()) {
//            std::cout << "day cost lower bound before we do additional day costs is " << day_cost_lower_bound
//                      << std::endl;
            for (unsigned int i = 0; i < NB_FAMILIES; i++)
                day_cost_lower_bound += get_additional_day_cost_lb(i);
        }
//        std::cout << "time spent on computing the per-family additional day costs: " << clock() - t0 << std::endl;
//        std::cout << "local lower bounds on the costs from days yield " << day_cost_lower_bound << std::endl;
        t0 = clock();
       // if (should_we_compute_day_costs_with_DP()) {
            compute_day_cost_DP_lb();
//            std::cout << "time spent on computing the day costs lower bound with DP: " << clock() - t0 << std::endl;
            day_cost_lower_bound = std::max(day_cost_lower_bound, day_cost_DP_lb);
        //}
//        std::cout << "day costs final lower bound is " << day_cost_lower_bound;

//        if (day_cost_lower_bound < 1500)
//            std::cout << " before being set to 1500";
//        std::cout << std::endl;
    }
}

bool Presets::should_we_compute_additional_day_costs() const {
    if (NB_FAMILIES - nb_assignments > 25) {
//        std::cout << "we do not compute the additional day costs because there are more than 25 families left" << std::endl;
        return false;
    }
    for (unsigned int i = 0; i < NB_DAYS; i++){
        if (preset_occupancy[i] + 10 > occupancy_ub[i] || occupancy_explicit_ub[i] - 10 < occupancy_lb[i])
            return true;
    }
    return false;
}

bool Presets::should_we_compute_day_costs_with_DP() const {
    if (NB_FAMILIES - nb_assignments > 60) {
//        std::cout << "we do not compute the day costs with DP because there are more than 60 families left" << std::endl;
        return false;
    }
    return true;
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
            current_day = (*family_data)[current_family][k];
            if (presets[current_family][k] == ALLOWED &&
                possible_family_cum_sizes[current_day] + (*family_data)[current_family][NB_CHOICES] <= occupancy_ub[current_day] - preset_occupancy[current_day]){
                possible_family_cum_sizes[current_day] += (*family_data)[current_family][NB_CHOICES];
                bottleneck_bounds[current_day].second++;
            }
        }
    }
    possible_family_cum_sizes = std::vector<unsigned int>(NB_DAYS, 0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        for (unsigned int k = 0; k < K_MAX; k++) {
            current_family = sorted_families[NB_FAMILIES - i - 1].first;
            current_day = (*family_data)[current_family][k];
            if (presets[current_family][k] == ALLOWED &&
                possible_family_cum_sizes[current_day] < occupancy_lb[current_day] - preset_occupancy[current_day]){
                possible_family_cum_sizes[current_day] += (*family_data)[current_family][NB_CHOICES];
                bottleneck_bounds[current_day].first++;
            }
        }
    }
//    std::cout << "bottleneck bounds (min and then max):" << std::endl;
//    for (unsigned int i = 0; i < NB_DAYS; i++) {
//        unsigned int n = nb_chiffres(bottleneck_bounds[i].first);
//        if (n >= 4)
//            throw std::logic_error("Why is this integer so big ?");
//        std::cout << bottleneck_bounds[i].first << std::string(5 - n, ' ');
//    }
//    std::cout << std::endl;
//    for (unsigned int i = 0; i < NB_DAYS; i++) {
//        unsigned int n = nb_chiffres(bottleneck_bounds[i].second);
//        if (n >= 4)
//            throw std::logic_error("Why is this integer so big ?");
//        std::cout << bottleneck_bounds[i].second << std::string(5 - n, ' ');
//    }
//    std::cout << std::endl;
}

void Presets::compute_k_fold_bottleneck_ub() {
    unsigned int current_fam, current_day;
    k_fold_bottleneck_bounds = std::vector<std::vector<unsigned int>>(NB_DAYS, std::vector<unsigned int>(K_MAX, 0));
    std::vector<std::vector<unsigned int>> pos_fam_cum_sizes(NB_DAYS, std::vector<unsigned int>(K_MAX, 0));
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        for (unsigned int k = 0; k < K_MAX; k++){
            current_fam = sorted_families[i].first;
            current_day = (*family_data)[current_fam][k];
            if (presets[current_fam][k] == ALLOWED &&
                pos_fam_cum_sizes[current_day][k] + (*family_data)[current_fam][NB_CHOICES] <= occupancy_ub[current_day] - preset_occupancy[current_day]){
                pos_fam_cum_sizes[current_day][k] += (*family_data)[current_fam][NB_CHOICES];
                k_fold_bottleneck_bounds[current_day][k]++;
            }
        }
    }
}

std::vector<unsigned int> Presets::get_largest_unassigned_families() const {
    std::vector<unsigned int> largest(0);
    unsigned int maxi = 1;
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        if (!is_already_assigned[i] && (*family_data)[i][NB_CHOICES] > maxi){
            maxi = (*family_data)[i][NB_CHOICES];
            largest.clear();
            largest.push_back(i);
        }
        if (!is_already_assigned[i] && (*family_data)[i][NB_CHOICES] == maxi)
            largest.push_back(i);
    }
    return largest;
}

uint_pair Presets::get_bounds_to_branch() const {
    unsigned int max_gap = 0;
    unsigned int lower_bound, i_max;
    std::vector<unsigned int> gaps;
    for (unsigned int i = 0; i < NB_DAYS; i++){
        unsigned int gap = occupancy_ub[i] - occupancy_lb[i];
        gaps.push_back(gap);
        if (gap > max_gap){
            max_gap = gap;
            i_max = i;
            lower_bound = occupancy_lb[i];
        }
    }
    std::cout << "gaps between occupancy bounds" << std::endl;
    print_nicely(gaps, 5);
    return uint_pair(i_max, lower_bound + floor(max_gap/2));
}

std::vector<unsigned int> Presets::get_all_families_assigned_to_day(const unsigned int &day) const {
    std::vector<unsigned int> res(0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        if (!is_already_assigned[i])
            continue;
        for (unsigned int k = 0; k < K_MAX; k++)
            if (presets[i][k] == COMPULSORY && (*family_data)[i][k] == day)
                res.push_back(i);
    }
    return res;
}

unsigned int Presets::get_presets_costs(const unsigned int &day) const {
    unsigned int sum_costs = 0;
    unsigned int sum_occup = 0;
    for (unsigned int i = 0; i < NB_FAMILIES && sum_occup < preset_occupancy[day]; i++){
        if (!is_already_assigned[i])
            continue;
        for (unsigned int k = 0; k < K_MAX; k++)
            if (presets[i][k] == COMPULSORY && (*family_data)[i][k] == day){
            sum_occup += (*family_data)[i][NB_CHOICES];
            sum_costs += (*family_data)[i][NB_CHOICES] * MARGINAL_COST[k] + CONSTANT_COST[k];
        }
    }
    // check the families we have covered account for the preset occupancy
    if (sum_occup == preset_occupancy[day])
        return sum_costs;
    else
        throw std::logic_error("There is a problem");
}

std::vector<uint_pair> Presets::get_assignments_by_default() const {
    std::vector<uint_pair> res;
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        if (is_already_assigned[i]) continue;
        unsigned int count = 0;
        unsigned int last_allowed;
        for (unsigned int k = 0; k < K_MAX && count <= 1; k++){
            if (presets[i][k] == ALLOWED){
                count++;
                last_allowed = k;
            }
        }
        if (count == 1)
            res.push_back(uint_pair(i, last_allowed));
    }
    return res;
}

unsigned int Presets::get_family_max_size_min_choices() const {
    std::vector<unsigned int> largest_families = get_largest_unassigned_families();
    unsigned int mini = K_MAX;
    unsigned int i_mini;
    for (unsigned int &i: largest_families){
        unsigned int count = 0;
        for (unsigned int k = 0; k < K_MAX; k++)
            if (presets[i][k] == ALLOWED)
                count++;
        if (count < mini){
            mini = count;
            i_mini = i;
        }
    }
    return i_mini;
}

void Presets::compute_day_cost_DP_lb() {
    std::vector<std::vector<unsigned int>> possible_quantities_per_day = get_possible_quantities_per_day();

    unsigned int s = 0;
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        s += (*family_data)[i][NB_CHOICES];

//    if (prescribed_occupancy_ub.size() > 0)
//        get_day_cost_DP_ultimate_lb(possible_quantities_per_day, s);

    day_cost_DP_lb =  get_day_cost_DP_lb(possible_quantities_per_day);
//    std::cout << "The bound on the day costs computed with DP is " << day_cost_DP_lb << std::endl;
}

std::vector<std::vector<unsigned int>> Presets::get_possible_quantities_per_day() const {
    std::vector<std::vector<unsigned int>> possible_sizes_per_day(NB_DAYS, std::vector<unsigned int>(0));
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        for (unsigned int k = 0; k < K_MAX; k++)
            if (presets[i][k] == ALLOWED)
                possible_sizes_per_day[(*family_data)[i][k]].push_back((*family_data)[i][NB_CHOICES]);
    std::vector<std::vector<unsigned int>> possible_quantities_per_day(NB_DAYS, std::vector<unsigned int>(0));
    for (unsigned int i = 0; i < NB_DAYS; i++) {
        unsigned int min_allowed = (unsigned int)(std::max(int(get_occupancy_lb(i)) - int(get_presets_occupancy(i)), 0));
        unsigned int max_allowed = get_occupancy_ub(i) - get_presets_occupancy(i);
        //if (max_allowed > 1000) // unsigned int overflow, doesn't matter because unfeasible.
          //  return;
        possible_quantities_per_day[i] = get_possible_quantities(possible_sizes_per_day[i], min_allowed, max_allowed);
    }

    for (unsigned int i = 0; i < NB_DAYS; i++)
        for (unsigned int m = 0; m < possible_quantities_per_day[i].size(); m++)
            possible_quantities_per_day[i][m] += preset_occupancy[i];

    return possible_quantities_per_day;
}

Presets &Presets::operator=(const Presets &other) {
    //family_data = other.family_data;
    presets = other.presets;
    is_already_assigned = other.is_already_assigned;
    preset_occupancy  = other.preset_occupancy;
    preset_cardinal   = other.preset_cardinal;
    occupancy_lb   = other.occupancy_lb;
    occupancy_ub    = other.occupancy_ub ;
    occupancy_explicit_ub = other.occupancy_explicit_ub;
    prescribed_occupancy_ub = other.prescribed_occupancy_ub;
    prescribed_occupancy_lb  = other.prescribed_occupancy_lb;
    bottleneck_bounds   = other.bottleneck_bounds;
    k_fold_bottleneck_bounds   = other.k_fold_bottleneck_bounds;
    day_costs_lb  = other.day_costs_lb;
    day_cost_DP_lb  = other.day_cost_DP_lb;
    presets_costs  = other.presets_costs;
    day_cost_lower_bound    = other.day_cost_lower_bound;
    is_feasible_    = other.is_feasible_;
    nb_assignments = other.nb_assignments;
    sorted_families = other.sorted_families;
}

bool Presets::are_there_possible_quantities_per_day() const {
    std::vector<std::vector<unsigned int>> possible_qties_per_day = get_possible_quantities_per_day();
    bool is_there_empty = false;
    for (unsigned int j = 0; j < NB_DAYS && !is_there_empty; j++)
        if (possible_qties_per_day[j].size() == 0)
            is_there_empty = true;
    return !is_there_empty;
}

void Presets::write_presets(const std::string &filename) {
    std::ofstream write(filename.c_str());
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        write << i << ",";
        for (unsigned int k = 0; k < K_MAX; k++){
            write << presets[i][k] << ",";
        }
        write << std::endl;
    }
}

