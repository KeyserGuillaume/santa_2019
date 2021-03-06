#include "Assignment.h"
#include <math.h>
#include <iostream>

bool Day::has_removable_family() const {
    if (N < MIN_NB_PEOPLE_PER_DAY + 10 || N > 300) {
        for (unsigned int i = 0; i < assigned_families.size(); i++) {
            if (N - assigned_families[i]->get_nb_people() >= MIN_NB_PEOPLE_PER_DAY &&
                N - assigned_families[i]->get_nb_people() <= MAX_NB_PEOPLE_PER_DAY)
                return true;
        }
        return false;
    }
    else{
        return true;
    }
}

Family *Day::get_random_family() {
    unsigned int j = rand() % assigned_families.size();
    Family* tmp = assigned_families[j];
    assigned_families[j] = assigned_families[assigned_families.size() - 1];
    assigned_families[assigned_families.size() - 1] = tmp;
    return tmp;
}


Family *Day::get_random_removable_family(const bool & is_Friday_like) {
    unsigned int j = 0;
    if (N < MIN_NB_PEOPLE_PER_DAY + 10 || N > 300) {
        std::vector<unsigned int> possible_indexes(0);
        for (unsigned int i = 0; i < assigned_families.size(); i++) {
            if (N - assigned_families[i]->get_nb_people() >= MIN_NB_PEOPLE_PER_DAY &&
                N - assigned_families[i]->get_nb_people() <= MAX_NB_PEOPLE_PER_DAY &&
                (!is_Friday_like || N - assigned_families[i]->get_nb_people() == MIN_NB_PEOPLE_PER_DAY)) {
                possible_indexes.push_back(i);
            }
        }
        j = possible_indexes[rand() % possible_indexes.size()];
    }
    else{
        j = rand() % assigned_families.size();
    }
    Family* tmp = assigned_families[j];
    assigned_families[j] = assigned_families[assigned_families.size() - 1];
    assigned_families[assigned_families.size() - 1] = tmp;
    return tmp;
}

void Day::compute_cost() {
    cost = is_feasible() ? ceil(1 * (N - MIN_NB_PEOPLE_PER_DAY)/400.*pow(N, 0.5+(abs(previous_day->get_N() - N)/50.))) : 0;
    //cost = 10*floor(cost/10.);
    // do not try to change 0 with cost (meaning the cost of an unfeasible day would be the cost when last it was feasible)
    // because then the twin path move will no longer keep the right count of the cost variation (because trying to insert
    // a family before removing it can make it temporarily feasible with a huge cost (ie if N = 127 and Friday_like),
    // and then when you set a family to fall back on a reasonable value, you did not expect to get the hugely negative
    // minus previous cost.
}

int Day::remove_family(Family *f) {
    for (int i = assigned_families.size() - 1; i >= 0; i--){
        if (assigned_families[i]->get_id() == f->get_id()){
            assigned_families[i] = assigned_families[assigned_families.size() - 1];
            assigned_families.pop_back();
            N -= f->get_nb_people();
            if (next_day->get_id() == id){
                unsigned int previous_cost = cost;
                compute_cost();
                return cost - previous_cost;
            }
            else{
                unsigned int previous_cost = cost + next_day->get_cost();
                compute_cost();
                next_day->compute_cost();
                return cost + next_day->get_cost() - previous_cost;
            }
        }
    }
    throw std::logic_error("Tried to remove a family from a day it was not assigned to");
}

int Day::add_family(Family *family) {
    assigned_families.push_back(family);
    N += family->get_nb_people();
    if (next_day->get_id() == id) {
        unsigned int previous_cost = cost;
        compute_cost();
        return cost - previous_cost;
    } else {
        unsigned int previous_cost = cost + next_day->get_cost();
        compute_cost();
        next_day->compute_cost();
        return cost + next_day->get_cost() - previous_cost;
    }
}

unsigned int Day::get_assignments_costs() const {
    unsigned int res = 0;
    for (unsigned int j = 0; j < assigned_families.size(); j++)
        res += get_ith_family(j)->get_cost();
    return res;
}

Family::Family(const unsigned int &id,
               const unsigned int &n_people,
               const std::vector<Day*> &preferred_days_,
               Day* assigned_day): id(id), n_people(n_people), assigned_day(assigned_day) {
    preferred_days.clear();
    for (unsigned int i = 0; i < preferred_days_.size(); i++){
        preferred_days.push_back(preferred_days_[i]);
    }
    compute_cost();
}

void Family::compute_cost() {
    k = 0;
    bool found = false;
    for (; k < NB_CHOICES && !found; k++)
        if (assigned_day->get_id() == preferred_days[k]->get_id()) {
            found = true;
            k--;
        }
    cost = CONSTANT_COST[k] + n_people*MARGINAL_COST[k];
}

std::pair<int, int> Family::set_assigned_day(Day* d) {
    //std::cout << "setting family " << id << " from day " << assigned_day->get_id() << " to day " << d->get_id()<<std::endl;
    int cost_var = d->add_family(this) + assigned_day->remove_family(this);
    unsigned int previous_cost = cost;
    assigned_day = d;
    compute_cost();
    return std::pair<int, int>(cost - previous_cost, cost_var);
}

Day *Family::get_best_possible_day() const {
    for (unsigned int i = 0; i < k; i++){
        Day* res = preferred_days[i];
        if (res->get_N() + n_people <= MAX_NB_PEOPLE_PER_DAY)
            return res;
    }
    return preferred_days[k];
}

Day *Family::get_random_preferred_day_within_threshold(const unsigned int &threshold) const {
    if (k == NB_CHOICES) return preferred_days[rand() % NB_CHOICES];
    for (unsigned int i = k + 1; i < NB_CHOICES; i++)
        if (CONSTANT_COST[i] + n_people * MARGINAL_COST[i] - cost > 1.5*threshold)
            return preferred_days[rand() % i];
    // we get here when the family has their last preferred day
    return preferred_days[rand() % NB_CHOICES];
}

bool Family::is_removable() {
    return (assigned_day->get_N() - n_people >= MIN_NB_PEOPLE_PER_DAY);
}

Assignment::Assignment(const std::vector<std::vector<unsigned int>> &family_data, const std::vector<unsigned int> &solution) {
    families = new Family[NB_FAMILIES];
    days = new Day[NB_DAYS];

    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        std::vector<Day*> preferred_days(0);
        for (unsigned int k = 0; k < NB_CHOICES; k++){
            preferred_days.push_back(days + family_data[i][k]);
        }
        families[i] = Family(
                i,
                family_data[i][NB_CHOICES],
                //std::vector<unsigned int>(std::vector<unsigned int>(family_data[i].begin(), family_data[i].end() - 1)),
                //std::vector<unsigned int>(family_data[i].begin(), family_data[i].end() - 1),
                preferred_days,
                days + solution[i]);
    }

    for (unsigned int j = 0; j < NB_DAYS - 1; j++)
        days[j].set_previous_day(days + j + 1);

    for (unsigned int j = 1; j < NB_DAYS; j++)
        days[j].set_next_day(days + j - 1);

    for (unsigned int j = 0; j < NB_DAYS; j++)
        days[j].set_id(j);

    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        days[solution[i]].add_family(families + i);

    get_cost(); // initialize the k of families
}

void Assignment::write_solution(const std::string &filename) const {
    std::vector<unsigned int> solution(NB_FAMILIES);
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        solution[i] = families[i].get_assigned_day()->get_id();
    write_solution_(solution, filename);
}

unsigned int Assignment::get_cost(){
    unsigned int res1 = 0, res2 = 0;
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        families[i].compute_cost();
        res1 += families[i].get_cost();
    }
    for (unsigned int j = 0; j < NB_DAYS; j++) {
        days[j].compute_cost();
        res2 += days[j].get_cost();
    }
    accounting_costs = res2;
    //std::cout << res1 << " " << res2 << std::endl;
    return res1 + res2;
}

void Assignment::check_solution_is_ok() {
    for (unsigned int j = 0; j < NB_DAYS; j++){
        if (!days[j].is_feasible()){
            throw std::logic_error("Day " + std::to_string(j) + " has bad number of people : " + std::to_string(days[j].get_N()));
        }
    }
}

void Assignment::stats() const {
    unsigned int res1 = 0, res2 = 0;
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        families[i].compute_cost();
        res1 += families[i].get_cost();
    }
    for (unsigned int j = 0; j < NB_DAYS; j++) {
        days[j].compute_cost();
        res2 += days[j].get_cost();
    }
    std::cout << "total cost: " << res1 + res2 << std::endl;
    std::cout << get_nb_Friday_like() << " Friday-likes in the solution" << std::endl;
    std::cout << "Families Costs : " << res1 << ";   Accounting Costs : " << res2 << std::endl << std::endl;

    std::vector<unsigned int> count(NB_CHOICES + 1, 0);
    std::vector<unsigned int> cost(NB_CHOICES + 1, 0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        count[families[i].get_k()]++;
        cost[families[i].get_k()]+=families[i].get_cost();
    }
    int max_preference = NB_CHOICES;
    std::cout << "  stats: histogram on the preference level" << std::endl;
    for (unsigned int i=0; i < max_preference + 1; i++)
        std::cout << i << "       ";
    std::cout << std::endl;
    for (unsigned int i=0; i < max_preference + 1; i++){
        std::cout << count[i] << std::string(7-nb_chiffres(count[i]) + nb_chiffres(i), ' ');
    }
    std::cout << std::endl;
    for (unsigned int i=0; i < max_preference + 1; i++){
        std::cout << cost[i] << std::string(7-nb_chiffres(cost[i]) + nb_chiffres(i), ' ');
    }
    std::cout << std::endl << std::endl;

    std::vector<unsigned int> count_size(9, 0);
    std::vector<unsigned int> cost_size(9, 0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        count_size[families[i].get_nb_people()]++;
        cost_size[families[i].get_nb_people()]+=families[i].get_cost();
    }
    std::cout << "  stats: histogram on the size of the family" << std::endl;
    for (unsigned int i=0; i < 9; i++)
        std::cout << i << "       ";
    std::cout << std::endl;
    for (unsigned int i=0; i < 9; i++){
        std::cout << count_size[i] << std::string(7-nb_chiffres(count_size[i]) + nb_chiffres(i), ' ');
    }
    std::cout << std::endl;
    for (unsigned int i=0; i < 9; i++){
        std::cout << cost_size[i] << std::string(7-nb_chiffres(cost_size[i]) + nb_chiffres(i), ' ');
    }
    std::cout << std::endl << std::endl;

    // the cost of the families of every day
    std::vector<unsigned int> tmp(NB_DAYS, 0);
    for (unsigned int i = 0; i < NB_DAYS; i++)
        tmp[i] = days[i].get_assignments_costs();
    print_nicely(tmp, 5);

    for (unsigned int i = 0; i < NB_DAYS; i++)
         tmp[i] = days[i].get_N();
    print_nicely(tmp, 5);


//    std::vector<unsigned int> count_crowdedness(301, 0);
//    std::vector<unsigned int> cost_crowdedness(301, 0);
//    for (unsigned int i = 0; i < NB_DAYS; i++) {
//        count_crowdedness[days[i].get_N()]++;
//        for (unsigned int m = 0; m < days[i].get_nb_families(); m++) {
//            cost_crowdedness[days[i].get_N()] += days[i].get_ith_family(m)->get_cost();
//        }
//    }
//    std::cout << "  stats: histogram on the crowdedness of the day" << std::endl;
//    for (unsigned int i=0; i < 301; i++)
//        if (count_crowdedness[i] > 0)
//            std::cout << i << "       ";
//    std::cout << std::endl;
//    for (unsigned int i=0; i < 301; i++)
//        if (count_crowdedness[i] > 0)
//            std::cout << count_crowdedness[i] << std::string(7-nb_chiffres(count_crowdedness[i]) + nb_chiffres(i), ' ');
//    std::cout << std::endl;
//    for (unsigned int i=0; i < 301; i++)
//        if (count_crowdedness[i] > 0)
//            std::cout << cost_crowdedness[i] << std::string(7-nb_chiffres(cost_crowdedness[i]) + nb_chiffres(i), ' ');
//    std::cout << std::endl;
}

Day *Assignment::get_random_Friday_like() const {
    std::vector<Day*> possible_results(0);
    for (unsigned int i = 0; i < NB_DAYS; i++)
        if (days[i].is_Friday_like())
            possible_results.push_back(days + i);
    return possible_results[rand() % possible_results.size()];
}

bool Assignment::has_Friday_like() const {
    for (unsigned int i = 0; i < NB_DAYS; i++)
        if (days[i].is_Friday_like())
            return true;
    return false;
}

unsigned int Assignment::get_nb_Friday_like() const {
    unsigned int count = 0;
    for (unsigned int i = 0; i < NB_DAYS; i++)
        if (days[i].is_Friday_like())
            count++;
    return count;
}

int Assignment::set_assigned_day(Family *f, Day *i) {
    unsigned int prev_acc_cost = accounting_costs;
    std::pair<int, int> my_pair = f->set_assigned_day(i);
    accounting_costs += my_pair.second;
    int accounting_costs_limit = 3000;
    if (accounting_costs < accounting_costs_limit && prev_acc_cost < accounting_costs_limit)
       return my_pair.first + my_pair.second;
    else if (accounting_costs < accounting_costs_limit)
        return my_pair.first + my_pair.second - 1000;
    else
        return my_pair.first + my_pair.second + 1000;
}

int Assignment::set_assigned_day(const unsigned int &i, const unsigned int &j) {
    return set_assigned_day(families + i, days + j);
}
