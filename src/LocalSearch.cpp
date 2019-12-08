#include "LocalSearch.h"

void LocalSearch::run_on_time_limit(const clock_t time_limit, const int &period_display) {
    int k0 = k;
    bool has_exceeded_time_limit = false;
    for (; !has_exceeded_time_limit; k++) {
        jump();
        if ((k - k0)%period_display==0) {
            display();
            A->check_solution_is_ok(); // raises errors if pb
        }
        if ((k - k0) % 100 == 0) {
            if (clock() > time_limit)
                has_exceeded_time_limit = true;
        }
    }
}

void LocalSearch::run(const long &nb_iteration, const int &period_display) {
    int k0 = k;
    for (; k < k0 + nb_iteration; k++){
        jump();
        if ((k - k0)%period_display==0) {
            display();
            A->check_solution_is_ok(); // raises errors if pb
        }
    }
}

void LocalSearch::jump(){
    //augmenting_path_move();
    twin_paths_move();
}

void LocalSearch::twin_paths_move(){
    int abort_threshold = 100;
    int current_cost_variation = 0;
    int temp_cost_variation;
    bool has_a_path_moved_to_first = false;
//    std::vector<bool> day_was_visited(NB_DAYS, false);
    std::vector<Family*> visited_families_A(0);
    std::vector<Family*> visited_families_B(0);
    std::vector<Day*> visited_days_A(0);
    std::vector<Day*> visited_days_B(0);
    std::vector<Day*> * visited_days;
    std::vector<Family*> * visited_families;
    Day* current_day;
    Family* current_family;

    Day* first_day = A->get_random_day();
    Family* big_family = first_day->get_random_family();

    if (big_family->get_k() == 0) return;

    current_day = big_family->get_random_improving_day();
    visited_days_A.push_back(current_day);
    visited_days_B.push_back(current_day);
    current_cost_variation += big_family->set_assigned_day(current_day);

    do{
        current_family = current_day->get_random_family();
    } while(current_family->get_id() == big_family->get_id());
    temp_cost_variation = current_family->set_assigned_day(first_day);
    if (temp_cost_variation < abort_threshold) {
        visited_days_A.push_back(first_day);
        current_cost_variation += temp_cost_variation;
        has_a_path_moved_to_first = true;
    }
    else {
        current_family->set_assigned_day(current_day);
        current_day = current_family->get_random_preferred_day_within_threshold(abort_threshold);
        current_cost_variation += current_family->set_assigned_day(current_day);
        visited_days_A.push_back(current_day);
        if (current_day->get_id() == first_day->get_id())
            has_a_path_moved_to_first = true;
    }
    visited_families_A.push_back(current_family);

    current_day = big_family->get_assigned_day();
    if (!current_day->has_removable_family()) {
        current_family->set_assigned_day(current_day);
        big_family->set_assigned_day(first_day);
        A->check_solution_is_ok();
        return;
    }
    unsigned int c = 0;
    do {current_family = current_day->get_random_removable_family();
    c++;
    } while (current_family->get_id() == big_family->get_id() && c < 10);
    temp_cost_variation = current_family->set_assigned_day(first_day);
    if (temp_cost_variation < abort_threshold && (!has_a_path_moved_to_first || first_day->is_feasible())) {
        visited_days_B.push_back(first_day);
        current_cost_variation += temp_cost_variation;
        has_a_path_moved_to_first = true;
    }
    else {
        current_family->set_assigned_day(current_day);
        current_day = current_family->get_random_preferred_day_within_threshold(abort_threshold);
        current_cost_variation += current_family->set_assigned_day(current_day);
        visited_days_B.push_back(current_day);
        if (current_day->get_id() == first_day->get_id())
            has_a_path_moved_to_first = true;
    }
    visited_families_B.push_back(current_family);

    bool do_path_A = false; // we really do start with path_A

    while(current_cost_variation < abort_threshold && (
            !first_day->is_feasible() ||
            !visited_days_A[visited_days_A.size() - 1]->is_feasible() ||
            !visited_days_B[visited_days_B.size() - 1]->is_feasible())){
        do_path_A = 1 - do_path_A;
        if (do_path_A){
            visited_days = &visited_days_A;
            visited_families = &visited_families_A;
        } else {
            visited_days = &visited_days_B;
            visited_families = &visited_families_B;
        }
        current_day = (*visited_days)[visited_days->size() - 1];
        if (current_day->get_id() == first_day->get_id())
            continue;

        if (current_day->is_feasible() && !has_a_path_moved_to_first && (
                (do_path_A && !visited_days_B[visited_days_B.size() - 1]->is_feasible()) ||
                (!do_path_A && !visited_days_A[visited_days_A.size() - 1]->is_feasible())))
            continue;

        if (current_day->has_removable_family())
           current_family = current_day->get_random_removable_family();
        else
            current_family = current_day->get_random_family();
        // check that assigning to first day is not strictly improving
        temp_cost_variation = current_family->set_assigned_day(first_day);
        if (temp_cost_variation < abort_threshold && (first_day->is_feasible() || !has_a_path_moved_to_first)) {
            current_cost_variation += temp_cost_variation;
            visited_days->push_back(first_day);
            has_a_path_moved_to_first = true;
        } else {
            current_family->set_assigned_day(current_day);
            current_day = current_family->get_random_preferred_day_within_threshold(abort_threshold);
            current_cost_variation += current_family->set_assigned_day(current_day);
            visited_days->push_back(current_day);
            if (current_day->get_id() == first_day->get_id())
                has_a_path_moved_to_first = true;
        }
        visited_families->push_back(current_family);
    }

    if (current_cost_variation > 0){
        // undo what we did
        unsigned int path_A_length = visited_families_A.size();
        unsigned int path_B_length = visited_families_B.size();
        for (int i = std::max(path_A_length - 1, path_B_length - 1); i >= 0; i--) {
            if (i < path_B_length)
                visited_families_B[i]->set_assigned_day(visited_days_B[i]);
            if (i < path_A_length)
                visited_families_A[i]->set_assigned_day(visited_days_A[i]);
            //std::cout << "setting family " << visited_families[i]->get_id() << " back to day " << visited_days[i]->get_id() << std::endl;
        }
        big_family->set_assigned_day(first_day);
    }

    A->check_solution_is_ok();
}

void LocalSearch::augmenting_path_move() {
    int abort_threshold = 100;

    std::vector<Day*> visited_days(0);
    std::vector<Family*> visited_families(0);
    bool ongoing = true;
    int current_cost_variation = 0;
    int temp_cost_variation;

    Day* current_day;
    do {
        current_day = A->get_random_day();
    } while(!current_day->has_removable_family());
    Day* first_day = current_day;
    Family* current_family = current_day->get_random_removable_family();

    while (ongoing){
        visited_days.push_back(current_day);
        if (visited_days.size() > 1) {
            // check that assigning to first day is not strictly improving
            temp_cost_variation = current_family->set_assigned_day(first_day);
            if (first_day->is_feasible() && current_cost_variation + temp_cost_variation < 0) {
                if (visited_days.size() + 1 < count_order_change.size())
                    count_order_change[visited_days.size() + 1]++;
                nb_successful_augmenting_path_2++;
                return;
            }
            current_family->set_assigned_day(current_day);
        }
        // move the family to randomly chosen, but really acceptable preferred day
        current_day = current_family->get_random_preferred_day_within_threshold(abort_threshold);
        current_cost_variation += current_family->set_assigned_day(current_day);
        //std::cout << "setting family " << current_family->get_id() << " to day " << current_day->get_id()<<std::endl;
        if (current_day->is_feasible() && current_cost_variation < 0) {
            //A->check_solution_is_ok();
            if (visited_days.size() < count_order_change.size())
                count_order_change[visited_days.size()]++;
            nb_successful_augmenting_path_0++;
            return;
        }
        if (current_cost_variation > abort_threshold){
            ongoing = false;
        }
        visited_families.push_back(current_family);
        // now the current day can exceed the limit of 300 so the next current_family
        // is chosen so that its removal would make the day feasible again
        current_family = current_day->get_random_removable_family();
    }

    Day* best_possible_fallback = current_family->get_best_possible_day();
    // the best possible could be the current assigned day in which case this manoeuver is useless
    if (best_possible_fallback->get_id() != current_day->get_id()){
        current_cost_variation += current_family->set_assigned_day(best_possible_fallback);
        if (current_cost_variation <= 0) {
            if (visited_days.size() + 1 < count_order_change.size())
                count_order_change[visited_days.size() + 1]++;
            nb_successful_augmenting_path_1++;
            return;
        }
        current_family->set_assigned_day(current_day);
    }

    // undo what we did
    for (int i = visited_families.size() - 1; i >= 0; i--) {
        visited_families[i]->set_assigned_day(visited_days[i]);
        //std::cout << "setting family " << visited_families[i]->get_id() << " back to day " << visited_days[i]->get_id() << std::endl;
    }
    //A->check_solution_is_ok();
}

bool is_day_in_cycle(unsigned int family_size, Day* current_day, unsigned int &res,
                     std::vector<bool> &day_was_visited,
                     std::vector<bool> &day_is_on_stack,
                     std::vector<Family *> &predecessor_family){
    if (day_was_visited[current_day->get_id()]) return false;
    day_is_on_stack[current_day->get_id()] = true;

    std::vector<Family*> families_to_visit(0);
    for (unsigned int i = 0; i < current_day->get_nb_families(); i++){
        Family* current_family = current_day->get_ith_family(i);
        if (current_family->get_nb_people() != family_size) continue;
        unsigned int k = current_family->get_k();
        for (unsigned int i = 0; i < k; i++) {
            Day* next_day = current_family->get_ith_preferred_day(i);
            if (day_is_on_stack[next_day->get_id()]){
                predecessor_family[next_day->get_id()] = current_family;
                res = current_day->get_id();
                return true;
            }
            if (day_was_visited[next_day->get_id()]) continue;
            predecessor_family[next_day->get_id()] = current_family;
            if (is_day_in_cycle(family_size, next_day, res, day_was_visited, day_is_on_stack, predecessor_family))
                return true;
        }
    }

    day_is_on_stack[current_day->get_id()] = false;
    day_was_visited[current_day->get_id()] = true;
    return false;
}

void LocalSearch::oriented_cycle_move() {
    unsigned int family_size = 2 + rand() % 7;

    std::vector<bool> day_was_visited(NB_DAYS, false);
    std::vector<bool> day_is_on_stack(NB_DAYS, false);
    std::vector<Family *> predecessor_family(NB_DAYS);
    Day* current_day;
    unsigned int res;

    bool found_unvisited_day, found_cycle = false;
    do {
        predecessor_family = std::vector<Family *>(NB_DAYS);
        // go to next day not yet visited
        found_unvisited_day = false;
        unsigned int j = 0;
        for (; j < NB_DAYS && !found_unvisited_day; j++)
            if (!day_was_visited[j])
                found_unvisited_day = true;
        current_day = A->get_ith_day(j - 1);
        if (is_day_in_cycle(family_size, current_day, res, day_was_visited, day_is_on_stack, predecessor_family))
            found_cycle = true;
    } while (found_unvisited_day && !found_cycle);

    if (!found_cycle) return;

    nb_successful_oriented_cycle++;
    current_day = A->get_ith_day(res);
    unsigned int last_id = current_day->get_id();
    unsigned int nb_swaps = 0;if (current_day->get_id() == predecessor_family[current_day->get_id()]->get_assigned_day()->get_id()) std::cout<<"grd"<<std::endl;
    do{nb_swaps++;
        Family* current_family = predecessor_family[current_day->get_id()];
        Day* tmp = current_family->get_assigned_day();
        current_family->set_assigned_day(current_day);
        current_day = tmp;
    } while(current_day->get_id() != last_id);std::cout<<nb_swaps<<std::endl;
}

void LocalSearch::stats() const {
    std::cout << "  stats: histogram on the number of assignments made by move augmenting_path" << std::endl;
    for (unsigned int i=0; i < count_order_change.size(); i++)
        std::cout << i << "       ";
    std::cout << std::endl;
    for (unsigned int i=0; i < count_order_change.size(); i++){
        std::cout << count_order_change[i] << std::string(7-nb_chiffres(count_order_change[i]) + nb_chiffres(i), ' ');
    }
    std::cout << std::endl << std::endl;
    std::cout << "nb of successful augmenting path moves (no repair): " << nb_successful_augmenting_path_0 << std::endl;
    std::cout << "nb of successful augmenting path moves (with repair): " << nb_successful_augmenting_path_1 << std::endl;
    std::cout << "nb of successful augmenting path moves (cycling): " << nb_successful_augmenting_path_2 << std::endl;
    std::cout << "nb of successful oriented cycle moves: " << nb_successful_oriented_cycle << std::endl;
}

