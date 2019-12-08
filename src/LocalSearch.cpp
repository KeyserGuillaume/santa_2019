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
    if (rand()%3 == 0)
        twin_paths_move();
    else
        augmenting_path_move();
}

void LocalSearch::twin_paths_move(){//std::cout<<"begin"<<std::endl;
    int abort_threshold = 100;
    int current_cost_variation = 0;
    int temp_cost_variation;
    unsigned int c = 0;
    bool A_is_at_first = false;
    bool B_is_at_first = false;
    bool A_is_at_Friday_like = false, B_is_at_Friday_like = false, big_is_at_Friday_like;
    bool abort = false;

    Day* last_day_A, *last_day_B;
    std::vector<Day*> visited_days;
    std::vector<Family*> visited_families;
    Day* current_day, *twin_day;
    Family* current_family;

    Day* first_day = A->get_random_day();
    Family* big_family = first_day->get_random_family();

    if (big_family->get_k() == 0) return;

    // save state
    visited_days.push_back(first_day);
    visited_families.push_back(big_family);
    // choose day
    do {
        current_day = big_family->get_random_improving_day();
    } while (current_day->get_id() == first_day->get_id());
    big_is_at_Friday_like = current_day->is_Friday_like();
    // assign
    current_cost_variation += big_family->set_assigned_day(current_day);

    // Initialize A
    // choose family
    do{
        current_family = current_day->get_random_family();
    } while(current_family->get_id() == big_family->get_id());
    // check whether this family likes first_day
    temp_cost_variation = current_family->set_assigned_day(first_day);
    if (temp_cost_variation < abort_threshold) {
        visited_families.push_back(current_family);
        visited_days.push_back(current_day);
        current_cost_variation += temp_cost_variation;
        last_day_A = first_day;
        A_is_at_first = true;
    }
    else {
        // undo
        current_family->set_assigned_day(current_day);
        // save state
        visited_families.push_back(current_family);
        visited_days.push_back(current_day);
        // choose day
        c = 0;
        do {
            current_day = current_family->get_random_preferred_day_within_threshold(abort_threshold);
            c++;
        } while (c < 10 && (current_day->get_id() == first_day->get_id() ||
                            current_day->get_id() == current_family->get_assigned_day()->get_id()));
        if (c == 10){
            for (int i = visited_families.size() - 1; i >= 0; i--)
                visited_families[i]->set_assigned_day(visited_days[i]);
            //A->check_solution_is_ok();
            return;
        }
        A_is_at_Friday_like = current_day->is_Friday_like();
        // assign
        current_cost_variation += current_family->set_assigned_day(current_day);
        last_day_A = current_day;
    }
    // initialize B
    // check a family can be removed to make the day feasible
    current_day = big_family->get_assigned_day();
    if (!current_day->has_removable_family()) {
        for (int i = visited_families.size() - 1; i >= 0; i--)
            visited_families[i]->set_assigned_day(visited_days[i]);
        //A->check_solution_is_ok();
        return;
    }
    // choose family
    do {
        current_family = current_day->get_random_removable_family(big_is_at_Friday_like);
        c++;
    } while (c < 10 && current_family->get_id() == big_family->get_id());
    if (c == 10){
        for (int i = visited_families.size() - 1; i >= 0; i--)
            visited_families[i]->set_assigned_day(visited_days[i]);
        //A->check_solution_is_ok();
        return;
    }
    // check whether this family likes first day
    temp_cost_variation = current_family->set_assigned_day(first_day);
    if ((A_is_at_first && first_day->is_feasible() && temp_cost_variation + current_cost_variation <= 0) ||
        (!A_is_at_first && temp_cost_variation < abort_threshold)) {
        visited_families.push_back(current_family);
        visited_days.push_back(current_day);
        current_cost_variation += temp_cost_variation;
        B_is_at_first = true;
    }
    else {
        // undo
        current_family->set_assigned_day(current_day);
        // save state
        visited_families.push_back(current_family);
        visited_days.push_back(current_day);
        // choose day
        do {
            current_day = current_family->get_random_preferred_day_within_threshold(abort_threshold);
            c++;
        } while (c < 10 && (current_day->get_id() == first_day->get_id() ||
                            current_day->get_id() == current_family->get_assigned_day()->get_id()));
        if (c == 10) {
            for (int i = visited_families.size() - 1; i >= 0; i--)
                visited_families[i]->set_assigned_day(visited_days[i]);
            //A->check_solution_is_ok();
            return;
        }
        // assign
        B_is_at_Friday_like = current_day->is_Friday_like();
        current_cost_variation += current_family->set_assigned_day(current_day);
        last_day_B = current_day;
    }

    bool do_path_A = false; // we really do start with path_A

    while((current_cost_variation < abort_threshold || A_is_at_Friday_like || B_is_at_Friday_like) &&
          !(last_day_A->is_feasible() && last_day_B->is_feasible() && first_day->is_feasible() && !A_is_at_Friday_like && !B_is_at_Friday_like) &&
          !abort){
        do_path_A = 1 - do_path_A;
        if (do_path_A && (A_is_at_first || (last_day_A->is_feasible() && !last_day_B->is_feasible() && !B_is_at_first))) continue;
        if (!do_path_A && (B_is_at_first || (last_day_B->is_feasible() && !last_day_A->is_feasible() && !A_is_at_first))) continue;
        if (do_path_A){
            current_day = last_day_A;
            twin_day = last_day_B;
        } else {
            current_day = last_day_B;
            twin_day = last_day_A;
        }

        // choose family
        if (twin_day->get_id() == current_day->get_id())
            current_family = current_day->get_random_family();
        else if (current_day->has_removable_family())
           current_family = current_day->get_random_removable_family();
        else{
            abort = true;
            break;
        }
        // check whether this family likes first day
        temp_cost_variation = current_family->set_assigned_day(first_day);
        if (((A_is_at_first || B_is_at_first) && first_day->is_feasible() && temp_cost_variation + current_cost_variation <= 0) ||
            (!A_is_at_first && !B_is_at_first && temp_cost_variation < abort_threshold)) {
            current_cost_variation += temp_cost_variation;
            visited_days.push_back(current_day);
            visited_families.push_back(current_family);
            if (do_path_A) {
                A_is_at_first = true;
                A_is_at_Friday_like = false;
            }
            else {
                B_is_at_first = true;
                B_is_at_Friday_like = false;
            }
            current_day = first_day;
        } else {
            // undo
            current_family->set_assigned_day(current_day);
            // save state
            visited_families.push_back(current_family);
            visited_days.push_back(current_day);
            // choose day
            c = 0;
            do {
                current_day = current_family->get_random_preferred_day_within_threshold(abort_threshold);
                c++;
            } while (c < 10 && (current_day->get_id() == first_day->get_id() ||
                                current_day->get_id() == current_family->get_assigned_day()->get_id()));
            if (do_path_A)
                A_is_at_Friday_like = current_day->is_Friday_like();
            else
                B_is_at_Friday_like = current_day->is_Friday_like();
            if (c == 10)
                abort = true;
            else
                current_cost_variation += current_family->set_assigned_day(current_day);
        }
        if (do_path_A)
            last_day_A = current_day;
        else
            last_day_B = current_day;
    }

    if (!first_day->is_feasible() || current_cost_variation > 0 || abort){
        // undo what we did
        //std::cout << "undo" << std::endl;
        for (int i = visited_families.size() - 1; i >= 0; i--)
            visited_families[i]->set_assigned_day(visited_days[i]);
        return;
    }

    count_size_moved_by_twin_paths[big_family->get_nb_people()]++;
    nb_successful_twin_paths++;
    //A->check_solution_is_ok();
}

void LocalSearch::augmenting_path_move() {
    int abort_threshold = 100;

    std::vector<Day*> visited_days(0);
    std::vector<Family*> visited_families(0);
    bool ongoing = true;
    bool is_Friday_like = false;
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
        // move the family to randomly chosen preferred day which they could accept
        current_day = current_family->get_random_preferred_day_within_threshold(abort_threshold);
        is_Friday_like = current_day->is_Friday_like();
        current_cost_variation += current_family->set_assigned_day(current_day);
        if (current_day->is_feasible() && current_cost_variation < 0) {
            //A->check_solution_is_ok();
            if (visited_days.size() < count_order_change.size())
                count_order_change[visited_days.size()]++;
            nb_successful_augmenting_path_0++;
            return;
        }
        if (current_cost_variation > abort_threshold && !is_Friday_like){
            ongoing = false;
        }
        visited_families.push_back(current_family);
        // now the current day can exceed the limit of 300 so the next current_family
        // is chosen so that its removal would make the day feasible again
        current_family = current_day->get_random_removable_family(is_Friday_like);
    }

    Day* best_possible_fallback = current_family->get_best_possible_day();
    // the best possible could be the current assigned day in which case this manoeuver is useless
    if ((visited_days.size() > 1 || best_possible_fallback->get_id() != first_day->get_id()) && best_possible_fallback->get_id() != current_day->get_id()){
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
    std::cout << "  stats: histogram on the size of families moved by twin paths" << std::endl;
    for (unsigned int i=2; i < count_size_moved_by_twin_paths.size(); i++)
        std::cout << i << "       ";
    std::cout << std::endl;
    for (unsigned int i=2; i < count_size_moved_by_twin_paths.size(); i++){
        std::cout << count_size_moved_by_twin_paths[i] << std::string(7-nb_chiffres(count_size_moved_by_twin_paths[i]) + nb_chiffres(i), ' ');
    }
    std::cout << std::endl << std::endl;
    std::cout << "nb of successful augmenting path moves (no repair): " << nb_successful_augmenting_path_0 << std::endl;
    std::cout << "nb of successful augmenting path moves (with repair): " << nb_successful_augmenting_path_1 << std::endl;
    std::cout << "nb of successful augmenting path moves (cycling): " << nb_successful_augmenting_path_2 << std::endl;
    std::cout << "nb of successful twin paths moves: " << nb_successful_twin_paths << std::endl;
    std::cout << "nb of successful oriented cycle moves: " << nb_successful_oriented_cycle << std::endl;
}

