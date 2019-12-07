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
    augmenting_path_move();
}

void LocalSearch::augmenting_path_move() {
    int abort_threshold = 100;

    std::vector<Day*> visited_days(0);
    std::vector<Family*> visited_families(0);
    bool ongoing = true;
    int current_cost_variation = 0;

    //Day* current_day = A->get_ith_day(0);
    Day* current_day;
    do {
        current_day = A->get_random_day();
    } while(!current_day->has_removable_family());
    Family* current_family = current_day->get_random_removable_family();

    while (ongoing){
        visited_days.push_back(current_day);
        current_day = current_family->get_random_preferred_day();
        current_cost_variation += current_family->set_assigned_day(current_day);
        //std::cout << "setting family " << current_family->get_id() << " to day " << current_day->get_id()<<std::endl;
        if (current_day->is_feasible() && current_cost_variation < 0) return;//{A->check_solution_is_ok();return;}
        if (current_cost_variation > abort_threshold){
            ongoing = false;
        }
        visited_families.push_back(current_family);
        current_family = current_day->get_random_removable_family();
    }

    //current_family.

    // undo what we did
    for (int i = visited_families.size() - 1; i >= 0; i--) {
        visited_families[i]->set_assigned_day(visited_days[i]);
        //std::cout << "setting family " << visited_families[i]->get_id() << " back to day " << visited_days[i]->get_id() << std::endl;
    }
    //A->check_solution_is_ok();
}

