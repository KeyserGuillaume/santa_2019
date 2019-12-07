#pragma once

#include <stdexcept>
#include <vector>

#include "constants.h"
#include "tools.h"

class Day;

class Family{
private:
    unsigned int n_people; // nb people in family
    unsigned int id;
    std::vector<Day*> preferred_days;
    unsigned int k; // the family got their k-th choice (begins at 0, 10 if not one of their choices)
    unsigned int cost;
    Day* assigned_day;
    void compute_cost();
public:
    Family(){}
    Family(const unsigned int &id, const unsigned int &n_people, const std::vector<Day*> &preferred_days, Day* assigned_day);
    unsigned int set_assigned_day(Day* i);
    Day* get_assigned_day() const{return assigned_day;}
    unsigned int get_nb_people() const{return n_people;}
    unsigned int get_id() const{return id;}
    unsigned int get_cost() const{return cost;}
    Day* get_ith_preferred_day(const unsigned int &i) const{return preferred_days[i];}
    Day* get_random_preferred_day() const{return preferred_days[rand()%preferred_days.size()];}
    unsigned int get_k()const{return k;}
};

class Day{
    unsigned int id;
    unsigned int N = 0;
    unsigned int cost = 0; // no meaning as long as there are less than 125 people
    Day* previous_day, *next_day;
    void compute_cost();
public:
    std::vector<Family*> assigned_families;
    Day(){next_day = this; previous_day = this; assigned_families.clear();}
    void set_id(unsigned int i){id = i;}
    void set_previous_day(Day* d){previous_day = d;}
    void set_next_day(Day* d){next_day = d;}
    unsigned int add_family(Family* family);
    bool has_removable_family() const;
    Family* get_random_removable_family();
    unsigned int remove_family(Family* f);
    unsigned int get_N()const{return N;}
    unsigned int get_cost() const{return cost;}
    unsigned int get_id(){return id;}
    bool is_feasible(){return N >= MIN_NB_PEOPLE_PER_DAY && N <= MAX_NB_PEOPLE_PER_DAY;}
};

class Assignment {
private:
    Family* families;
    Day* days;
public:
    Assignment(const std::vector<std::vector<unsigned int>> &family_data, const std::vector<unsigned int> &solution);
    void write_solution(const std::string &filename) const;
    void stats() const;
    void check_solution_is_ok();
    Day* get_random_day(){return days + rand()%NB_DAYS;}
    unsigned int get_cost()const;
    Day* get_ith_day(const unsigned int &i) const{return days + i;}
    Family* get_ith_family(const unsigned int &i) const{return families + i;}
    ~Assignment(){delete[] families; delete[] days;}
};

