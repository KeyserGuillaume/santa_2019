#pragma once

#include "constants.h"
#include "tools.h"
#include "Presets.h"

#include <ilcplex/ilocplex.h>

class Column{
    char* selected_choices; // values of k of all families
    unsigned int cum_nb_choices[NB_FAMILIES + 1]; // family i has its choices between indexes cum_nb_choices[i] and cum_nb_choices[i + 1]
    double cost;
    char* nb_copies;
public:
    double get_cost() const {return std::max(-1000000., cost);}
    unsigned int get_family_use_nb(const unsigned int& i) const {return cum_nb_choices[i + 1] - cum_nb_choices[i];}
    IloNumArray get_family_use_nb(IloEnv &env) const {
        IloNumArray res(env, NB_FAMILIES);
        for (unsigned int i = 0; i < NB_FAMILIES; i++)
            res[i] = get_family_use_nb(i);
        return res;
    }
    unsigned int get_sum_family_use_nb() const {
        unsigned int res = 0;
        for (unsigned int i = 0; i < NB_FAMILIES; i++)
            res += get_family_use_nb(i);
        return res;
    }
    Column(){
        selected_choices = new char[0];
        nb_copies = new char[1];
        *nb_copies = 1;
    }
    // constructor for a column from a general solution (provided in terms of choices, not days)
    Column(const std::vector<std::vector<char>> &solution, const double& cost) : cost(cost){
        unsigned int nb_choices = 0;
        for (unsigned int i = 0; i < NB_FAMILIES; i++) {
            cum_nb_choices[i] = nb_choices;
            nb_choices += solution[i].size();
        }
        cum_nb_choices[NB_FAMILIES] = nb_choices;
        std::cout << "c = " << nb_choices << std::endl;
        selected_choices = new char[nb_choices];
        unsigned int c = 0;
        for (unsigned int i = 0; i < NB_FAMILIES; i++){
            for (unsigned int l = 0; l < solution[i].size(); l++){
                selected_choices[c] = solution[i][l];
                c++;
            }
        }
        nb_copies = new char[1];
        *nb_copies = 1;
    }
    // constructor for the column where all families have zero choices save family i which chose its favourite day
    Column(const unsigned int &i){
        selected_choices = new char[1];
        selected_choices[0] = 0;
        for (unsigned int i0 = 0; i0 < i + 1; i0++)
            cum_nb_choices[i0] = 0;
        for (unsigned int i0 = i + 1; i0 < NB_FAMILIES + 1; i0++)
            cum_nb_choices[i0] = 1;
        // give it a relatively high cost
        cost = LARGE_UPPER_BOUND;
        nb_copies = new char[1];
        *nb_copies = 1;
    }
    // had to define destructor to manage heap memory then even redefine the copy constructor to deal with vector reallocation...
    ~Column(){
        if (nb_copies == 0)
            delete[] selected_choices;
        nb_copies--;
    }
    Column(const Column& c){
        selected_choices = c.selected_choices;
        for (unsigned int i = 0; i < NB_FAMILIES + 1; i++)
            cum_nb_choices[i] = c.cum_nb_choices[i];
        cost = c.cost;
        nb_copies = c.nb_copies;
        *nb_copies ++;
    }
};

void test_stuff(Presets& presets, const std::vector<unsigned int> &initial_solution, const double& cost);
void tiny_test();