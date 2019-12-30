#pragma once

#include <iostream>
#include <vector>
#include <queue>
#include "constants.h"
#include <stdexcept>

#include "tools.h"
#include "Presets.h"

// Vertex must know about Arc and vice-versa
class Arc;


class Vertex {
private:
    unsigned int id;
    std::vector<Arc*> delta_p;
    std::vector<Arc*> delta_m;
public:
    Vertex(){}
    Vertex(unsigned int id): id(id){}
    void add_incoming_arc(Arc* a){
        delta_m.push_back(a);
    }
    void add_outgoing_arc(Arc* a){
        delta_p.push_back(a);
    }
    inline unsigned int get_id() const{return id;}
    unsigned int get_delta_p_size() const{return delta_p.size();}
    unsigned int get_delta_m_size() const{return delta_m.size();}
    Arc* get_ith_incoming(const unsigned int& i){return delta_m[i];}
    Arc* get_ith_outgoing(const unsigned int& i){return delta_p[i];}
};

class Arc{
private:
    unsigned int id;
    Vertex *u, *v;
    unsigned int flow = 0;
    unsigned int capa;
    int cost;
public:
    Arc(){}
    Arc(unsigned int id, Vertex* u, Vertex* v, unsigned int capa, int cost): id(id), u(u), v(v), capa(capa), cost(cost){}//if (capa < 1){throw std::logic_error("shit");}}
    void fill_vertex_deltas(){
        v->add_incoming_arc(this);
        u->add_outgoing_arc(this);
    }
    bool has_flow(){return flow > 0;}
    bool has_capacity_left(){return flow < capa;}
    Vertex* get_u() const {return u;}
    Vertex* get_v() const {return v;}
    unsigned int get_capa() const{return capa;}
    int get_cost() const{return cost;}
    unsigned int get_id()const{return id;}
    unsigned int get_flow()const{return flow;}
    void add_to_flow(int f){flow += f;}
    void clear_flow(){flow = 0;}
};


class Graph{
    // not a general class to encode a graph, just something for the graph I need
private:
    Vertex* V;
    Arc* A;
    unsigned int n_V, n_A;
    unsigned int flow = 0;
    unsigned int max_possible_flow;
    std::vector<int> distances;
    std::vector<Arc*> predecessor;
    std::vector<unsigned int> family_indexes, day_indexes;
    std::vector<int> full_family_indexes;
    unsigned int source_index, sink_index;
    int last_path_value;
    const Presets &presets;
    //std::vector<PriorPath> prior_paths;
    std::priority_queue<PriorPath, std::vector<PriorPath>, PriorPath_compare> prior_paths;
    int path_index_to_try = -1;

public:
    bool debug = false;
    Graph(const Presets &presets, const bool &toy);
    Graph(const Presets &presets);
    ~Graph(){
        delete [] V;
        delete [] A;
    }
    unsigned int get_const_cost(const unsigned int &i, const unsigned int &k) const;
    unsigned int get_marg_cost(const unsigned int &i, const unsigned int &k) const;

    bool find_and_apply_augmenting_path();
    void add_flow_to_augmenting_path(std::vector<Arc*> path, const unsigned int& additional_flow);
    void compute_max_flow_min_cost();
    std::vector<Arc*> get_shortest_path();
    void add_all_prior_paths();
    unsigned int add_obvious_flow(
            const unsigned int &family_index,
            const std::vector<unsigned int> &turns,
            const unsigned int &flow_quantity = 0);
    void update_distances();
    void apply_Bellman_Ford(std::queue<Vertex*> &Q);
    PriorPath process_prior_path(const PriorPath& p);
    void inspect_distances() const;
    unsigned int get_ith_day_flow(unsigned int i) const;
    unsigned int get_ith_day_capa(unsigned int i) const;

    void clear_flow();
    void add_flow_for_assigning_family(const unsigned int &family_idx, const unsigned int &k);
    void set_prior_paths(const FamilyDistribution &my_distribution);
    FamilyDistribution get_distribution();
    Presets get_solution();

    unsigned int get_current_flow() const;
    int get_flow_cost() const;
    int get_true_flow_cost() const;
    void init_distances_and_predecessors();
    void check_flow(const bool &check_maximal_flow = true);
    void check_day_costs_are_ok() const;
    bool day_costs_are_ok(const bool & throw_error = false) const;
    bool is_flow_maximal(const bool &throw_error = false) const;
    void show_distances() const;
    void show_schedule() const;
    void show_dispersion() const;
    std::vector<unsigned int> get_family_dispersion() const;
    std::vector<unsigned int> get_day_occupancy() const;
    uint_pair get_most_dispersed_family() const;
    unsigned int get_largest_least_dispersed_family() const;
    std::vector<float> get_real_day_costs() const;
    int get_overload_family() const;
};

