#pragma once

#include "tools.h"
#include "Flow.h"
#include "LocalSearch.h"
#include "Presets.h"
#include "day_sub_problem.h"

unsigned int brute_force(Presets &presets,
                         unsigned int& nb_nodes,
                         clock_t &time_graph,
                         clock_t &time_bounds,
                         clock_t &time_anticipating,
                         std::vector<unsigned int> reference_costs,
                         const FamilyDistribution &distribution);

float brute_force(Presets &presets,
                  unsigned int& nb_nodes,
                  clock_t &time_lagrangian_lb,
                  clock_t &time_bounds,
                  RLLB& rllb,
                  float &best_cost);


class Node{
private:
    Presets presets;
    float* best_cost;
    RLLB rllb;
    unsigned int level;
    float cost;
    Node* parent_node;
    int branching_family = -1;
    unsigned int nb_open_child_nodes = 0;
    bool is_closed = false;
    // below is a parameter
    bool last_minute_optimization = false;
    bool no_optimization = true;
public:
    std::vector<Node> child_nodes;
    float get_cost() const { return cost; }
    unsigned int get_level() const { return level; }
    void make_children(clock_t &time_lagrangian_lb, clock_t &time_bounds);
    void close(const bool &from_branching = true);
    void notify_child_closed_from_branching_on_family(const unsigned int &i);
    void check_branching_on_family(const unsigned int& i);
    bool operator< (const Node& other) const { return (level < other.level || (level == other.level && cost > other.cost)); }
    bool operator> (const Node& other) const { return (level > other.level || (level == other.level && cost < other.cost)); }
    Node () {}
    //Node (const Node& other): presets(other.presets), cost(other.cost), level(other.level), rllb(other.rllb){}
    Node (const Presets &presets, const RLLB &rllb, float* best_cost);
    Node get_child_node();
    friend class NodePtr;
    friend void carry_out_tests(Presets &presets, RLLB& rllb);
};

class NodePtr{
public:
    Node* node;
    NodePtr(Node* node):node(node) {}
    bool operator< (const NodePtr& other) const { return (node->level < other.node->level || (node->level == other.node->level && node->cost > other.node->cost)); }
    bool operator> (const NodePtr& other) const { return (node->level > other.node->level || (node->level == other.node->level && node->cost < other.node->cost)); }
    NodePtr() {}
};

float oriented_branch_and_bound(Presets &presets,
                                unsigned int& nb_nodes,
                                clock_t &time_lagrangian_lb,
                                clock_t &time_bounds,
                                RLLB& rllb,
                                float *best_cost);

void compute_true_day_occupancy_bounds(Presets& presets, const RLLB& rllb);
void compute_compulsory_assignations(Presets& presets, RLLB& rllb);
void compute_forbidden_assignations(Presets &presets, RLLB& rllb);
