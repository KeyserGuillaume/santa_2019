#include "BranchAndBound.h"
#include <queue>

std::vector<uint_pair> do_assignments_by_default(Presets& presets){
    std::vector<uint_pair> assignments_by_default = presets.get_assignments_by_default();
    for (unsigned int m = 0; m < assignments_by_default.size(); m++)
        presets.assign_family(assignments_by_default[m].first, assignments_by_default[m].second, false);
    if (assignments_by_default.size() > 0)
        presets.compute_all_bounds(); presets.compute_feasibility();
    return assignments_by_default;
}

void undo_assignments_by_default(Presets& presets, const std::vector<uint_pair> &assignments_by_default){
    for (unsigned int m = 0; m < assignments_by_default.size(); m++) {
        presets.deassign_family(assignments_by_default[m].first, false);
        for (unsigned int k = 0; k < K_MAX; k++)
            if (k != assignments_by_default[m].second)
                presets.forbid_assignment(assignments_by_default[m].first, k, false);
    }
    if (assignments_by_default.size() > 0)
        presets.compute_all_bounds(); presets.compute_feasibility();
}

// uses the lb given by flow. gives optimal solution for 20 families I think
unsigned int brute_force(Presets &presets,
                         unsigned int& nb_nodes,
                         clock_t &time_graph,
                         clock_t &time_bounds,
                         clock_t &time_anticipating,
                         std::vector<unsigned int> reference_costs,
                         const FamilyDistribution &distribution){
    // b&b with lb computed by flow + presets
    unsigned int cost = presets.get_presets_costs() + presets.get_day_cost_lb();
    if (cost > BEST_SOLUTION)
        return cost;
    if (presets.get_nb_assignments() == NB_FAMILIES){
        presets.write_solution("flow_solution");
        return cost;
    }

    // do the default assignments (families which only have one choice left)
    std::vector<uint_pair> assignments_by_default = do_assignments_by_default(presets);
    if (!presets.is_feasible()) {
        undo_assignments_by_default(presets, assignments_by_default);
        return BEST_SOLUTION + 1;
    }

    nb_nodes++;
    clock_t t0 = clock();
    Graph G = Graph(presets);
    G.set_prior_paths(distribution);
    G.compute_max_flow_min_cost();
    time_graph += clock() - t0;
    std::cout << "flow is " << G.get_current_flow() << " with current cost " << G.get_flow_cost()
              << " and with true cost " << G.get_true_flow_cost() << std::endl;
    cost = G.get_true_flow_cost();
    if (!G.is_flow_maximal() || cost > BEST_SOLUTION) {
        // indicate that there is no acceptable solution from this preset
        undo_assignments_by_default(presets, assignments_by_default);
        return BEST_SOLUTION + 1;
    }
    if (G.get_current_flow() == 0) {
        G.check_day_costs_are_ok();
        presets.write_solution("flow_solution");
        undo_assignments_by_default(presets, assignments_by_default);
        return cost;
    }

    if (G.get_most_dispersed_family().second == 1 && G.day_costs_are_ok()) {
        G.check_day_costs_are_ok();
        G.get_solution().write_solution("flow_solution");
        undo_assignments_by_default(presets, assignments_by_default);
        return cost;
    }

    unsigned int current_family = (unsigned int)(presets.get_largest_unassigned_families()[0]);

    FamilyDistribution new_distribution = G.get_distribution();
    preset old_preset = presets.get_preset(current_family);
    unsigned int min_cost = BEST_SOLUTION + 1;
    for (unsigned int k_assign = 0; k_assign < K_MAX; k_assign++){
        if (old_preset[k_assign] == FORBIDDEN)
            continue;

        // we branch on current_family : here we assign it to k_assign.
        t0 = clock();
        presets.assign_family(current_family, k_assign);
        time_bounds += clock() - t0;

        // only do the next parts if it looks feasible
        if (presets.is_feasible()) {
            // anticipation stuff
            t0 = clock();
            std::vector<std::pair<unsigned int, preset>> saves_of_presets(0);
            std::vector<std::vector<unsigned int>> assignations(K_MAX, std::vector<unsigned int>(0));
            std::vector<std::vector<unsigned int>> counter_assignations(K_MAX, std::vector<unsigned int>(0));
            for (unsigned int k = 0; k < K_MAX; k++){
                if (old_preset[k] == FORBIDDEN)
                    continue;
                unsigned int current_day = presets.get_family_data(current_family, k);
                anticipate_on_day(presets, current_day, reference_costs[current_day], assignations[k], counter_assignations[k]);
                for (unsigned int& family_index: assignations[k]){
                    saves_of_presets.push_back(std::pair<unsigned int, preset>(family_index, presets.get_preset(family_index)));
                    presets.assign_family(family_index, current_day, false, false);
                }
                for (unsigned int& family_index: counter_assignations[k]){
                    presets.forbid_assignment(family_index, current_day, false, false);
                }
                if (assignations[k].size() > 0 || counter_assignations[k].size() > 0)
                    presets.compute_all_bounds(); presets.compute_feasibility();
            }
            time_anticipating += clock() - t0;

            // we have assigned current_family, and hopefully we've anticipated some assignations
            // and counter_assignations. Now we call the function to go down in the exploration tree.
            cost = brute_force(presets, nb_nodes, time_graph, time_bounds, time_anticipating, reference_costs, new_distribution);
            if (cost < min_cost) {
                min_cost = cost;
            }

            // now we will remove the anticipated stuff, it is no longer valid when we assign to another day.
            t0 = clock();
            for (unsigned int m = 0; m < saves_of_presets.size(); m++)
                presets.set_preset(saves_of_presets[m].first, saves_of_presets[m].second, false);
            for (unsigned int k = 0; k < K_MAX; k++)
                for (unsigned int s = 0; s < counter_assignations[k].size(); s++)
                    presets.enable_assignment(counter_assignations[k][s], presets.get_family_data(current_family, k), false, false);
            presets.compute_all_bounds(); presets.compute_feasibility();
            time_anticipating += clock() - t0;
        }

    }

    // now we put the preset for current_family back to what it was before
    t0 = clock();
    presets.set_preset(current_family, old_preset);
    time_bounds += clock() - t0;

    // finally we do the same for the assignments by default done at the beginning
    undo_assignments_by_default(presets, assignments_by_default);
    return min_cost;
}

// uses lb of RLLB. used to give optimal solution for up to 1000 families, not sure now since the bound has really improved
float brute_force(Presets &presets,
                         unsigned int& nb_nodes,
                         clock_t &time_lagrangian_lb,
                         clock_t &time_bounds,
                         RLLB& rllb,
                         float &best_cost){
    // b&b with lower bound computed with lagrangian relaxation
    // do the default assignments (families which only have one choice left)
    std::vector<uint_pair> assignments_by_default = do_assignments_by_default(presets);
    float cost = presets.get_presets_costs() + presets.get_day_cost_lb();
    if (!presets.is_feasible() || cost >= best_cost) {
        undo_assignments_by_default(presets, assignments_by_default);
        return best_cost + 1;
    }
    if (presets.get_nb_assignments() == NB_FAMILIES){
        presets.write_solution("lagrangian_bound_solution");
        undo_assignments_by_default(presets, assignments_by_default);
        return cost;
    }

    nb_nodes++;
    clock_t t0 = clock();
    if (!rllb.is_primal_feasible())
        rllb.optimize_lambda(presets, false, best_cost);
    cost = rllb.get_lb();
    // If we are primal feasible there is no need to continue because the solution is optimal.
    if (rllb.is_primal_feasible() && cost < best_cost) {
        write_solution_(rllb.get_solution(presets), "../../solutions/rllb_lb_" + std::to_string(cost));
        best_cost = cost;
        return cost;
    }
    if (cost >= best_cost) {
        // indicate that there is no acceptable solution from this preset
        undo_assignments_by_default(presets, assignments_by_default);
        return best_cost + 1;
    }

    time_lagrangian_lb += clock() - t0;

    //unsigned int branching_family = (unsigned int)(presets.get_largest_unassigned_families()[0]);
    unsigned int branching_family = rllb.suggest_branching_family(presets);

    preset old_preset = presets.get_preset(branching_family);
    unsigned int min_cost = best_cost + 1;
    for (unsigned int k_assign = 0; k_assign < K_MAX; k_assign++){
        if (old_preset[k_assign] == FORBIDDEN)
            continue;

        // we branch on branching_family : here we assign it to k_assign.
        t0 = clock();
        presets.assign_family(branching_family, k_assign);
        time_bounds += clock() - t0;

        // only do the next parts if it looks feasible
        if (presets.is_feasible()) {
            RLLB new_rllb = rllb;
            new_rllb.notify_assignment(presets, branching_family, k_assign);
            if (new_rllb.get_lb() < best_cost) {
                cost = brute_force(presets, nb_nodes, time_lagrangian_lb, time_bounds, new_rllb, best_cost);
                if (cost < min_cost)
                    min_cost = cost;
            }
        }
    }

    // now we put the preset for branching_family back to what it was before
    t0 = clock();
    presets.set_preset(branching_family, old_preset);
    time_bounds += clock() - t0;

    // finally we do the same for the assignments by default done at the beginning
    undo_assignments_by_default(presets, assignments_by_default);
    return min_cost;
}




Node::Node(const Presets &presets, const RLLB &rllb, float* best_cost): presets(presets), rllb(rllb), best_cost(best_cost) {
    level = 0;
    cost = rllb.get_lb();
}

Node Node::get_child_node() {
    Node child_node(presets, rllb, best_cost);
    child_node.level = level + 1;
    child_node.parent_node = this;
    return child_node;
}

void Node::make_children(clock_t &time_lagrangian_lb, clock_t &time_bounds) {
    if (is_closed)
        return;
    child_nodes = std::vector<Node>(0);
    clock_t t0;
    // if this node's lb was not computed we compute it and we close it if the lb is too low.
    // here we turn the algorith into a semi-greedy oriented search descent by optimizing lambda every tenth level
    if ((last_minute_optimization && !no_optimization) || (level%1 == 0 && no_optimization)){
        t0 = clock();
        rllb.optimize_lambda(presets, false, *best_cost);
        time_lagrangian_lb += clock() - t0;
        cost = rllb.get_lb();
        if (cost > *best_cost) {
            close(false);
            return;
        }
        if (rllb.is_primal_feasible()) {
            write_solution_(rllb.get_solution(presets), "../../solutions/rllb_lb_" + std::to_string(cost));
            *best_cost = cost;
            return;
        }
    }
    std::cout << "Expanding a node at level " << level << " and cost " << cost << std::endl;
    t0 = clock();
    branching_family = rllb.suggest_branching_family(presets);
//    branching_family = rllb.suggest_best_branching_family(presets); // really not looking good
    preset old_preset = presets.get_preset(branching_family);
    for (unsigned int k_assign = 0; k_assign < K_MAX; k_assign++){
        if (old_preset[k_assign] == FORBIDDEN)
            continue;

        // we branch on branching_family : here we assign it to k_assign.
        Node child_node = get_child_node();
        t0 = clock();
        child_node.presets.assign_family(branching_family, k_assign);

        // if it is unfeasible, stop this child node
        if (!child_node.presets.is_feasible()) continue;

        child_node.rllb.notify_assignment(child_node.presets, branching_family, k_assign);
        if (!last_minute_optimization && !no_optimization) child_node.rllb.optimize_lambda(child_node.presets, false, *best_cost);

        child_node.cost = child_node.rllb.get_lb();
        if (child_node.cost >= *best_cost) continue;

        if (child_node.rllb.is_primal_feasible()) {
            write_solution_(child_node.rllb.get_solution(child_node.presets), "../../solutions/rllb_lb_" + std::to_string(child_node.cost));
            *best_cost = cost;
            continue;
        }

        child_nodes.push_back(child_node);
        nb_open_child_nodes++;
    }
    if (child_nodes.size() == 0)
        close();
}

void Node::check_branching_on_family(const unsigned int &i) {
    std::cout << "check branching on family " << i << " of node of level " << level << std::endl;
    preset old_preset = presets.get_preset(i);
    bool found = false;
    for (unsigned int k_assign = 0; k_assign < K_MAX && !found; k_assign++) {
        if (old_preset[k_assign] == FORBIDDEN)
            continue;
        // we branch on branching_family : here we assign it to k_assign.
        Node child_node = get_child_node();
//        t0 = clock();
        child_node.presets.assign_family(i, k_assign);
//        time_bounds += clock() - t0;
        // if it is unfeasible, stop this child node
        if (!child_node.presets.is_feasible()) continue;
//        t0 = clock();
        child_node.rllb.notify_assignment(child_node.presets, i, k_assign);
//        if (!last_minute_optimization) child_node.rllb.optimize_lambda(child_node.presets, false, best_cost);
//        time_lagrangian_lb += clock() - t0;
        float child_node_cost = child_node.rllb.get_lb();
        if (child_node_cost >= *best_cost) continue;
        found = true;
    }
    if (!found){
        for (unsigned int m = 0; m < child_nodes.size(); m++)
            child_nodes[m].close(false);
        close(false);
        if (level > 0)
            parent_node->notify_child_closed_from_branching_on_family(i);
    }
}

void Node::close(const bool & from_branching) {
    is_closed = true;
    std::cout << "closing node of level " << level << " and cost " << cost << (from_branching ? " from branching " : " ") << std::endl;
    if (from_branching && level > 0)
        parent_node->notify_child_closed_from_branching_on_family(branching_family);
}

void Node::notify_child_closed_from_branching_on_family(const unsigned int &i) {
    nb_open_child_nodes--;
    if (nb_open_child_nodes == 0)
        close(false);
    check_branching_on_family(i);
}

// An implementation of branch&bound with a semi-greedy method (always explore first the nodes with max depth and min lb)
// a modification (optimize the lambda every 500th / 300th / 100th / ... / 10th level depth) has yielded very good feasible solutions
// has memory problems though... so watch the RAM and be ready to end the process !!!
// the pb is that the nodes that are created are never removed from memory
// because the tree + priority queue structure stinks.
float oriented_branch_and_bound(Presets &presets, unsigned int &nb_nodes, clock_t &time_lagrangian_lb, clock_t &time_bounds,
                                RLLB &rllb, float *best_cost) {
    std::priority_queue<NodePtr> search_tree;
    Node root = Node(presets, rllb, best_cost);
    NodePtr first = NodePtr(&root);
    search_tree.push(first);
    while (!search_tree.empty()){
        first = search_tree.top();
        search_tree.pop();
        nb_nodes++;
        first.node->make_children(time_lagrangian_lb, time_bounds);
        for (unsigned int i = 0; i < first.node->child_nodes.size(); i++)
            search_tree.push(NodePtr(&first.node->child_nodes[i]));
//        // just checking that the queue's order is consistent
//        std::vector<NodePtr> tmp_storage(0);
//        NodePtr nptr;
//        NodePtr prev_nptr = search_tree.top();
//        while (!search_tree.empty()) {
//            nptr = search_tree.top();
//            if (nptr > prev_nptr) throw std::logic_error("kh");
//            prev_nptr = nptr;
//            tmp_storage.push_back(nptr);
//            search_tree.pop();
//        }
//        for (unsigned int m = 0; m < tmp_storage.size(); m++)
//            search_tree.push(tmp_storage[m]);
    }
    return *best_cost;
}

// A function that determines optimal bounds on the day occupancies.
// Could improve the bound (by diverting work from the lambda) but mostly
// accelerates the lb computation.
// a method of rllb sharing the same name gives much faster results, but equivalent to
// the result of this algorithm without re-optimizing on the lambdas.
// (better try both algorithms and analyze the results)
void compute_true_day_occupancy_bounds(Presets& presets, const RLLB& rllb) {
    RLLB tmp_rllb;
    float cutting_bound = 68900;
    float cutting_point = 1/8.;
    std::vector<bool> pruned_from_below(NB_DAYS, true);
    std::vector<bool> pruned_from_above(NB_DAYS, true);
    for (unsigned int s = 0; s < 100; s++){
        for (unsigned int j = 0; j < NB_DAYS; j++){
            bool exceeds_bound;
            presets.compute_all_bounds(false);
            unsigned int a = presets.get_occupancy_lb(j);
            unsigned int b = presets.get_occupancy_ub(j);
            if (pruned_from_below[j]) {
                presets.prescribe_occupancy_ub(j, (unsigned int) (a * (1 - cutting_point) + b * cutting_point));
                presets.compute_all_bounds(false);
                if (presets.is_feasible()) {
                    tmp_rllb = RLLB(presets, rllb);
                    tmp_rllb.optimize_lambda(presets, false, cutting_bound);
                    exceeds_bound = tmp_rllb.get_lb() >= cutting_bound;
                } else {
                    exceeds_bound = true;
                }
                presets.pop_last_occupancy_ub_prescription();
                if (exceeds_bound) {
                    presets.prescribe_occupancy_lb(j, (unsigned int) (a * (1 - cutting_point) + b * cutting_point) + 1);
                    presets.compute_all_bounds(false);
                    tmp_rllb = RLLB(presets, rllb);
                    if (tmp_rllb.get_lb() > cutting_bound) throw std::logic_error("smthg is wrong");
                    std::cout << "day " << j << ":  " << presets.get_occupancy_lb(j) << " <= N <= "
                              << presets.get_occupancy_ub(j) << std::endl;
                } else {
                    pruned_from_below[j] = false;
                }
            }
            if (pruned_from_above[j]) {
                presets.prescribe_occupancy_lb(j, (unsigned int) (a * cutting_point + b * (1 - cutting_point)));
                presets.compute_all_bounds();
                if (presets.is_feasible()) {
                    tmp_rllb = RLLB(presets, rllb);
                    tmp_rllb.optimize_lambda(presets, false, cutting_bound);
                    exceeds_bound = tmp_rllb.get_lb() >= cutting_bound;
                } else {
                    exceeds_bound = true;
                }
                presets.pop_last_occupancy_lb_prescription();
                if (exceeds_bound) {
                    presets.prescribe_occupancy_ub(j, (unsigned int) (a * cutting_point + b * (1 - cutting_point)) - 1);
                    presets.compute_all_bounds();
                    tmp_rllb = RLLB(presets, rllb);
                    if (tmp_rllb.get_lb() > cutting_bound) throw std::logic_error("smthg is wrong");
                    std::cout << "day " << j << ":  " << presets.get_occupancy_lb(j) << " <= N <= "
                              << presets.get_occupancy_ub(j) << std::endl;
                } else {
                    pruned_from_above[j] = false;
                }
            }
        }
        bool found = false;
        for (unsigned int j = 0; j < NB_DAYS; j++) {
            std::cout << "day " << j << ":  " << presets.get_occupancy_lb(j) << " <= N <= " << presets.get_occupancy_ub(j) << std::endl;
            found = found || pruned_from_below[j] || pruned_from_above[j];
        }
        if (!found)
            break;
    }
}

// the bound was good enough (68830 I think for assumed optimal solution 68888) that we could find ~2000 families for whom only
// one of their choices could preserve optimality. This in turn improved the lb to 68848. Pb is that it takes a long time.
void compute_compulsory_assignations(Presets &presets, RLLB& rllb) {
    unsigned int nb_possible_assign;
    unsigned int last_possible_assign;
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        if (presets.is_family_alr_assigned(i)) continue;
        preset old_preset = presets.get_preset(i);
        RLLB specific_rllb;
        nb_possible_assign = 0;
        for (unsigned int k_assign = 0; k_assign < K_MAX && nb_possible_assign <= 1; k_assign++) {
            specific_rllb = rllb;
            if (old_preset[k_assign] == FORBIDDEN)
                continue;
            presets.assign_family(i, k_assign);
            if (!presets.is_feasible()) continue;
            specific_rllb.notify_assignment(presets, i, k_assign);
            if (specific_rllb.get_lb() < 68900) {
                nb_possible_assign++;
                last_possible_assign = k_assign;
            }
        }
        if (nb_possible_assign == 0)
            throw std::logic_error("jhvg");
        else if (nb_possible_assign == 1){
            presets.assign_family(i, last_possible_assign);
            rllb.notify_assignment(presets, i, last_possible_assign);
            std::cout << i << " was assigned " << std::endl;
        }
        else
            presets.set_preset(i, old_preset);
    }
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        if (presets.is_family_alr_assigned(i))
            for (unsigned int k = 0; k < K_MAX; k++)
                if (presets[i][k] == COMPULSORY)
                    std::cout << i << " " << k << std::endl;
}

// the bound was good enough (68848 for assumed optimal solution 68888) that we could reduce possible choices to
// 2 or 3 of them. This made the lb computation quicker
void compute_forbidden_assignations(Presets &presets, RLLB& rllb) {
    std::vector<uint_pair> forbidden_assignations(0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        if (presets.is_family_alr_assigned(i)) continue;
        preset old_preset = presets.get_preset(i);
        RLLB specific_rllb;
        for (unsigned int k_assign = 0; k_assign < K_MAX; k_assign++) {
            specific_rllb = rllb;
            if (old_preset[k_assign] == FORBIDDEN)
                continue;
            presets.assign_family(i, k_assign);
            if (!presets.is_feasible()){
                std::cout << "forbid from assigning family " << i << " to choice " << k_assign << std::endl;
                forbidden_assignations.push_back(uint_pair(i, k_assign));
                continue;
            }
            specific_rllb.notify_assignment(presets, i, k_assign);
            if (specific_rllb.get_lb() > 68900) {
                std::cout << "forbid from assigning family " << i << " to choice " << k_assign << std::endl;
                forbidden_assignations.push_back(uint_pair(i, k_assign));
            }
        }
        presets.set_preset(i, old_preset);
    }
    for (unsigned int m = 0; m < forbidden_assignations.size(); m++)
        presets.forbid_assignment(forbidden_assignations[m].first, forbidden_assignations[m].second);
}
