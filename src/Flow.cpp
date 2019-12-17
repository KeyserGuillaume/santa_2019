#include "Flow.h"

#include "tools.h"



Graph::Graph() {
    n_V = 5;
    V = new Vertex [n_V];
    n_A = 5;
    A = new Arc [n_A];
    unsigned int next_vertex_index = 0,
            next_arc_index = 0;
    V[next_vertex_index] = Vertex(next_vertex_index);
    next_vertex_index++;
    V[next_vertex_index] = Vertex(next_vertex_index);
    next_vertex_index++;
    V[next_vertex_index] = Vertex(next_vertex_index);
    next_vertex_index++;
    V[next_vertex_index] = Vertex(next_vertex_index);
    next_vertex_index++;
    V[next_vertex_index] = Vertex(next_vertex_index);
    A[0] = Arc(0, V + 0, V + 1, 2, 1);
    A[1] = Arc(1, V + 0, V + 2, 2, 2);
    A[2] = Arc(2, V + 1, V + 3, 2, 1);
    A[3] = Arc(3, V + 2, V + 3, 2, -20);
    A[4] = Arc(4, V + 3, V + 4, 3, 0);
    source_index = 0;
    sink_index = 4;
    // build the deltas
    for (unsigned int i = 0; i < n_A; i++){
        A[i].fill_vertex_deltas();
    }
    init_distances_and_predecessors();
    max_possible_flow = 3;
}

int Graph::init_distances_and_predecessors() {
    distances = std::vector<int> (n_V, 10*UPPER_BOUND);
    predecessor = std::vector<Arc*>(n_V);
    distances[source_index] = 0;
    std::queue<Vertex*> Q;
    Q.push(V + source_index);
    // I assume this is not an infinite loop (no absorbing circuit)
    while (!Q.empty()) {
        apply_Bellman_Ford(Q);
    }
}

Graph::Graph(const std::vector<std::vector<unsigned int>> &family_data, std::vector<preset> presets) {
    // build the graph from right to left
    unsigned int k_max = 4;

    unsigned int nb_already_assigned = 0;
    unsigned int nb_forbidden_assignments = 0;
    std::vector<bool> is_already_assigned (NB_FAMILIES, false);
    std::vector<unsigned int> presets_schedule(NB_DAYS, 0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        for (unsigned int k = 0; k < k_max; k++) {
            if (presets[i][k] == COMPULSORY) {
                presets_schedule[family_data[i][k]] += family_data[i][NB_CHOICES];
                presets_costs += CONSTANT_COST[k] + family_data[i][NB_CHOICES] * MARGINAL_COST[k];
                is_already_assigned[i] = true;
                nb_already_assigned++;
            }
            if (presets[i][k] != ALLOWED)
                nb_forbidden_assignments++;
        }
    }
    n_V = 1 + NB_FAMILIES - nb_already_assigned + NB_FAMILIES*k_max - nb_forbidden_assignments + 2*NB_DAYS + 1;
    V = new Vertex [n_V];
    n_A = NB_FAMILIES - nb_already_assigned + 3*NB_FAMILIES*k_max - 3*nb_forbidden_assignments + NB_DAYS*3;
    A = new Arc [n_A];
    unsigned int next_vertex_index = 0, next_arc_index = 0;
    family_indexes.clear(); day_indexes.clear();

    // the sink
    V[next_vertex_index] = Vertex(next_vertex_index);
    sink_index = next_vertex_index;
    next_vertex_index++;

    // the days
    for (unsigned int i = 0; i < NB_DAYS; i++){
        if (presets_schedule[i] > 300) throw std::logic_error("Why do I have too many people assigned to this day ?");
        unsigned int min_allowed = (unsigned int)(std::max(125 - int(presets_schedule[i]), 0));
        unsigned int max_allowed = MAX_NB_PEOPLE_PER_DAY - presets_schedule[i];
        V[next_vertex_index] = Vertex(next_vertex_index);
        A[next_arc_index] = Arc(
                next_arc_index,
                V + next_vertex_index,
                V,
                max_allowed - min_allowed,
//                MAX_NB_PEOPLE_PER_DAY - MIN_NB_PEOPLE_PER_DAY,
                0);
        next_vertex_index++;
        next_arc_index++;
        V[next_vertex_index] = Vertex(next_vertex_index);
        A[next_arc_index] = Arc(next_arc_index, V + next_vertex_index, V, min_allowed, -UPPER_BOUND);
//        A[next_arc_index] = Arc(next_arc_index, V + next_vertex_index, V, MIN_NB_PEOPLE_PER_DAY, -UPPER_BOUND);
//        A[next_arc_index] = Arc(next_arc_index, V + next_vertex_index, V, schedule[i] - epsilon, -UPPER_BOUND);
        next_arc_index++;
        A[next_arc_index] = Arc(
                next_arc_index,
                V + next_vertex_index,
                V + next_vertex_index - 1,
//                2*epsilon,
 //               MAX_NB_PEOPLE_PER_DAY - MIN_NB_PEOPLE_PER_DAY,
                max_allowed - min_allowed,
                0);
        next_arc_index++;
        day_indexes.push_back(next_vertex_index);
        next_vertex_index++;
    }

    // the families
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        if (is_already_assigned[i]) continue;
        V[next_vertex_index] = Vertex(next_vertex_index);
        family_indexes.push_back(next_vertex_index);
        unsigned int family_index = next_vertex_index;
        next_vertex_index++;
        for (unsigned int k = 0; k < k_max; k++){
            if (presets[i][k] != ALLOWED) continue;
            unsigned int j = family_data[i][k];
            V[next_vertex_index] = Vertex(next_vertex_index);
            A[next_arc_index] = Arc(next_arc_index, V + next_vertex_index, V + day_indexes[j], family_data[i][NB_CHOICES] - 1, 0);
            next_arc_index++;
            A[next_arc_index] = Arc(next_arc_index, V + family_index, V + next_vertex_index, family_data[i][NB_CHOICES] - 1, MARGINAL_COST[k]);
            next_arc_index++;
            A[next_arc_index] = Arc(next_arc_index, V + family_index, V + day_indexes[j], 1, CONSTANT_COST[k]);
            next_arc_index++;
            next_vertex_index++;
        }
    }

    // the source
    V[next_vertex_index] = Vertex(next_vertex_index);
    source_index = next_vertex_index;
    for (unsigned int i = 0; i < family_indexes.size(); i++){
        A[next_arc_index] = Arc(next_arc_index, V + next_vertex_index, V + family_indexes[i], family_data[i][NB_CHOICES], 0);
        next_arc_index++;
    }

    std::cout << next_arc_index << " " << n_A << " " << next_vertex_index << " " << n_V << std::endl;

    // build the deltas
    for (unsigned int i = 0; i < n_A; i++){
        A[i].fill_vertex_deltas();
    }
    init_distances_and_predecessors();

    max_possible_flow = 0;
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        if (!is_already_assigned[i])
            max_possible_flow += family_data[i][NB_CHOICES];
}

//std::vector<unsigned int> Graph::get_solution() const {
//    std::vector<unsigned int> solution (NB_FAMILIES, 101);
//    for (unsigned int i = 0; i < NB_FAMILIES; i++){
//        Vertex* u = V + family_indexes[i];
//        for (unsigned int j = 0; j < u->get_delta_p_size(); j++){
//            Arc* a = u->get_ith_outgoing(j);
//            if (a->has_flow()){
//                unsigned int day_id = a->get_v()->get_id();
//                for (unsigned int k = 0; k < day_indexes.size(); k++){
//                    if (day_indexes[k] == day_id){
//                        solution[i] = k;
//                        break;
//                    }
//                }
//                break;
//            }
//        }
//    }
//    return solution;
//}

bool Graph::find_and_apply_augmenting_path() {
    std::vector<Arc*> path = get_shortest_path();
    if (path.size() == 0) return false;

    unsigned int additional_flow = 10000; // big
    unsigned int current_index = source_index;
    for (unsigned int i = 0; i < path.size(); i++){
        Arc* a = path[i];
        if (a->get_v()->get_id() == current_index){
            current_index = a->get_u()->get_id();
            additional_flow = std::min(additional_flow, a->get_flow());
        }
        else {
            current_index = a->get_v()->get_id();
            additional_flow = std::min(additional_flow, a->get_capa() - a->get_flow());
        }
    }

    current_index = source_index;
    for (unsigned int i = 0; i < path.size(); i++){
        Arc* a = path[i];
        if (a->get_v()->get_id() == current_index){
            current_index = a->get_u()->get_id();
            a->add_to_flow(-additional_flow);
        }
        else {
            current_index = a->get_v()->get_id();
            a->add_to_flow(additional_flow);
        }
    }

    flow += additional_flow;

    return true;
}

void Graph::compute_max_flow_min_cost() {
    unsigned int i = 0;
    while(find_and_apply_augmenting_path()){
        i++;
        if (i%30 == 0)
            std::cout << "flow is " << get_current_flow() << " with cost " << get_flow_cost() << " and paths have length " << get_shortest_path().size() << std::endl;
    }
}

unsigned int Graph::get_current_flow() const {
    unsigned int f1 = 0;
    unsigned int f2 = 0;
    Vertex* source = V + source_index;
    Vertex* sink = V + sink_index;
    unsigned int d1 = source->get_delta_p_size();
    unsigned int d2 = sink->get_delta_m_size();
    for (unsigned int i = 0; i < d1; i++){
        f1 += source->get_ith_outgoing(i)->get_flow();
    }
    for (unsigned int i = 0; i < d2; i++){
        f2 += sink->get_ith_incoming(i)->get_flow();
    }
    if (f1 != f2) throw std::logic_error("Why is my flow inconsistent ?");
    return f1;
}

void Graph::update_distances() {
    std::queue<Vertex *> Q;
    std::vector<Vertex*> reset_vertexes (0);
    std::vector<bool> inside_queue(n_V, false);
    Arc *a;
    Vertex *u, *v, *w;
    unsigned int size_d_plus;
    unsigned int size_d_minus;
    for (unsigned int i = 0; i < n_V; i++){
        if (i == source_index) continue;
        if (distances[i] == 10 * UPPER_BOUND) continue;
        w = V + i;
        a = predecessor[i];
        if ((a->get_v()->get_id() == i && !a->has_capacity_left()) || (a->get_u()->get_id() == i && !a->has_flow())) {
            Q.push(w);
        }
    }
    while(!Q.empty()){
        w = Q.front();
        Q.pop();
        unsigned int i = w->get_id();
        inside_queue[i] = false;
        if (i == source_index) continue;
        if (distances[i] == 10 * UPPER_BOUND) continue;
        a = predecessor[i];
        if ((a->get_v()->get_id() == i && !a->has_capacity_left()) ||
            (a->get_u()->get_id() == i && !a->has_flow()) ||
            (distances[a->get_u()->get_id()] == 10 * UPPER_BOUND) ||
            (distances[a->get_v()->get_id()] == 10 * UPPER_BOUND)) {
            distances[i] = 10 * UPPER_BOUND;
            reset_vertexes.push_back(w);
            size_d_plus = w->get_delta_p_size();
            for (unsigned int j = 0; j < size_d_plus; j++) {
                v = w->get_ith_outgoing(j)->get_v();
                if (!inside_queue[v->get_id()]) {
                    Q.push(v);
                    inside_queue[v->get_id()] = true;
                }
            }
            size_d_minus = w->get_delta_m_size();
            for (unsigned int j = 0; j < size_d_minus; j++) {
                u = w->get_ith_incoming(j)->get_u();
                if (!inside_queue[u->get_id()]) {
                    Q.push(u);
                    inside_queue[u->get_id()] = true;
                }
            }
        }
    }

    // print percentage of distances that need to be recomputed
    //if (rand()%100 == 0) std::cout << 100*reset_vertexes.size()/float(n_V)<<std::endl;

    for (unsigned int i = 0; i < reset_vertexes.size(); i++){
        w = reset_vertexes[i];
        size_d_plus = w->get_delta_p_size();
        for (unsigned int j = 0; j < size_d_plus; j++) {
            a = w->get_ith_outgoing(j);
            v = a->get_v();
            if (!inside_queue[v->get_id()] && distances[v->get_id()] != UPPER_BOUND && a->has_flow()) {
                Q.push(v);
                inside_queue[v->get_id()] = true;
            }
        }
        size_d_minus = w->get_delta_m_size();
        for (unsigned int j = 0; j < size_d_minus; j++) {
            a = w->get_ith_incoming(j);
            u = a->get_u();
            if (!inside_queue[u->get_id()] && distances[u->get_id()] != UPPER_BOUND && a->has_capacity_left()) {
                Q.push(u);
                inside_queue[u->get_id()] = true;
            }
        }
    }

    // I assume this is not an infinite loop (no absorbing circuit)
    while (!Q.empty()) {
        apply_Bellman_Ford(Q);
    }
}


std::vector<Arc*> Graph::get_shortest_path(){
    if (flow < 0.8*max_possible_flow) // not a good strategy
        update_distances();
    else
        init_distances_and_predecessors();

    Vertex *u, *v, *w; Arc* a;
    unsigned int size_d_plus;
    unsigned int size_d_minus;

    if (distances[sink_index] == 10*UPPER_BOUND) return std::vector<Arc*> (0);

    std::vector<Arc*> shortest_path_inverse(0);
    w = V + sink_index;
    while (w->get_id() != source_index){
        a = predecessor[w->get_id()];
        if (a->get_u()->get_id() == w->get_id())
            w = a->get_v();
        else
            w = a->get_u();
        shortest_path_inverse.push_back(a);
    }
    unsigned int path_length = shortest_path_inverse.size();
    std::vector<Arc*> shortest_path(0);
    for (unsigned int i = 0; i < path_length; i++){
        shortest_path.push_back(shortest_path_inverse[path_length - i  - 1]);
    }
    return shortest_path;
}


void Graph::apply_Bellman_Ford(std::queue<Vertex*> &Q) {
    Arc *a;
    Vertex *u, *v, *w;
    unsigned int size_d_plus;
    unsigned int size_d_minus;
    w = Q.front();
    Q.pop();
    if (distances[w->get_id()] == 10 * UPPER_BOUND) return;
    size_d_plus = w->get_delta_p_size();
    for (unsigned int i = 0; i < size_d_plus; i++) {
        a = w->get_ith_outgoing(i);
        v = a->get_v();
        if (a->has_capacity_left() && distances[w->get_id()] + a->get_cost() < distances[v->get_id()]) {
            distances[v->get_id()] = distances[w->get_id()] + a->get_cost();
            Q.push(v);
            predecessor[v->get_id()] = a;
        }
    }
    size_d_minus = w->get_delta_m_size();
    for (unsigned int i = 0; i < size_d_minus; i++) {
        a = w->get_ith_incoming(i);
        u = a->get_u();
        if (a->has_flow() && distances[w->get_id()] - a->get_cost() < distances[u->get_id()]) {
            distances[u->get_id()] = distances[w->get_id()] - a->get_cost();
            Q.push(u);
            predecessor[u->get_id()] = a;
        }
    }
}

int Graph::get_flow_cost() const {
    int res = 0;
    for (unsigned int i = 0; i < n_A; i++)
        res += (A + i)->get_flow() * (A + i)->get_cost();
    return res;
}

int Graph::get_true_flow_cost() const {
    unsigned int res = 0;
    for (unsigned int i = 0; i < family_indexes.size(); i++){
        Vertex* u = V + family_indexes[i];
        unsigned int d_out = u->get_delta_p_size();
        for (unsigned int k = 0; k < d_out; k++){
            Arc* a = u->get_ith_outgoing(k);
            res += a->get_flow() * a->get_cost();
        }
    }
    //std::cout << presets_costs << std::endl;
    return res + presets_costs;
}

void Graph::show_distances() {
    for (unsigned int i = 0; i < n_V; i++){
        std::cout << i << " " << distances[i] << std::endl;
    }
}

void Graph::show_schedule() {
    std::cout << "flows:"<<std::endl;
    for (unsigned int i = 0; i < NB_DAYS; i++) {
        std::cout << (V + day_indexes[i])->get_ith_outgoing(0)->get_flow() + (V + day_indexes[i])->get_ith_outgoing(1)->get_flow() << "  ";
    }
    std::cout << std::endl << std::endl << "capas:"<<std::endl;
    for (unsigned int i = 0; i < NB_DAYS; i++) {
        std::cout << (V + day_indexes[i])->get_ith_outgoing(0)->get_capa() + (V + day_indexes[i])->get_ith_outgoing(1)->get_capa() << "  ";
    }
    std::cout << std::endl;
}

void Graph::check_flow() {
    for (unsigned int i = 0; i < n_V; i++){
        if (i == source_index || i == sink_index) continue;
        unsigned int f1 = 0, f2 = 0;
        Vertex* w = V + i;
        unsigned int size_d_plus = w->get_delta_p_size();
        for (unsigned int j = 0; j < size_d_plus; j++) {
            f1 += w->get_ith_outgoing(j)->get_flow();
        }
        unsigned int size_d_minus = w->get_delta_m_size();
        for (unsigned int j = 0; j < size_d_minus; j++) {
            f2 += w->get_ith_incoming(j)->get_flow();
        }
        if (f1 != f2)
            throw std::logic_error("Why is my flow not ok ?");
    }

    for (unsigned int i = 0; i < n_A; i++){
        Arc* a = A + i;
        if (a->get_flow() > a->get_capa())
            throw std::logic_error("Why is my flow greater than my capa ?");
    }
}
