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
}

int Graph::init_distances_and_predecessors() {
    distances = std::vector<int> (n_V, 10*UPPER_BOUND);
    predecessor = std::vector<Arc*>(n_V);
}

Graph::Graph(const std::vector<std::vector<unsigned int>> &family_data) {
    // build the graph from right to left
    unsigned int k_max = 10;
    n_V = 1 + NB_FAMILIES + NB_DAYS*NB_FAMILIES + 2*NB_DAYS + 1;
    V = new Vertex [n_V];
    n_A = NB_FAMILIES + NB_FAMILIES*k_max*3 + NB_DAYS*3;
    A = new Arc [n_A];
    unsigned int next_vertex_index = 0,
            next_arc_index = 0;
    family_indexes.clear(); day_indexes.clear();

    // the sink
    V[next_vertex_index] = Vertex(next_vertex_index);
    sink_index = next_vertex_index;
    next_vertex_index++;

    // the days
    for (unsigned int i = 0; i < NB_DAYS; i++){
        V[next_vertex_index] = Vertex(next_vertex_index);
        A[next_arc_index] = Arc(
                next_arc_index,
                V + next_vertex_index,
                V,
                MAX_NB_PEOPLE_PER_DAY - MIN_NB_PEOPLE_PER_DAY,
                0);
        next_vertex_index++;
        next_arc_index++;
        V[next_vertex_index] = Vertex(next_vertex_index);
        A[next_arc_index] = Arc(next_arc_index, V + next_vertex_index, V, MIN_NB_PEOPLE_PER_DAY, -UPPER_BOUND);
//        A[next_arc_index] = Arc(next_arc_index, V + next_vertex_index, V, schedule[i] - epsilon, -UPPER_BOUND);
        next_arc_index++;
        A[next_arc_index] = Arc(
                next_arc_index,
                V + next_vertex_index,
                V + next_vertex_index - 1,
//                2*epsilon,
                MAX_NB_PEOPLE_PER_DAY - MIN_NB_PEOPLE_PER_DAY,
                0);
        next_arc_index++;
        day_indexes.push_back(next_vertex_index);
        next_vertex_index++;
    }

    // the families
    //unsigned int choice_number;
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        V[next_vertex_index] = Vertex(next_vertex_index);
        family_indexes.push_back(next_vertex_index);
        next_vertex_index++;
        for (unsigned int k = 0; k < k_max; k++){
            unsigned int j = family_data[i][k];
            V[next_vertex_index] = Vertex(next_vertex_index);
            A[next_arc_index] = Arc(next_arc_index, V + next_vertex_index, V + day_indexes[j], family_data[i][NB_CHOICES] - 1, 0);
            next_arc_index++;
            A[next_arc_index] = Arc(next_arc_index, V + family_indexes[i], V + next_vertex_index, family_data[i][NB_CHOICES] - 1, MARGINAL_COST[k]);
            next_arc_index++;
            A[next_arc_index] = Arc(next_arc_index, V + family_indexes[i], V + day_indexes[j], 1, CONSTANT_COST[k]);
            next_arc_index++;
            next_vertex_index++;
        }
//        for (unsigned int j = 0; j < NB_DAYS; j++){
//            choice_number = 10;
//            for (unsigned int k = 0; k < NB_CHOICES; k++)
//                if (family_data[i][k] == j)
//                    choice_number = k;
//            V[next_vertex_index] = Vertex(next_vertex_index);
//            A[next_arc_index] = Arc(next_arc_index, V + next_vertex_index, V + day_indexes[j], family_data[i][NB_CHOICES] - 1, 0);
//            next_arc_index++;
//            A[next_arc_index] = Arc(next_arc_index, V + family_indexes[i], V + next_vertex_index, family_data[i][NB_CHOICES] - 1, MARGINAL_COST[choice_number]);
//            next_arc_index++;
//            A[next_arc_index] = Arc(next_arc_index, V + family_indexes[i], V + day_indexes[j], 1, CONSTANT_COST[choice_number]);
//            next_arc_index++;
//            next_vertex_index++;
//        }
    }

    // the source
    V[next_vertex_index] = Vertex(next_vertex_index);
    source_index = next_vertex_index;
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        A[next_arc_index] = Arc(next_arc_index, V + next_vertex_index, V + family_indexes[i], family_data[i][NB_CHOICES], 0);
        next_arc_index++;
    }

    // build the deltas
    for (unsigned int i = 0; i < n_A; i++){
        A[i].fill_vertex_deltas();
    }
    init_distances_and_predecessors();

    std::cout << next_arc_index << " " << n_A << " " << next_vertex_index << " " << n_V << std::endl;
}

std::vector<unsigned int> Graph::get_solution() const {
    std::vector<unsigned int> solution (NB_FAMILIES, 101);
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        Vertex* u = V + family_indexes[i];
        for (unsigned int j = 0; j < u->get_delta_p_size(); j++){
            Arc* a = u->get_ith_outgoing(j);
            if (a->has_flow()){
                unsigned int day_id = a->get_v()->get_id();
                for (unsigned int k = 0; k < day_indexes.size(); k++){
                    if (day_indexes[k] == day_id){
                        solution[i] = k;
                        break;
                    }
                }
                break;
            }
        }
    }
    return solution;
}

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

    return true;
}

void Graph::compute_max_flow_min_cost() {
    unsigned int i = 0;
    while(find_and_apply_augmenting_path()){
        i++;
        if (i%30 == 0)
            std::cout << "flow is " << get_current_flow() << " with cost " << get_flow_cost() << std::endl;
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

std::vector<Arc *> Graph::get_shortest_path() const {
    distances[source_index] = 0;
    std::queue<Vertex*> Q;
    Q.push(V + source_index);
    Arc* a; Vertex *u, *v, *w;
    unsigned int size_d_plus;
    unsigned int size_d_minus;
    // I assume this is not an infinite loop (no absorbing circuit)
    while (!Q.empty()){
        w = Q.front();
        Q.pop();
        size_d_plus = w->get_delta_p_size();
        for (unsigned int i = 0; i < size_d_plus; i++){
            a = w->get_ith_outgoing(i);
            v = a->get_v();
            if (a->has_capacity_left() && distances[w->get_id()] + a->get_cost() < distances[v->get_id()]){
                distances[v->get_id()] = distances[w->get_id()] + a->get_cost();
                Q.push(v);
                predecessor[v->get_id()] = a;
            }
        }
        size_d_minus = w->get_delta_m_size();
        for (unsigned int i = 0; i < size_d_minus; i++){
            a = w->get_ith_incoming(i);
            u = a->get_u();
            if (a->has_flow() && distances[w->get_id()] - a->get_cost() < distances[u->get_id()]){
                distances[u->get_id()] = distances[w->get_id()] - a->get_cost();
                Q.push(u);
                predecessor[u->get_id()] = a;
            }
        }
    }

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

int Graph::get_flow_cost() const {
    int res = 0;
    for (unsigned int i = 0; i < n_A; i++)
        res += (A + i)->get_flow() * (A + i)->get_cost();
    return res;
}

int Graph::get_true_flow_cost() const {
    unsigned int res = 0;
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        Vertex* u = V + family_indexes[i];
        unsigned int d_out = u->get_delta_p_size();
        for (unsigned int k = 0; k < d_out; k++){
            Arc* a = u->get_ith_outgoing(k);
            res += a->get_flow() * a->get_cost();
        }
    }
    return res;
}

