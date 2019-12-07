#include "Flow.h"

#include "tools.h"

Graph::Graph(const std::vector<unsigned int> &schedule, unsigned int epsilon) {
    // build the graph from right to left
    n_V = 1 + NB_FAMILIES + NB_DAYS*NB_FAMILIES + 2*NB_DAYS + 1;
    V = new Vertex [n_V];
    n_A = NB_FAMILIES + NB_FAMILIES*NB_DAYS*3 + NB_DAYS*3;
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
        A[next_arc_index] = Arc(next_arc_index, V + next_vertex_index, V, schedule[i] - epsilon, -UPPER_BOUND);
        next_arc_index++;
        A[next_arc_index] = Arc(
                next_arc_index,
                V + next_vertex_index,
                V + next_vertex_index - 1,
                2*epsilon,
                0);
        next_arc_index++;
        day_indexes.push_back(next_vertex_index);
        next_vertex_index++;
    }

    // the families
    std::vector<std::vector<unsigned int>> family_data = read_instance(INSTANCE_PATH);
    unsigned int choice_number;
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        V[next_vertex_index] = Vertex(next_vertex_index);
        family_indexes.push_back(next_vertex_index);
        next_vertex_index++;
        for (unsigned int j = 0; j < NB_DAYS; j++){
            choice_number = 10;
            for (unsigned int k = 0; k < NB_CHOICES; k++)
                if (family_data[i][k] == j)
                    choice_number = k;
            V[next_vertex_index] = Vertex(next_vertex_index);
            A[next_arc_index] = Arc(next_arc_index, V + next_vertex_index, V + day_indexes[j], family_data[i][NB_CHOICES] - 1, 0);
            next_arc_index++;
            A[next_arc_index] = Arc(next_arc_index, V + family_indexes[i], V + next_vertex_index, family_data[i][NB_CHOICES] - 1, MARGINAL_COST[choice_number]);
            next_arc_index++;
            A[next_arc_index] = Arc(next_arc_index, V + family_indexes[i], V + day_indexes[j], 1, CONSTANT_COST[choice_number]);
            next_arc_index++;
            next_vertex_index++;
        }
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
    std::vector<unsigned int> distances(n_V, 10*UPPER_BOUND);
    distances[source_index] = 0;
    std::vector<unsigned int> predecessor(n_V);
    //std::vector<bool> has_changed(n_V, true);
    bool changed = true;
    for (unsigned int i = 0; changed && i < n_V; i++){
        changed = false;
        for (unsigned int j = 0; j < n_A; j++){
            Arc* a = A + j;
            if (a->has_capacity_left() && distances[a->get_u()->get_id()] + a->get_cost() < distances[a->get_v()->get_id()]){
                changed = true;
                distances[a->get_v()->get_id()] = distances[a->get_u()->get_id()] + a->get_cost();
                predecessor[a->get_v()->get_id()] = j;
            }
            if (a->has_flow() && distances[a->get_v()->get_id()] - a->get_cost() < distances[a->get_u()->get_id()]){
                changed = true;
                distances[a->get_u()->get_id()] = distances[a->get_v()->get_id()] - a->get_cost();
                predecessor[a->get_u()->get_id()] = j;
            }
        }
    }

    if (distances[sink_index] == 10*UPPER_BOUND) return false;

    unsigned int current_index = sink_index;
    unsigned int additional_flow = 10000; // big
    while (current_index != source_index){
        Arc* a = A + predecessor[current_index];
        if (a->get_v()->get_id() == current_index){
            current_index = a->get_u()->get_id();
            additional_flow = std::min(additional_flow, a->get_capa() - a->get_flow());
        }
        else {
            current_index = a->get_v()->get_id();
            additional_flow = std::min(additional_flow, a->get_flow());
        }
    }
    current_index = sink_index;
    while (current_index != source_index){
        Arc* a = A + predecessor[current_index];
        if (a->get_v()->get_id() == current_index){
            current_index = a->get_u()->get_id();
            a->add_to_flow(additional_flow);
        }
        else{
            current_index = a->get_v()->get_id();
            a->add_to_flow(-additional_flow);
        }

    }
    std::cout<<"ok " << additional_flow << std::endl;
    return true;
}

void Graph::compute_max_flow_min_cost() {
    while(find_and_apply_augmenting_path()){}
}
