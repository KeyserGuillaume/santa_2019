#include "Flow.h"


Graph::Graph() {
    n_V = 5;
    V = new Vertex [n_V];
    n_A = 5;
    A = new Arc [n_A];
    unsigned int next_vertex_index = 0;
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

void Graph::init_distances_and_predecessors() {
    distances = std::vector<int> (n_V, 10*UPPER_BOUND);
    predecessor = std::vector<Arc*>(n_V);
    distances[source_index] = 0;
    std::queue<Vertex*> Q;
    Q.push(V + source_index);
    // I assume this is not an infinite loop (no absorbing circuit)
    while (!Q.empty()) {
        apply_Bellman_Ford(Q);
    }
    if (distances[sink_index] == 0)
        can_we_add_zeros_valued_flows = true;
}

void Graph::get_affluence_bounds(const std::vector<std::vector<unsigned int>> &family_data, const std::vector<preset> &presets, std::vector<unsigned int> &lower_bounds, std::vector<unsigned int> &upper_bounds) const {
    std::vector<unsigned int> preset_lower_bound(NB_DAYS, 0);
    std::vector<unsigned int> preset_upper_bound(NB_DAYS, 0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        for (unsigned int k = 0; k < K_MAX; k++) {
            if (presets[i][k] == COMPULSORY) {
                preset_lower_bound[family_data[i][k]] += family_data[i][NB_CHOICES];
                preset_upper_bound[family_data[i][k]] += family_data[i][NB_CHOICES];
            }
            if (presets[i][k] == ALLOWED) {
                preset_upper_bound[family_data[i][k]] += family_data[i][NB_CHOICES];
            }
        }
    }

    for (unsigned int i = 0; i < NB_DAYS; i++){
        upper_bounds[i] = std::min(MAX_NB_PEOPLE_PER_DAY, preset_upper_bound[i]);
        lower_bounds[i] = std::max(MIN_NB_PEOPLE_PER_DAY, preset_lower_bound[i]);
    }

    for (unsigned int i = 0; i < NB_DAYS; i++){
        if (i < NB_DAYS - 1 && lower_bounds[i] > MIN_NB_PEOPLE_PER_DAY && upper_bounds[i] <= lower_bounds[i + 1]){
            //std::cout << upper_bounds[i] << " " << lower_bounds[i + 1]<<std::endl;
            unsigned int a;
            for (unsigned int j = 0; j < 10; j++) {
                a = lower_bounds[i];
                lower_bounds[i] = std::max(a, (unsigned int)(std::max(0., ceil(lower_bounds[i + 1] - 50./log(a)*(log(400.*1000./float(a - MIN_NB_PEOPLE_PER_DAY) - 0.5))))));
            }
            a = lower_bounds[i];
            upper_bounds[i + 1] = std::min(upper_bounds[i + 1], (unsigned int)(floor(upper_bounds[i] + 50./log(a)*(log(400.*1000./float(a - MIN_NB_PEOPLE_PER_DAY) - 0.5)))));
            //std::cout << preset_lower_bound[i] << " " << lower_bounds[i] << " " << upper_bounds[i] << " " << lower_bounds[i + 1] << " " << upper_bounds[i + 1] << " " << preset_upper_bound[i + 1] << std::endl;
        }
    }
}

std::vector<uint_pair> Graph::get_bottleneck_bounds(const std::vector<std::vector<unsigned int>> &family_data,
                                                    const std::vector<preset> &presets,
                                                    const std::vector<unsigned int> &presets_schedule,
                                                    std::vector<unsigned int> lower_bounds,
                                                    std::vector<unsigned int> upper_bounds) const{
    // this gives the bottleneck bound. To know how many families can be accomodated during the day, you need to add the
    // nb of families already put on that day by the preset.
    // first is the min number of families that should be assigned to the day, second is the max number.
    std::vector<uint_pair> sorted_families(0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        sorted_families.push_back(uint_pair(i, family_data[i][NB_CHOICES]));
    }
    sort_by_second(sorted_families);
    unsigned int current_family, current_day;
    std::vector<uint_pair> bottleneck_bounds(NB_DAYS, uint_pair(0, 0));
    std::vector<unsigned int> possible_family_cum_sizes(NB_DAYS, 0);
    std::vector<bool> day_has_its_value(NB_DAYS, false);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        for (unsigned int k = 0; k < K_MAX; k++) {
            current_family = sorted_families[i].first;
            current_day = family_data[current_family][k];
            if (presets[current_family][k] == ALLOWED && !day_has_its_value[current_day]){
                possible_family_cum_sizes[current_day] += family_data[current_family][NB_CHOICES];
                if (possible_family_cum_sizes[current_day] > upper_bounds[current_day] - presets_schedule[current_day])
                    day_has_its_value[current_day] = true;
                else  // there is a good reason why this is not identical to the stuff below
                    bottleneck_bounds[current_day].second++;
            }
        }
    }
    for (unsigned int i = 0; i < NB_DAYS; i++)
        std::cout << bottleneck_bounds[i].first << std::string(4 - nb_chiffres(bottleneck_bounds[i].first),' ');
    std::cout << std::endl;
    possible_family_cum_sizes = std::vector<unsigned int>(NB_DAYS, 0);
    day_has_its_value = std::vector<bool> (NB_DAYS, false);
    for (unsigned int i = 0; i < NB_FAMILIES; i++) {
        for (unsigned int k = 0; k < K_MAX; k++) {
            current_family = sorted_families[NB_FAMILIES - i - 1].first;
            current_day = family_data[current_family][k];
            if (presets[current_family][k] == ALLOWED && !day_has_its_value[current_day]){
                possible_family_cum_sizes[current_day] += family_data[current_family][NB_CHOICES];
                bottleneck_bounds[current_day].first++;
                if (possible_family_cum_sizes[current_day] >= lower_bounds[current_day] - presets_schedule[current_day])
                    day_has_its_value[current_day] = true;
            }
        }
    }
    for (unsigned int i = 0; i < NB_DAYS; i++)
        std::cout << bottleneck_bounds[i].first << std::string(4 - nb_chiffres(bottleneck_bounds[i].first),' ') << bottleneck_bounds[i].second << std::endl;
    return bottleneck_bounds;
}


Graph::Graph(const std::vector<std::vector<unsigned int>> &family_data, std::vector<preset> presets) {
    // build the graph from right to left

    unsigned int nb_already_assigned = 0;
    unsigned int nb_forbidden_assignments = 0;
    std::vector<bool> is_already_assigned (NB_FAMILIES, false);
    std::vector<unsigned int> lower_bounds(NB_DAYS, 0);
    std::vector<unsigned int> upper_bounds(NB_DAYS, 0);
    get_affluence_bounds(family_data, presets, lower_bounds, upper_bounds);

    std::vector<unsigned int> presets_schedule(NB_DAYS, 0);
    std::vector<unsigned int> presets_cardinal(NB_DAYS, 0);
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        for (unsigned int k = 0; k < K_MAX; k++) {
            if (presets[i][k] == COMPULSORY) {
                presets_costs += CONSTANT_COST[k] + family_data[i][NB_CHOICES] * MARGINAL_COST[k];
                is_already_assigned[i] = true;
                nb_already_assigned++;
                presets_schedule[family_data[i][k]] += family_data[i][NB_CHOICES];
                presets_cardinal[family_data[i][k]] += 1;
            }
            if (presets[i][k] != ALLOWED)
                nb_forbidden_assignments++;
        }
    }
    n_V = 1 + 3*(NB_FAMILIES - nb_already_assigned) + 2*NB_DAYS + 1;
    V = new Vertex [n_V];
    n_A = 3*(NB_FAMILIES - nb_already_assigned) + 2*NB_FAMILIES*K_MAX - 2*nb_forbidden_assignments + 4*NB_DAYS;
    A = new Arc [n_A];
    unsigned int next_vertex_index = 0, next_arc_index = 0;
    family_indexes.clear(); full_family_indexes = std::vector<int>(NB_FAMILIES, -1);
    day_indexes.clear();

    // the sink
    V[next_vertex_index] = Vertex(next_vertex_index);
    sink_index = next_vertex_index;
    next_vertex_index++;
    
    std::vector<unsigned int> bottleneck_indexes(0);
    std::vector<uint_pair> bottleneck_bounds = get_bottleneck_bounds(family_data, presets, presets_schedule, lower_bounds, upper_bounds);

    // the days
    for (unsigned int i = 0; i < NB_DAYS; i++){
        if (presets_schedule[i] > MAX_NB_PEOPLE_PER_DAY) throw std::logic_error("Why do I have too many people assigned to this day ?");
        unsigned int min_allowed = (unsigned int)(std::max(int(lower_bounds[i]) - int(presets_schedule[i]), 0));
        unsigned int max_allowed = upper_bounds[i] - presets_schedule[i];
//        A[next_arc_index] = Arc(next_arc_index, V + next_vertex_index, V, MAX_NB_PEOPLE_PER_DAY - MIN_NB_PEOPLE_PER_DAY, 0);
        V[next_vertex_index] = Vertex(next_vertex_index);
        day_indexes.push_back(next_vertex_index);
        A[next_arc_index] = Arc(next_arc_index, V + day_indexes[i], V + sink_index, min_allowed, -UPPER_BOUND);
        next_arc_index++;
        A[next_arc_index] = Arc(next_arc_index, V + day_indexes[i], V + sink_index, max_allowed - min_allowed, 0);
//                2*epsilon, MAX_NB_PEOPLE_PER_DAY - MIN_NB_PEOPLE_PER_DAY,
        next_arc_index++;
        next_vertex_index++;
        V[next_vertex_index] = Vertex(next_vertex_index);
        bottleneck_indexes.push_back(next_vertex_index);
        next_vertex_index++;
        A[next_arc_index] = Arc(next_arc_index, V + bottleneck_indexes[i], V + day_indexes[i], bottleneck_bounds[i].first, -UPPER_BOUND);
        next_arc_index++;
        A[next_arc_index] = Arc(next_arc_index, V + bottleneck_indexes[i], V + day_indexes[i], bottleneck_bounds[i].second - bottleneck_bounds[i].first, 0);
        next_arc_index++;
        std::cout << presets_cardinal[i] << std::string(4 - nb_chiffres(presets_cardinal[i]), ' ');
    }

    // the families
    for (unsigned int i = 0; i < NB_FAMILIES; i++){
        if (is_already_assigned[i]) continue;
        V[next_vertex_index] = Vertex(next_vertex_index);
        unsigned int family_const_cost_middleman = next_vertex_index;
        next_vertex_index++;
        V[next_vertex_index] = Vertex(next_vertex_index);
        unsigned int family_marg_cost_middleman = next_vertex_index;
        next_vertex_index++;
        V[next_vertex_index] = Vertex(next_vertex_index);
        family_indexes.push_back(next_vertex_index);
        full_family_indexes[i] = next_vertex_index;
        unsigned int family_index = next_vertex_index;
        next_vertex_index++;
        A[next_arc_index] = Arc(next_arc_index, V + family_index, V + family_const_cost_middleman, 1, 0);
        next_arc_index++;
        A[next_arc_index] = Arc(next_arc_index, V + family_index, V + family_marg_cost_middleman, family_data[i][NB_CHOICES] - 1, 0);
        next_arc_index++;
        for (unsigned int k = 0; k < K_MAX; k++){
            if (presets[i][k] != ALLOWED) continue;
            unsigned int j = family_data[i][k];
            A[next_arc_index] = Arc(next_arc_index, V + family_const_cost_middleman, V + bottleneck_indexes[j], 1, CONSTANT_COST[k]);
            next_arc_index++;
            A[next_arc_index] = Arc(next_arc_index, V + family_marg_cost_middleman, V + day_indexes[j], family_data[i][NB_CHOICES] - 1, MARGINAL_COST[k]);
            next_arc_index++;
        }
    }

    // the source
    V[next_vertex_index] = Vertex(next_vertex_index);
    source_index = next_vertex_index;
    next_vertex_index++;
    for (unsigned int i = 0; i < full_family_indexes.size(); i++){
        if (full_family_indexes[i] == -1) continue;
        A[next_arc_index] = Arc(next_arc_index, V + source_index, V + full_family_indexes[i], family_data[i][NB_CHOICES], 0);
        next_arc_index++;
    }

    std::cout << next_arc_index << " " << n_A << " " << next_vertex_index << " " << n_V << std::endl;

    // build the deltas
    for (unsigned int i = 0; i < n_A; i++){
        A[i].fill_vertex_deltas();
    }

    max_possible_flow = 0;
    for (unsigned int i = 0; i < NB_FAMILIES; i++)
        if (!is_already_assigned[i])
            max_possible_flow += family_data[i][NB_CHOICES];

    add_m2M_valued_flows();
    check_flow(false);
    init_distances_and_predecessors();
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
        if (can_we_add_zeros_valued_flows)
            add_zero_valued_flows();
        if (can_we_add_mM_valued_flows)
            add_mM_valued_flows();
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
    // if a vertex is at finite distance from source and the arc that made this distance possible is now saturated,
    // then this vertex needs an update on its distance
    for (unsigned int i = 0; i < n_V; i++){
        if (i == source_index) continue;
        if (distances[i] == 10 * UPPER_BOUND) continue;
        w = V + i;
        a = predecessor[i];
        if ((a->get_v()->get_id() == i && !a->has_capacity_left()) || (a->get_u()->get_id() == i && !a->has_flow())) {
            Q.push(w);
        }
    }
    // Q contains vertexes whose distance may not be supported. If current vertex distance is not supported, set
    // distance to infinity and push all neighbors in Q.
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

    // now Q is the bellman-ford queue (not exactly bellman-ford because it loops forever if any absorbing circuit...
    // we push in Q all vertices who have a neighbor whose distance needs to be recomputed
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

    if (distances[sink_index] == 0)
        can_we_add_zeros_valued_flows = true;

    if (distances[sink_index] == -UPPER_BOUND)
        can_we_add_mM_valued_flows = true;

    Vertex *w; Arc* a;

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
        Vertex* family_const_cost_middleman = (V + family_indexes[i])->get_ith_outgoing(0)->get_v();
        Vertex* family_marg_cost_middleman = (V + family_indexes[i])->get_ith_outgoing(1)->get_v();
        unsigned int d_out = family_const_cost_middleman->get_delta_p_size();
        for (unsigned int k = 0; k < d_out; k++){
            Arc* a = family_const_cost_middleman->get_ith_outgoing(k);
            res += a->get_flow() * a->get_cost();
            a = family_marg_cost_middleman->get_ith_outgoing(k);
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

void Graph::check_flow(const bool &check_minimal_flow) {
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

    for (unsigned int i = 0; i < full_family_indexes.size(); i++)
        if (full_family_indexes[i] != -1 && (V + full_family_indexes[i])->get_ith_incoming(0)->has_capacity_left() && check_minimal_flow)
            throw std::logic_error("Why do I have this capacity left ?");

    for (unsigned int i = 0; i < n_A; i++){
        Arc* a = A + i;
        if (a->get_flow() > a->get_capa())
            throw std::logic_error("Why is my flow greater than my capa ?");
        if (a->get_cost() == -UPPER_BOUND && a->has_capacity_left() && check_minimal_flow)
            throw std::logic_error("Why do I have this capacity left ?");
    }

}

void Graph::add_m2M_valued_flows() {

    // flows of value - 2*UPPER_BOUND
    for (unsigned int i = 0; i < family_indexes.size(); i++)
        add_obvious_flow(family_indexes[i], std::vector<unsigned int>({0, 0, 0, 0}));
}

void Graph::add_mM_valued_flows() {
    // flows of value - UPPER_BOUND
    for (unsigned int i = 0; i < family_indexes.size(); i++)
        add_obvious_flow(family_indexes[i], std::vector<unsigned int>({0, 0, 0, 1}));
    for (unsigned int i = 0; i < family_indexes.size(); i++)
        add_obvious_flow(family_indexes[i], std::vector<unsigned int>({0, 0, 1, 0}));
    for (unsigned int i = 0; i < family_indexes.size(); i++)
        add_obvious_flow(family_indexes[i], std::vector<unsigned int>({1, 0, 0}));
    for (unsigned int i = 0; i < family_indexes.size(); i++)
        add_obvious_flow(family_indexes[i], std::vector<unsigned int>({1, 1, 0}));

}

void Graph::add_zero_valued_flows() {
    // flows of value 0
    for (unsigned int i = 0; i < family_indexes.size(); i++)
        add_obvious_flow(family_indexes[i], std::vector<unsigned int> ({0, 0, 1, 1}));
    for (unsigned int i = 0; i < family_indexes.size(); i++)
        add_obvious_flow(family_indexes[i], std::vector<unsigned int> ({1, 0, 1}));
    for (unsigned int i = 0; i < family_indexes.size(); i++)
        add_obvious_flow(family_indexes[i], std::vector<unsigned int> ({1, 1, 1}));
}

void Graph::add_obvious_flow(const unsigned int &family_index, const std::vector<unsigned int> &turns, const unsigned int &flow_quantity, const bool &force) {
    std::vector<Arc*> path(0);
    Vertex* w = V + family_index;
    Arc* a = w->get_ith_incoming(0);
    if (!a->has_capacity_left()) return;
    path.push_back(a);
    unsigned int additional_flow = a->get_capa() - a->get_flow();
    for (unsigned int i = 0; i < turns.size(); i++){
        a = w->get_ith_outgoing(turns[i]);
        if ((!force && a->get_cost() > 0) || !a->has_capacity_left()) return;
        path.push_back(a);
        additional_flow = std::min(additional_flow, a->get_capa() - a->get_flow());
        w = a->get_v();
    }
    if (flow_quantity)
        additional_flow = flow_quantity;
    for (unsigned int i = 0; i < path.size(); i++)
        path[i]->add_to_flow(additional_flow);
}

void Graph::clear_flow() {
    for (unsigned int i = 0; i < n_A; i++)
        (A + i)->clear_flow();
}

void Graph::add_flow_for_assigning_family(const unsigned int &family_idx, const unsigned int &k, const unsigned int &n_people) {
    //std::cout << "This add_flow_for_assigning_family function will give unreliable results if the presets contain any FORBIDDEN values" << std::endl;
    if (k >= K_MAX) throw std::invalid_argument("k is greater than K_MAX");
    if (full_family_indexes[family_idx] == 1) throw std::invalid_argument("You can't assign flow for a family already assigned");
    if ((V + full_family_indexes[family_idx])->get_ith_outgoing(k + 1)->get_v()->get_ith_outgoing(0)->has_capacity_left())
        add_obvious_flow(full_family_indexes[family_idx], std::vector<unsigned int>({0, k, 0, 0}), 1, true);
    else
        add_obvious_flow(full_family_indexes[family_idx], std::vector<unsigned int>({0, k, 0, 1, 0}), 1, true);
    Arc* a = (V + full_family_indexes[family_idx])->get_ith_outgoing(2*k + 1)->get_v()->get_ith_outgoing(0);
    unsigned int value = std::min(n_people - 1, a->get_capa() - a->get_flow());
    if (a->has_capacity_left())
        add_obvious_flow(full_family_indexes[family_idx], std::vector<unsigned int>({1, k, 0}), value, true);
    add_obvious_flow(full_family_indexes[family_idx], std::vector<unsigned int>({1, k, 1, 0}), n_people - 1 - value, true);
    Vertex* w = V + full_family_indexes[family_idx];
   // if (family_idx == 10)
   //     std::cout << n_people << " " << value << " " <<a->get_capa() << " " << a->get_flow() << std::endl;
}

void Graph::inspect_distances() const {
    std::cout << "source distance: " << distances[source_index] << std::endl;
    std::cout << "family distances: ";
    for (unsigned int i = 0; i < family_indexes.size(); i++)
        std::cout << distances[family_indexes[i]] << std::string(10 - nb_chiffres(family_indexes[i]),' ');
    std::cout << std::endl << "day distances: ";
    for (unsigned int i = 0; i < day_indexes.size(); i++)
        std::cout << distances[day_indexes[i]] << std::string(10 - nb_chiffres(day_indexes[i]),' ');
    std::cout << "sink distance: " << distances[sink_index] << std::endl;
}
