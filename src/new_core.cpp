#include <queue>
#include <set>
#include <tuple>
#include <string.h>
#include "new_core.hpp"
#include "robin_hood.h"

algo::Core::Core(void)
{

}

/**
 * The function `rndm` generates a random number between `a` and `b` using a power law distribution
 * with exponent `g`, and returns the rounded result.
 * 
 * @param a The parameter "a" is an integer value.
 * @param b The parameter "b" represents the upper bound of the range from which the random number will
 * be generated.
 * @param g The parameter "g" represents the exponent used in the calculation. It is used to raise the
 * values of "a" and "b" to the power of "g" in the formula.
 * @param size The "size" parameter is not used in the given code snippet. It is not clear what it is
 * intended for.
 * 
 * @return a rounded value calculated using the given parameters.
 */
int algo::rndm(int a, int b, double g) {
    std::random_device eng;
    using dist_params = typename std::uniform_int_distribution<int>::param_type;
    std::uniform_int_distribution<int> dist (0, INT_MAX);
    double r = (double) dist(eng) / INT_MAX;
    double ag = pow(a, g);
    double bg = pow(b, g);
    return round(pow((ag + (bg - ag)*r), (double) (1./g)));
}

/**
 * The function `init_data` reads data from a file, processes it, and initializes the `graph` object
 * with the data. Works for datasets: WQ, SO and LK.
 * 
 * @param src The parameter `src` is a `std::string` that represents the file path of the input file.
 * 
 * @return Nothing is being returned.
 */
void algo::Core::init_data(std::string const &src)
{
    std::ifstream file(src);
    char trash;
    int max_src = 0, max_dst = 0;

    if (!file.is_open())
        throw std::runtime_error("couldn't open file:" + src);
    while (!file.eof()) {
        algo::edge e; double ts; file >> e.src >> e.dst >> trash >> ts;
        e.start_time = (int) ts;
        e.src--; e.dst--;
        max_src = std::max(max_src, e.src);
        max_dst = std::max(max_dst, e.dst);
        e.end_time = e.start_time + rndm(3600, 172800, -2.5); // 1 hour to 48 hours
        edges.push_back(e);
    }
    std::sort(edges.begin(), edges.end(), [](edge const &a, edge const &b) {
        return a.start_time < b.start_time;
    });

    graph.init(max_src, max_dst);

    for (int i = 0; i < edges.size() ; i++) {
        graph.upper[edges[i].src].edges.push_back(i);
        graph.lower[edges[i].dst].edges.push_back(i);
    }
}

/**
 * The function "get_vertex_order" returns a vector of integers representing the order of vertices in a
 * graph, sorted in descending order based on the number of edges each vertex has.
 * 
 * @param graph temporal bipartie graph
 * 
 * @return a std::vector<int> containing the order of vertices in the graph.
 */
std::vector<int> algo::Core::get_vertex_order(algo::Graph const &graph)
{
    std::vector<int> order;

    for (int i = 0; i < graph.upper.size(); i++)
        order.push_back(i);

    std::sort(order.begin(), order.end(), [&](int a, int b) -> bool {
        return graph.upper[a].edges.size() > graph.upper[b].edges.size();
    });
    return order;
}

/**
 * The function "get_partitions_of_node" takes a vector of edges and a node, and returns a vector of
 * partitions of the node's edges, the edges in each partition overlap with at least one other edge in the partition
 * and edges from different partitions don't overlap.
 * 
 * @param edges The parameter `edges` is a vector of `algo::edge` objects.
 * @param node The "node" parameter is of type "algo::Node&", which means it is a reference to an
 * object of type "algo::Node".
 * 
 * @return a vector of partitions, where each partition is represented by a vector of edges.
 */
std::vector<algo::partition> algo::Core::get_partitions_of_node(std::vector<algo::edge> const &edges, algo::Node &node)
{
    if (node.edges.size() == 0) {
        return {};
    }
    std::vector<algo::partition> partitions = {{node.edges[0]}};
    int actual_end_time = edges[partitions[0][0]].end_time;

    for (int i = 1 ; i < node.edges.size() ; i++) {
        if (edges[node.edges[i]].start_time < actual_end_time)
            partitions.rbegin()->push_back(node.edges[i]);
        else
            partitions.push_back({node.edges[i]});
        actual_end_time = std::max(actual_end_time, edges[node.edges[i]].end_time);
    }
    return partitions;
}

/**
 * The function `TBP_build` builds a TBP_index based on a given vertex order and graph.
 * 
 * @param vertex_order
 * @param graph temporal bipartie graph
 * 
 * @return a TBP-index.
 */
algo::TBP_index algo::Core::TBP_build(std::vector<int> const &vertex_order, algo::Graph &graph)
{
    algo::TBP_index index(graph.upper.size(), graph.upper.size());

    for (int i = 0; i < graph.lower.size() ; i++) {
        graph.lower[i].partitions = get_partitions_of_node(edges, graph.lower[i]);
    }
    for (int i = 0; i < vertex_order.size() ; i++) {
        std::vector<bool> visited(edges.size());
        std::vector<int> minT(graph.upper.size(), INT_MAX);
        {
            auto cmp = [](entry const &a, entry const &b) {return a.start_time < b.start_time;};
            std::priority_queue<algo::entry, std::vector<algo::entry>, decltype(cmp)> q(cmp);
            auto &uk = graph.upper[vertex_order[i]];

            minT[vertex_order[i]] = 0;

            for (int j = uk.edges.size() - 1 ; j >= 0 ; j--) {
                if (j && edges[uk.edges[j]].start_time == edges[uk.edges[j - 1]].start_time) continue;
                q.push({vertex_order[i], edges[uk.edges[j]].start_time, edges[uk.edges[j]].start_time});

                while (!q.empty()) {
                    auto w = q.top(); q.pop();
                    if (w.node != uk.id) {
                        if (TBP_query(index, uk, graph.upper[w.node], w.start_time, w.end_time)) {
                            continue;
                        }
                        index.L_in[w.node].push_back({uk.id, w.start_time, w.end_time});
                    }

                    auto &w_edges = graph.upper[w.node].edges;
                    for (int k = 0 ; k < w_edges.size() ; k++) {
                        if (!(!visited[w_edges[k]] && edges[w_edges[k]].start_time >= w.end_time)) continue;
                        visited[w_edges[k]] = true;
                        algo::partition const *partition = find_partition(graph.lower[edges[w_edges[k]].dst].partitions, w_edges[k]);
                        if (!partition) continue;
                        for (int l = 0; l < partition->size() ; l++) {
                            auto &edge_id = (*partition)[l];
                            if (edge_id == w_edges[k]) {
                                continue;
                            }
                            if (graph.upper[edges[edge_id].src].order <= i) {
                                continue;
                            }
                            if (edges[edge_id].end_time >= minT[edges[edge_id].src]) {
                                continue;
                            }
                            if (edge_colliding(w_edges[k], edge_id)) {
                                q.push({edges[edge_id].src, w.start_time, edges[edge_id].end_time});
                                minT[edges[edge_id].src] = edges[edge_id].end_time;
                            }
                        }
                    }
                }
            }
        }

        std::vector<bool> visited1(edges.size());
        std::vector<int> maxT(graph.upper.size());
        {
            auto cmp = [](entry const &a, entry const &b) {return a.end_time > b.end_time;};
            std::priority_queue<algo::entry, std::vector<algo::entry>, decltype(cmp)> q(cmp);
            auto &uk = graph.upper[vertex_order[i]];

            maxT[vertex_order[i]] = INT_MAX;
            std::set<int> sorter;
            for(int j = 0; j<uk.edges.size(); j++){
                sorter.insert(edges[uk.edges[j]].end_time);
            }
            std::set<int>::iterator it;
            int saved_te = *sorter.begin();
            for (it = sorter.begin(); it!= sorter.end(); it++){
                int te = *it;
                if (te == saved_te && it != sorter.begin()) continue;
                q.push({vertex_order[i], te, te});
                saved_te = te;

                while (!q.empty()) {
                    auto w = q.top(); q.pop();
                    if (w.node != uk.id) {
                        if (TBP_query(index, graph.upper[w.node], uk, w.start_time, w.end_time)) continue;
                        index.L_out[w.node].push_back({uk.id, w.start_time, w.end_time});
                    }

                    auto &w_edges = graph.upper[w.node].edges;
                    for (int k = 0 ; k < w_edges.size() ; k++) {
                        if (!(!visited1[w_edges[k]] && edges[w_edges[k]].end_time <= w.start_time)) continue;
                        visited1[w_edges[k]] = true;
                        algo::partition const *partition = find_partition(graph.lower[edges[w_edges[k]].dst].partitions, w_edges[k]);
                        if (!partition) continue;
                        for (int l = 0; l < partition->size() ; l++) {
                            auto &edge_id = (*partition)[l];
                            if (edge_id == w_edges[k]) continue;
                            if (graph.upper[edges[edge_id].src].order <= i) continue;
                            if (edges[edge_id].start_time <= maxT[edges[edge_id].src]) continue;
                            if (edge_colliding(w_edges[k], edge_id)) {
                                q.push({edges[edge_id].src, edges[edge_id].start_time, w.end_time});
                                maxT[edges[edge_id].src] = edges[edge_id].start_time;
                            }
                        }
                    }
                }
            }
        }
    }
    return index;
}

/**
 * The function "find_partition" searches for a specific edge ID within a vector of partitions and
 * returns a pointer to the partition that contains the edge ID, or nullptr if the edge ID is not
 * found.
 * 
 * @param partitions A vector of algo::partition objects, representing a collection of partitions.
 * @param edge_id The `edge_id` parameter is an integer representing the ID of an edge.
 * 
 * @return a pointer to a constant `algo::partition` object.
 */
algo::partition const *algo::find_partition(std::vector<algo::partition> const &partitions, int edge_id)
{
    for (int i = 0; i < partitions.size() ; i++)
        for (int j = 0; j < partitions[i].size() ; j++)
            if (partitions[i][j] == edge_id)
                return &(partitions[i]);
    return nullptr;
}

/**
 * The TBP_query function checks if there is a valid path between two nodes within a specified time
 * range.
 * 
 * @param index The parameter "index" is of type TBP_index const&, which is a constant reference to an
 * object of type TBP_index.
 * @param u The parameter `u` is of type `Node` and represents a node in the graph.
 * @param w The parameter "w" in the TBP_query function represents a node in the graph.
 * @param s The parameter "s" represents the start time for the query. It is used to filter the edges
 * based on their start time.
 * @param e The parameter "e" in the TBP_query function represents the end time for the query. It is
 * used to check if there is a valid path between nodes u and w within the specified time range (from s
 * to e).
 * 
 * @return a boolean value.
 */
bool algo::Core::TBP_query(TBP_index const &index, Node const &u, Node const &w, int s, int e)
{
    if (index.L_in[w.id].empty()) {
        for (int i = 0; i < index.L_out[u.id].size() ; i++) {
            if (index.L_out[u.id][i].node == w.id && index.L_out[u.id][i].start_time >= s && index.L_out[u.id][i].end_time <= e)
                return true;
        }
        return false;
    } else if (index.L_out[u.id].empty()) {
        for (int i = 0; i < index.L_in[w.id].size() ; i++) {
            if (index.L_in[w.id][i].node == u.id && index.L_in[w.id][i].start_time >= s && index.L_in[w.id][i].end_time <= e)
                return true;
        }
        return false;
    }

    int i = 0, j = 0;
    while ((i + j + 2) <= (index.L_out[u.id].size() + index.L_in[w.id].size())) {
        auto &a = index.L_out[u.id][i];
        auto &b = index.L_in[w.id][j];
        if (a.node == w.id && a.start_time >= s && a.end_time <= e) return true;
        if (b.node == u.id && b.start_time >= s && b.end_time <= e) return true;
        if (a.node == b.node && a.start_time >= s && a.end_time <= b.start_time && b.end_time <= e) return true;

        if (i == index.L_out[u.id].size() - 1)
            j++;
        else if (j == index.L_in[w.id].size() - 1)
            i++;
        else {
            if (graph.upper[a.node].order < graph.upper[b.node].order)
                i++;
            else if (graph.upper[a.node].order > graph.upper[b.node].order) {
                j++;
            } else {
                i++; j++;
            }
        }
    }
    return false;
}

/**
 * The function checks if two edges are colliding based on their start and end times.
 * 
 * @param a The parameter "a" represents the index of the first edge in the edges array.
 * @param b The parameter "b" in the edge_colliding function represents the index of the second edge in
 * the edges array.
 * 
 * @return a boolean value.
 */
bool algo::Core::edge_colliding(int a, int b)
{
    return (std::min(edges[a].end_time, edges[b].end_time) > std::max(edges[a].start_time, edges[b].start_time));
}

/**
 * The function `inverted_in_label_set` takes a 2D vector `L_in` and returns a new 2D vector `invL_in`
 * where the entries are inverted based on the node value.
 * 
 * @param L_in L_in is a vector of vectors of algo::entry objects. Each inner vector represents a set
 * of entries, and each entry contains information about a node, start time, and end time.
 * 
 * @return a vector of vectors of `algo::entry` objects.
 */
std::vector<std::vector<algo::entry>> algo::Core::inverted_in_label_set(std::vector<std::vector<algo::entry>> L_in)
{
    std::vector<std::vector<algo::entry>> invL_in(L_in.size());
    for (int i = 0; i < L_in.size(); i++) {
        for (int j = 0; j < L_in[i].size(); j++) {
            invL_in[L_in[i][j].node].push_back({i, L_in[i][j].start_time, L_in[i][j].end_time});
        }
    }
    return invL_in;
}

/**
 * The function SSRQ calculates the number of nodes reachable from a given node in a directed graph.
 * 
 * @param node The "node" parameter is a reference to an object of type "Node". It represents a node in
 * a graph.
 * @param L_out L_out is a vector of vectors of algo::entry objects. It represents the outgoing edges
 * from each node in the graph. Each inner vector contains the entries for a specific node, where each
 * entry represents an outgoing edge from that node.
 * @param invL_in The parameter `invL_in` is a reference to a vector of vectors of `algo::entry`
 * objects. It represents the inverse adjacency list of a graph, where each inner vector contains the
 * entries for a particular node. Each `entry` object represents an edge in the graph and contains
 * information such
 * 
 * @return the size of the vector "reach".
 */
int algo::Core::SSRQ(Node &node, std::vector<std::vector<algo::entry>> &L_out, std::vector<std::vector<algo::entry>> &invL_in)
{
    std::vector<bool> list_reach(graph.upper.size(), false);
    std::vector<int> reach;
    for (auto entry : L_out[node.id]) {
        int idu = entry.node;
        if (!list_reach[idu]) {
            reach.emplace_back(idu);
            list_reach[idu] = true;
        }
        for (auto entry1 : invL_in[idu]) {
            if (entry.end_time <= entry1.start_time && entry1.node != node.id) {
                if (!list_reach[entry1.node]) {
                    reach.emplace_back(entry1.node);
                    list_reach[entry1.node] = true;
                }
            }
        }
    }
    for (auto entry : invL_in[node.id]) {
        int idu = entry.node;
        if (!list_reach[idu]) {
            reach.emplace_back(idu);
            list_reach[idu] = true;
        }
    }
    return reach.size();
}

/**
 * The function maxR_baseline calculates the maximum reachability value for a given graph using the
 * SSRQ algorithm.
 * 
 * @param L_out L_out is a reference to a vector of vectors of algo::entry objects. It is an output
 * parameter that will store the reachability information for each node in the graph.
 * @param invL_in The parameter `invL_in` is a reference to a vector of vectors of `algo::entry`
 * objects. It is used as an input to the `maxR_baseline` function.
 * 
 * @return an integer value, which represents the maximum reach value.
 */
int algo::Core::maxR_baseline(std::vector<std::vector<algo::entry>> &L_out, std::vector<std::vector<algo::entry>> &invL_in)
{
    int maxR = 0;
    for (int i = 0; i < graph.upper.size(); i++) {
        int reach = SSRQ(graph.upper[i], L_out, invL_in);
        maxR = std::max(maxR, reach);
    }
    return maxR;
}

/**
 * The function `maxReach` takes a graph as input and returns the maximum reach and a list of reach
 * values for each vertex in the graph.
 * 
 * @param graph The "graph" parameter is an object of type algo::Graph, which represents a graph data
 * structure. It is passed by reference to the maxReach function.
 * 
 * @return a tuple containing two values: an integer representing the maximum reach and a vector of
 * integers representing the reach for each vertex in the graph.
 */
std::tuple<int, std::vector<int>> algo::Core::maxReach(algo::Graph &graph)
{
    int maxReach = 0;
    std::vector<int> list_reach;
    auto order = get_vertex_order(graph);
    for (int i = 0 ; i < order.size() ; i++) {
        graph.upper[order[i]].order = i;
    }
    algo::TBP_index index = TBP_build(order, graph);
    std::vector<std::vector<algo::entry>> invL_in = inverted_in_label_set(index.L_in);
    for (auto u : graph.upper) {
        int reach = SSRQ(u, index.L_out, invL_in);
        list_reach.push_back(reach);
        maxReach = std::max(reach, maxReach);
    }
    return std::make_tuple(maxReach, list_reach);
}

/**
 * The function removes a vertex from a graph by deleting all edges connected to that vertex.
 * 
 * @param node The parameter "node" is a reference to an object of type "algo::Node".
 * @param graph The "graph" parameter is a reference to an object of type "algo::Graph". It represents
 * the graph from which the vertex needs to be removed.
 * 
 * @return a graph object.
 */
algo::Graph algo::Core::remove_vertex(algo::Node &node, algo::Graph &graph)
{
    for (auto edge1 : node.edges){
        algo::Node &node1 = graph.upper[edges[edge1].src];
        auto position = std::find(node1.edges.begin(), node1.edges.end(), edge1);
        if (position != node1.edges.end()) {
            node1.edges.erase(position);
        }
    }
    node.edges.clear();
    return graph;
}

/**
 * The function "choose_vertex" selects a vertex from a graph that minimizes the maximum reachability
 * of the remaining graph after removing the selected vertex.
 * 
 * @param graph The "graph" parameter is an object of type "algo::Graph", which represents a graph data
 * structure. It contains a vector called "lower" which represents the lower level of the graph. Each
 * element in the "lower" vector is an object of type "algo::Vertex", which represents a
 * 
 * @return an integer value, which represents the chosen vertex.
 */
int algo::Core::choose_vertex(algo::Graph &graph)
{
    int minimalMaxReach = INT_MAX;
    int vertex = 0;
    for (int i = 0; i < graph.lower.size(); i++){
        algo::Graph G1 = graph;
        if (!G1.lower[i].edges.empty()){
            G1 = remove_vertex(G1.lower[i],G1);
            std::tuple<int, std::vector<int>> reach = maxReach(G1);
            
            if (std::get<0>(reach) < minimalMaxReach){
                minimalMaxReach = std::get<0>(reach);
                vertex = i;
            }
        }
    }
    return vertex;
}

/**
 * The function `greedy1` implements a greedy algorithm to remove `k` vertices from a graph, returning
 * the list of removed vertices.
 * 
 * @param graph The "graph" parameter is an object of type algo::Graph, which represents a graph data
 * structure. It contains information about the vertices and edges of the graph.
 * @param k The parameter "k" represents the number of vertices to be removed from the graph.
 * 
 * @return The function `greedy1` returns a vector of integers, which represents the vertices that have
 * been removed from the graph.
 */
std::vector<int> algo::Core::greedy1(algo::Graph &graph, int k)
{
    std::vector<int> vertices_removed;
    algo::Graph G1 = graph;
    for(int i = 0; i < k; i++){
        int idv = choose_vertex(G1);
        vertices_removed.push_back(idv);
        G1 = remove_vertex(G1.lower[idv], G1);
    }
    return vertices_removed;
}

/**
 * The function `sort_te` is used to insert a `wedge` object into a sorted vector `list_w` based on the
 * `end_time` property of the `wedge` object.
 * 
 * @param list_w A reference to a vector of wedge objects. This vector represents a list of wedges that
 * need to be sorted.
 * @param wedge The `wedge` parameter is a reference to a `wedge` object.
 */
void algo::Core::sort_te(std::vector<wedge> &list_w, wedge &wedge)
{
    int low = 0;
    int high = list_w.size();
    while (low < high) {
        int mid = (low + high) / 2;
        if (wedge.end_time < list_w[mid].end_time) {
            high = mid;
        } else {
            low = mid + 1;
        }
    }
    list_w.insert(list_w.begin()+low, wedge);
}

/**
 * The function checks if one wedge is a dominant wedge of another wedge.
 * 
 * @param wedge1 The first wedge object, which represents a time interval with a start time and an end
 * time.
 * @param wedge2 The `wedge2` parameter is an object of type `wedge`.
 * 
 * @return A boolean value is being returned.
 */
bool algo::Core::dominant_wedge(algo::wedge &wedge1, algo::wedge &wedge2)
{
    if (wedge2.start_time <= wedge1.start_time && wedge2.end_time >= wedge1.end_time)
        return true;
    return false;
}


/**
 * The function `necessary_wedges` takes a list of wedges and a single wedge as input, and returns a
 * new list of wedges that are necessary based on certain conditions.
 * 
 * @param list_w A reference to a vector of algo::wedge objects.
 * @param wedge The parameter "wedge" is of type "algo::wedge", which is a user-defined data type. It
 * represents a wedge object.
 * 
 * @return a vector of algo::wedge objects.
 */
std::vector<algo::wedge> algo::Core::necessary_wedges(std::vector<algo::wedge> &list_w, algo::wedge &wedge)
{
    std::vector<algo::wedge> new_list_w = list_w;
    for (int i = 0; i < list_w.size(); i++) {
        if (dominant_wedge(list_w[i], wedge)) {
            return list_w;
        }
        if (dominant_wedge(wedge, list_w[i])) {
            int pos;
            for (int j = 0; j < new_list_w.size(); j++) {
                if (new_list_w[j].start_time == list_w[i].start_time && new_list_w[j].end_time == list_w[i].end_time) {
                    pos = j;
                    break;
                }
            }
            new_list_w.erase(new_list_w.begin() + pos);
        }
    }
    sort_te(new_list_w, wedge);
    return new_list_w;
}

/**
 * The function `WTB_index_build` builds a WTB index for a graph by iterating through the
 * lower level nodes and their edges, and creating a map of vertex pairs to a vector of wedges.
 * 
 * @param graph The `graph` parameter is an object of type `algo::Graph`, which represents a graph data
 * structure. It contains a vector `lower` which stores nodes, and each node contains a vector `edges`
 * which stores edges.
 * @param start The "start" parameter is an integer representing the starting time for filtering the
 * edges. It is used to check if the end time of an edge is within the specified range.
 * @param end The "end" parameter is an integer value representing the end time for filtering edges. It
 * is used to check if the end time of an edge falls within a specified range.
 * 
 * @return a vector of maps. Each map contains a pair of vertices as the key and a vector of wedges as
 * the value.
 */
std::vector<std::map<algo::vertices_pair, std::vector<algo::wedge>, algo::pairCompare>> algo::Core::WTB_index_build(algo::Graph &graph, int const &start, int const &end)
{
    std::vector<std::map<algo::vertices_pair, std::vector<algo::wedge>, algo::pairCompare>> list_dict(graph.lower.size());
    for (int i = 0; i < graph.lower.size(); i++) {
        algo::Node v = graph.lower[i];
        for (int j = 0; j < v.edges.size(); j++) {
            algo::edge edge1 = edges[v.edges[j]];
            int src = edge1.src;
            for (int k = j+1; k < v.edges.size(); k++) {
                algo::edge edge2 = edges[v.edges[k]];
                int dst = edge2.src;
                if (src != dst && edge_colliding(v.edges[j], v.edges[k])) {
                    if (edge2.end_time >= start && edge2.end_time <= end) {
                        algo::wedge w = {edge1.start_time, edge2.end_time};
                        std::vector<algo::wedge> list_w;
                        std::map<algo::vertices_pair, std::vector<algo::wedge>, algo::pairCompare>::iterator it = list_dict[i].find({src,dst});
                        if (it == list_dict[i].end()) {
                            list_dict[i][{src,dst}] = necessary_wedges(list_w, w);
                        } else {
                            list_w = it->second;
                            list_dict[i][{src,dst}] = necessary_wedges(list_w, w);
                        }
                    }
                    if (edge1.end_time >= start && edge1.end_time <= end) {
                        algo::wedge w1 = {edge2.start_time, edge1.end_time};
                        std::vector<algo::wedge> list_w1;
                        std::map<algo::vertices_pair, std::vector<algo::wedge>, algo::pairCompare>::iterator it1 = list_dict[i].find({dst,src});
                        if (it1 != list_dict[i].end()) {
                            list_w1 = it1->second;
                        }
                        list_dict[i][{dst,src}] = necessary_wedges(list_w1, w1);
                    }
                }
            }
        }
    }
    return list_dict;
}

/**
 * The function `sort_list_entries` sorts the WTB-index based on end times.
 * 
 * @param vector The parameter `list_dict` is a vector of maps. Each map in the vector represents a
 * dictionary where the key is of type `algo::vertices_pair` and the value is a vector of `algo::wedge`
 * objects.
 * 
 * @return a sorted vector of `algo::third_entry` objects.
 */
std::vector<algo::third_entry> algo::Core::sort_list_entries(std::vector<std::map<algo::vertices_pair, std::vector<algo::wedge>, algo::pairCompare>> &list_dict)
{
    std::vector<algo::third_entry> list_entries_ordered;
    for (int i = 0; i < list_dict.size(); i++) {
        if (!list_dict[i].empty()) {
            for (auto it = list_dict[i].begin(); it != list_dict[i].end(); it++) {
                algo::vertices_pair key = it->first;
                std::vector<algo::wedge> val = it->second;
                for (int j = 0; j < val.size() ; j++) {
                    list_entries_ordered.push_back({key.src, key.dst, val[j].start_time, val[j].end_time, i});
                }
            }
        }
    }
    std::sort(list_entries_ordered.begin(), list_entries_ordered.end(),[](const algo::third_entry &a, const algo::third_entry &b)
    {
        return a.end_time < b.end_time;
    });
    return list_entries_ordered;
}

/**
 * The function `choose_vertex_WTB` first constructs the lists L_reach and L_count for a WTB-index
 * and then finds the best vertex to remove using L_count.
 * 
 * @param vector `list_dict` is a WTB-index (sorted).
 * @param size The parameter `size` is the number of upper vertices of the graph.
 * 
 * @return an integer value, which represents the ID of the chosen vertex.
 */
int algo::Core::choose_vertex_WTB(std::vector<std::map<algo::vertices_pair, std::vector<algo::wedge>, algo::pairCompare>> &list_dict, int &size)
{
    std::vector<algo::third_entry> list_entries_ordered = sort_list_entries(list_dict);
    std::vector<int> list_max_reach;
    for (int i = 0; i < list_dict.size(); i++) {
        if (list_dict[i].empty())
            list_max_reach.push_back(INT_MAX);
        else {
            list_max_reach.push_back(0);
        }
    }
    std::vector<robin_hood::unordered_flat_map<int, algo::list_removal>> list_reach(size);
    std::vector<algo::count_removal> count_reach(size);
    std::vector<int> count_wedges(list_dict.size(), 0);

    for (auto entry : list_entries_ordered) {
        count_wedges[entry.vertex_id]++;
        if (list_reach[entry.dst].count(entry.src) == 0) {
            list_reach[entry.dst][entry.src] = {entry.end_time, {{entry.vertex_id, INT_MAX}}};
            count_reach[entry.src].reach++;
            count_reach[entry.src].removal[entry.vertex_id]++;
        } else {
            for (auto it = list_reach[entry.dst][entry.src].removal.begin(); it != list_reach[entry.dst][entry.src].removal.end(); it++) {
                if (it->first != entry.vertex_id && it->second == INT_MAX) {
                    it->second = entry.end_time;
                    count_reach[entry.src].removal[it->first]--;
                }
            }
        }
        for (auto it = list_reach[entry.src].begin(); it != list_reach[entry.src].end(); it++) {
            if (it->second.reach > entry.start_time || it->first == entry.dst) {
                continue;
            } else {
                if (list_reach[entry.dst].count(it->first) == 0) {
                    robin_hood::unordered_flat_map<int, int> list_v = list_reach[entry.src][it->first].removal;
                    for (auto it1 = list_reach[entry.src][it->first].removal.begin(); it1 != list_reach[entry.src][it->first].removal.end(); it1++) {
                        if (it1->second <= entry.start_time) list_v.erase(it1->first);
                        else {
                            list_v[it1->first] = INT_MAX;
                            count_reach[it->first].removal[it1->first]++;
                        }
                    }
                    list_v[entry.vertex_id] = INT_MAX;
                    if (it->second.removal.count(entry.vertex_id) == 0 
                    || it->second.removal[entry.vertex_id] <= entry.start_time) {
                        count_reach[it->first].removal[entry.vertex_id]++;
                    }
                    list_reach[entry.dst][it->first] = {entry.end_time, list_v};
                    list_v.clear();
                    count_reach[it->first].reach++;
                } else {
                    for (auto it2 = list_reach[entry.dst][it->first].removal.begin(); it2 != list_reach[entry.dst][it->first].removal.end(); it2++) {
                        if (it2->first != entry.vertex_id &&  it2->second == INT_MAX) {
                            if (list_reach[entry.src][it->first].removal.count(it2->first) == 0 || list_reach[entry.src][it->first].removal[it2->first] <= entry.start_time) {
                                it2->second = entry.end_time;
                                count_reach[it->first].removal[it2->first]--;
                            }
                        }
                    }
                }
            }
        }
    }
    for (int i = 0; i < count_reach.size(); i++) {
        if (count_reach[i].reach == 0 || count_reach[i].removal.empty()) continue;
        robin_hood::unordered_flat_map<int, short> l_v = count_reach[i].removal;
        for (int j = 0; j < list_max_reach.size(); j++) {
            int reach = count_reach[i].reach;
            if (l_v.count(j) != 0)
                reach -= l_v[j];
            if (reach > list_max_reach[j])
                list_max_reach[j] = reach;
        }
    }
    int idv = 0;
    int minMaxReach = INT_MAX;
    int curr_wedges = 0;
    for (int i = 0; i < list_max_reach.size(); i++) {
        if (list_max_reach[i] < minMaxReach || (list_max_reach[i] == minMaxReach && count_wedges[i] > curr_wedges)) {
            idv = i;
            minMaxReach = list_max_reach[i];
            curr_wedges = count_wedges[i];
        }
    }
    return idv;
}

/**
 * The function `greedy_WTB` implements a greedy algorithm to remove `k` vertices from a graph, returning
 * the list of removed vertices.
 * 
 * @param graph Temporal bipartite graph.
 * @param k The parameter "k" represents the number of vertices to be removed from the graph.
 * @param start
 * @param end
 * 
 * @return Vector of the removed vertices.
 */
std::vector<int> algo::Core::greedy_WTB(algo::Graph &graph, int k, int const &start, int const &end)
{
    std::vector<int> vertices_removed;
    std::vector<std::map<algo::vertices_pair, std::vector<algo::wedge>, algo::pairCompare>> list_dict = WTB_index_build(graph, start, end);
    int size = graph.upper.size();
    for(int i = 0; i < k; i++){
        int idv = choose_vertex_WTB(list_dict, size);
        vertices_removed.push_back(idv);
        list_dict[idv].clear();
    }
    return vertices_removed;
}

/**
 * The function "list_position" takes a vector of "third_entry" objects and returns two vectors of
 * vectors, where the first vector contains the positions of each source node in the input vector, and
 * the second vector contains the positions of each destination node in the input vector.
 * 
 * @param LEO A WTB-index (sorted).
 * 
 * @return The function `list_position` returns a `std::tuple` containing two elements:
 * 1. `src_positions`: A `std::vector` of `std::vector<int>`, where each inner vector represents the
 * positions of the source nodes in the `LEO` vector.
 * 2. `dst_positions`: A `std::vector` of `std::vector<int>`, where each inner vector
 */
std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>> algo::Core::list_position(std::vector<algo::third_entry> &LEO)
{
    std::vector<std::vector<int>> src_positions(graph.upper.size());
    std::vector<std::vector<int>> dst_positions(graph.upper.size());
    for (int i = 0; i < LEO.size(); i++) {
        src_positions[LEO[i].src].emplace_back(i);
        dst_positions[LEO[i].dst].emplace_back(i);
    }
    return std::make_tuple(src_positions, dst_positions);
}

/**
 * The function SPRQ_static checks if there is a path between two nodes in a graph based on the given
 * source and destination positions.
 * 
 * @param LEO A vector of algo::third_entry objects, representing a list of entries.
 * @param src_pos src_pos is a vector of vectors of integers that represents the positions of source
 * nodes in the LEO vector. Each inner vector corresponds to a source node and contains the positions
 * of its entries in the LEO vector.
 * @param dst_pos dst_pos is a vector of vectors that represents the positions of destinations. Each
 * inner vector corresponds to a destination and contains the positions of that destination in the LEO
 * vector.
 * @param idu The parameter "idu" represents the source node ID in the graph.
 * @param idw The parameter "idw" represents the ID of the destination node in the graph.
 * 
 * @return a boolean value.
 */
bool algo::Core::SPRQ_static(std::vector<algo::third_entry> &LEO, std::vector<std::vector<int>> &src_pos, std::vector<std::vector<int>> &dst_pos, int &idu, int &idw)
{
    int first_pos;
    int last_pos;
    if (dst_pos[idw].empty()) {
        return false;
    } else {
        last_pos = dst_pos[idw].back();
    }
    if (src_pos[idu].empty()) {
        return false;
    } else {
        first_pos = src_pos[idu][0];
    }
    std::vector<int> list_reach(graph.upper.size());
    for (int i = first_pos; i < last_pos+1; i++) {
        algo::third_entry entry = LEO[i];
        if (entry.src == idu) {
            if (entry.dst == idw) return true;
            if (list_reach[entry.dst] == 0) {
                list_reach[entry.dst] = entry.end_time;
            }
        } else {
            if (list_reach[entry.src] != 0) {
                if (list_reach[entry.src] > entry.start_time) continue;
                if (entry.dst == idw) return true;
                if (idu != entry.dst && list_reach[entry.dst] == 0) {
                    list_reach[entry.dst] = entry.end_time;
                }
            }
        }
    }
    return false;
}

/**
 * The function `SSRQ_static` takes in a list of entries, a list of positions, and an ID, and returns
 * a list of reachable nodes.
 * 
 * @param LEO A vector of algo::third_entry objects, representing a list of entries.
 * @param src_pos src_pos is a reference to a vector of vectors of integers. It represents the
 * positions of source nodes in the LEO vector. Each inner vector corresponds to a source node and
 * contains the positions of its entries in the LEO vector.
 * @param idu The parameter "idu" is an integer representing the source node ID.
 * 
 * @return a vector of integers.
 */
std::vector<int> algo::Core::SSRQ_static(std::vector<algo::third_entry> &LEO, std::vector<std::vector<int>> &src_pos, int &idu)
{
    int first_pos;
    if (src_pos[idu].empty()) {
        return {};
    } else {
        first_pos = src_pos[idu][0];
    }
    std::vector<int> list_reach(graph.upper.size());
    std::vector<int> reach;
    for (int i = first_pos; i < LEO.size(); i++) {
        algo::third_entry entry = LEO[i];
        if (entry.src == idu) {
            if (list_reach[entry.dst] == 0) {
                list_reach[entry.dst] = entry.end_time;
                reach.emplace_back(entry.dst);
            }
        } else {
            if (list_reach[entry.src] != 0) {
                if (list_reach[entry.src] > entry.start_time) continue;
                if (idu != entry.dst && list_reach[entry.dst] == 0) {
                    list_reach[entry.dst] = entry.end_time;
                    reach.emplace_back(entry.dst);
                }
            }
        }
    }
    return reach;
}

/**
 * The function `max_reach_static` calculates the maximum reachability of nodes in a graph.
 * 
 * @param LEO A WTB-index (sorted).
 * 
 * @return the maximum reachability value, which is an integer.
 */
int algo::Core::max_reach_static(std::vector<algo::third_entry> &LEO, std::vector<std::vector<int>> &src_pos)
{
    //int epoc = 0;
    int maxR = 0;
    for (int i = 0; i < graph.upper.size(); i++) {
        //if (epoc % 10000 == 0 && epoc != 0) std::cout << "Ep = " << epoc << std::endl;
        //epoc++;
        int reach = SSRQ_static(LEO, src_pos, i).size();
        maxR = std::max(maxR, reach);
    }
    return maxR;
}

/**
 * The function `vertex_removal` removes all wedges associated to the specified vertex.
 * 
 * @param LEO A WTB-index (sorted).
 * @param vertex_id The `vertex_id` parameter is an integer that represents the ID of the vertex that
 * needs to be removed from the vector `LEO`.
 * @param layer The "layer" parameter is a boolean flag that determines whether the vertex is in the
 * upper or lower layer.
 * 
 * @return the modified vector `LEO` after removing entries.
 */
std::vector<algo::third_entry> algo::Core::vertex_removal(std::vector<algo::third_entry> &LEO, int vertex_id, bool layer)
{
    if (layer) {
        for (auto it = LEO.begin(); it != LEO.end();) {
            if (it->src == vertex_id || it->dst == vertex_id) {
                it = LEO.erase(it);
            } else {
                it++;
            }
        }
    } else {
        for (auto it = LEO.begin(); it != LEO.end();) {
            if (it->vertex_id == vertex_id) {
                it = LEO.erase(it);
            } else {
                it++;
            }
        }
    }
    return LEO;
}

/**
 * The function `edge_removal` removes wedges associated to a specific edge.
 * 
 * @param LEO A WTB-index (sorted).
 * @param edge The specific edge.
 * 
 * @return the modified vector `LEO` after removing entries.
 */
std::vector<algo::third_entry> algo::Core::edge_removal(std::vector<algo::third_entry> &LEO, algo::edge &edge)
{
    for (auto it = LEO.begin(); it != LEO.end();) {
        if (it->src == edge.src && it->vertex_id == edge.dst && it->start_time == edge.start_time) {
            it = LEO.erase(it);
        }
        else if (it->dst == edge.src && it->vertex_id == edge.dst && it->end_time == edge.end_time) {
            it = LEO.erase(it);
        } else {
            it++;
        }
    }
    return LEO;
}

/**
 * The function checks if two edges are colliding based on their start and end times.
 * 
 * @param e1
 * @param e2
 * 
 * @return a boolean value.
 */
bool algo::Core::edge_colliding2(algo::edge const &e1, algo::edge const &e2)
{
    return (std::min(e1.end_time, e2.end_time) > std::max(e1.start_time, e2.start_time));
}

/**
 * The function "entry_insertion" inserts a given entry into a sorted vector of entries in a way that
 * maintains the sorted order.
 * 
 * @param LEO.
 * @param entry
 */
void algo::Core::entry_insertion(std::vector<algo::third_entry> &LEO, algo::third_entry &entry)
{
    int low = 0;
    int high = LEO.size();
    while (low < high) {
        int mid = (low + high) / 2;
        if (entry.end_time < LEO[mid].end_time) {
            high = mid;
        } else {
            low = mid + 1;
        }
    }
    LEO.insert(LEO.begin()+low, entry);
}

/**
 * The function `edge_addition` add all newly created wedges to the index.
 * 
 * @param LEO.
 * @param edge
 * 
 * @return
 */
std::vector<algo::third_entry> algo::Core::edge_addition(std::vector<algo::third_entry> &LEO, algo::edge &edge)
{
    for (auto ide : graph.lower[edge.dst].edges) {
        algo::edge edge1 = edges[ide];
        if (edge.src != edge1.src && edge_colliding2(edge, edge1)) {
            algo::third_entry entry = {edge.src, edge1.src, edge.start_time, edge1.end_time, edge.dst};
            entry_insertion(LEO, entry);
            algo::third_entry entry1 = {edge1.src, edge.src, edge1.start_time, edge.end_time, edge.dst};
            entry_insertion(LEO, entry1);
        }
    }
    return LEO;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * The following functions are the test functions for the different operations (SPRQ, SSRQ, MRQ...) for both
 * the TBP-index and the WTB-index.
 */
int algo::Core::test_construction_baseline(algo::Graph graph, algo::partition order)
{
    auto start = std::chrono::high_resolution_clock::now();

    auto index = TBP_build(order, graph);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "TBP - Construction: " << duration.count() << " microseconds" << std::endl;
    return 0;
}

int algo::Core::test_construction_WTB(algo::Graph graph)
{
    auto start = std::chrono::high_resolution_clock::now();

    auto index = WTB_index_build(graph, 0, INT_MAX);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "WTB - Construction: " << duration.count() << " microseconds" << std::endl;
    return 0;
}

int algo::Core::test_SPRQ_baseline(algo::TBP_index const &index, std::vector<int> test, std::vector<int> test2)
{
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<bool> eval;
    for (int i = 0; i < test.size(); i++) {
        bool e = TBP_query(index, graph.upper[test[i]], graph.upper[test2[i]], 0, INT_MAX);
        eval.push_back(e);
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "TBP - SPRQ for " << test.size() << " queries: " << duration.count() << " microseconds" << std::endl;
    return 0;
}

int algo::Core::test_SPRQ_WTB(std::vector<algo::third_entry> &LEO, std::vector<std::vector<int>> &src_pos, std::vector<std::vector<int>> &dst_pos, std::vector<int> test, std::vector<int> test2)
{
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<bool> eval;
    for (int i = 0; i < test.size(); i++) {
        bool e = SPRQ_static(LEO, src_pos, dst_pos, test[i], test2[i]);
        eval.push_back(e);
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "WTB - SPRQ for " << test.size() << " queries: " << duration.count() << " microseconds" << std::endl;
    return 0;
}

int algo::Core::test_SSRQ_baseline(std::vector<std::vector<algo::entry>> &L_out, std::vector<std::vector<algo::entry>> &invL_in, std::vector<int> test)
{
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<int> eval;
    for (int i = 0; i < test.size(); i++) {
        int e = SSRQ(graph.upper[test[i]], L_out, invL_in);
        eval.emplace_back(e);
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "TBP - SSRQ for " << test.size() << " queries: " << duration.count() << " microseconds" << std::endl;
    return 0;
}

int algo::Core::test_SSRQ_WTB(std::vector<algo::third_entry> &LEO, std::vector<std::vector<int>> &src_pos, std::vector<int> test)
{
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<int> eval;
    for (int i = 0; i < test.size(); i++) {
        int e = SSRQ_static(LEO, src_pos, test[i]).size();
        eval.emplace_back(e);
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "WTB - SSRQ for " << test.size() << " queries: " << duration.count() << " microseconds" << std::endl;
    return 0;
}

int algo::Core::test_maxR_baseline(std::vector<std::vector<algo::entry>> &L_out, std::vector<std::vector<algo::entry>> &invL_in)
{
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<int> eval;
    for (int i = 0; i < 1; i++) {
        int e = maxR_baseline(L_out, invL_in);
        eval.emplace_back(e);
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "TBP - MaxReach: " << duration.count() << " microseconds" << std::endl;
    return 0;
}

int algo::Core::test_maxR_WTB(std::vector<algo::third_entry> &LEO, std::vector<std::vector<int>> &src_pos)
{
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<int> eval;
    for (int i = 0; i < 1; i++) {
        int e = max_reach_static(LEO, src_pos);
        eval.emplace_back(e);
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "WTB - MaxReach: " << duration.count() << " microseconds" << std::endl;
    return 0;
}

int algo::Core::test_VRCQ_baseline(algo::Graph &graph)
{
    auto start = std::chrono::high_resolution_clock::now();

    int minimalMaxReach = INT_MAX;
    int vertex = 0;
    for (int i = 0; i < graph.lower.size(); i++){
        algo::Graph G1 = graph;
        if (!G1.lower[i].edges.empty()){
            G1 = remove_vertex(G1.lower[i],G1);
            std::tuple<int, std::vector<int>> reach = maxReach(G1);
            if (std::get<0>(reach) < minimalMaxReach){
                minimalMaxReach = std::get<0>(reach);
                vertex = i;
            }
        }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "TBP - VRCQ: " << duration.count() << " microseconds" << std::endl;
    return vertex;
}

int algo::Core::test_VRCQ_WTB(std::vector<std::map<algo::vertices_pair, std::vector<algo::wedge>, algo::pairCompare>> &list_dict, int size)
{
    auto start = std::chrono::high_resolution_clock::now();

    int id = choose_vertex_WTB(list_dict, size);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "WTB - VRCQ: " << duration.count() << " microseconds" << std::endl;
    return 0;
}

int algo::Core::test_vertex_removal(std::vector<algo::third_entry> &LEO, std::vector<int> test, std::vector<int> test2)
{
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<algo::third_entry> new_LEO = LEO;
    for (int i = 0; i < test.size()/2; i++){
        vertex_removal(new_LEO, test[i], 0);
    }
    for (int i = 0; i < test2.size()/2; i++){
        vertex_removal(new_LEO, test2[i], 1);
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "WTB - " << test.size()/2 + test2.size()/2 << " vertex removals: " << duration.count() << " microseconds" << std::endl;
    return 0;
}

int algo::Core::test_edge_removal(std::vector<algo::third_entry> &LEO, std::vector<int> test_edges)
{
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<algo::third_entry> new_LEO = LEO;
    for (int i = 0; i < test_edges.size(); i++){
        edge_removal(new_LEO, edges[test_edges[i]]);
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "WTB - " << test_edges.size() << " edge removals: " << duration.count() << " microseconds" << std::endl;
    return 0;
}

int algo::Core::test_edge_addition(std::vector<algo::third_entry> &LEO, std::vector<algo::edge> test_new_edges)
{
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<algo::third_entry> new_LEO = LEO;
    for (int i = 0; i < test_new_edges.size(); i++){
        edge_addition(new_LEO, test_new_edges[i]);
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "WTB - " << test_new_edges.size() << " edge additions: " << duration.count() << " microseconds" << std::endl;
    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * Main function to time the operations
*/

int algo::Core::main(int ac, char **av)
{
    if (ac != 2)
        throw std::runtime_error("Need to specify dataset src");

    init_data(av[1]);
    
    std::vector<int> test;
    std::vector<int> test2;
    std::vector<int> test3;
    std::vector<int> test_edges;
    std::vector<algo::edge> test_new_edges;
    std::random_device eng;
    using dist_params = typename std::uniform_int_distribution<int>::param_type;
    std::uniform_int_distribution<int> dist (0, graph.upper.size()-1);
    std::uniform_int_distribution<int> dist2 (0, graph.lower.size()-1);
    std::uniform_int_distribution<int> dist_edges (0, edges.size()-1);
    std::uniform_int_distribution<int> dist_start_time (1100000000, 1400000000);

    for (int i = 0; i < 1000; i++) {
        int a = dist(eng);
        int b = dist(eng);
        int c = dist2(eng);
        int d = dist_edges(eng);
        while (b == a) {
            b = dist(eng);
        }
        test.emplace_back(a);
        test2.emplace_back(b);
        test3.emplace_back(c);
        test_edges.emplace_back(d);
        int start = dist_start_time(eng);
        int end = start + rndm(3600, 172800, -2.5);
        algo::edge new_edge = {a, c, start, end};
        test_new_edges.emplace_back(new_edge);
    }

    ///////////////////////////
    // TBP

    auto order = get_vertex_order(graph);
    for (int i = 0 ; i < order.size() ; i++) {
        graph.upper[order[i]].order = i;
    }

    test_construction_baseline(graph, order);

    auto index = TBP_build(order, graph);
    auto invL_in = inverted_in_label_set(index.L_in);

    test_SPRQ_baseline(index, test, test2);
    test_SSRQ_baseline(index.L_out, invL_in, test);
    test_maxR_baseline(index.L_out, invL_in);
    //test_VRCQ_baseline(graph);

    ////////////////////////////
    // WTB

    test_construction_WTB(graph);

    auto list_dict1 = WTB_index_build(graph, 0, INT_MAX);
    std::vector<algo::third_entry> LEO = sort_list_entries(list_dict1);
    std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>> posit = list_position(LEO);

    test_SPRQ_WTB(LEO, std::get<0>(posit), std::get<1>(posit), test, test2);
    test_SSRQ_WTB(LEO, std::get<0>(posit), test);
    test_maxR_WTB(LEO, std::get<0>(posit));
    test_VRCQ_WTB(list_dict1, graph.upper.size());

    test_vertex_removal(LEO, test, test3);
    test_edge_removal(LEO, test_edges);
    test_edge_addition(LEO, test_new_edges);

    return 0;
}