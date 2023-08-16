#ifndef CORE
#define CORE

#include <iostream>
#include <random>
#include <algorithm>
#include <vector>
#include <numeric>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <map>
#include <unordered_map>
#include "robin_hood.h"
#include <limits.h>

namespace algo {

    /**
     * The above code defines a struct called "edge" with four integer fields and a print function.
     * @property {int} src - The source vertex of the edge.
     * @property {int} start_time - The start time of the edge, indicating when the edge becomes active
     * or starts.
     */
    struct edge {
        int src, dst;
        int start_time, end_time;

        void print(void) {
            std::cout << src << " " << dst << " " << start_time << " " << end_time << std::endl;
        }
    };

    typedef std::vector<int> partition;

    /* The Node class represents a node in a graph, with an ID, order, edges, and partitions. */
    class Node {
    public:
        Node(int ID) {
            id = ID;
            order = 0;
        }

        int id;
        int order;
        std::vector<int> edges;
        std::vector<algo::partition> partitions;
    };

    /* The Graph class represents a graph with two sets of nodes, upper and lower. */
    class Graph {
    public:
        Graph(void) {}

        void init(int nb_src, int nb_dst){
            for (int i = 0; i <= nb_src; i++)
                upper.emplace_back(i);
            for (int i = 0; i <= nb_dst; i++)
                lower.emplace_back(i);
        }

        std::vector<Node> upper;
        std::vector<Node> lower;
    };

    
    /**
     * The above code defines a struct named "entry" with three integer fields: "node", "start_time",
     * and "end_time".
     * @property {int} node - An integer representing a node in a graph or network.
     * @property {int} start_time - The start time of an entry, which represents the time at which a
     * certain event or task begins. It is an integer value.
     * @property {int} end_time - The end_time property represents the time at which a certain event or
     * task ends. It is an integer value that indicates the specific time, usually measured in units
     * such as seconds or minutes.
     */
    struct entry {
        int node;
        int start_time;
        int end_time;
    };

    /**
     * The struct "wedge" represents a time interval with a start time and an end time.
     * @property {int} start_time - An integer representing the start time of the wedge.
     * @property {int} end_time - The end_time property represents the end time of a wedge.
     */
    struct wedge {
        int start_time;
        int end_time;
    };

    /**
     * The struct "wedge1" represents a wedge with a start time, end time, and a value.
     * @property {int} start_time - The start time of the wedge1 struct. It represents the time at
     * which the wedge1 starts.
     * @property {int} end_time - The end_time property represents the end time of a wedge.
     * @property {int} v - The property "v" represents the value associated with the wedge.
     */
    struct wedge1 {
        int start_time;
        int end_time;
        int v;
    };

    /**
     * The struct "vertices_pair" represents a pair of integers, typically used to represent a source
     * and destination vertex in a graph.
     * @property {int} src - An integer representing the source vertex in a graph or a pair of vertices
     * in a graph.
     * @property {int} dst - The "dst" property in the struct represents the destination vertex in a
     * graph. It is used to store the index or identifier of the vertex that is connected to the source
     * vertex.
     */
    struct vertices_pair {
        int src;
        int dst;
    };

    /**
     * The `pairCompare` struct is a comparison operator used to compare `algo::vertices_pair` objects
     * based on their `src` and `dst` members.
     * 
     */
    struct pairCompare {
        bool operator() (algo::vertices_pair const &u1, algo::vertices_pair const &u2) const
        {
            return (u1.src < u2.src || u1.src == u2.src && u1.dst < u2.dst);
        }
    };

    /**
     * The above code defines a struct named "third_entry" with five integer fields.
     * @property {int} src - An integer representing the source vertex of an edge in a graph.
     * @property {int} start_time - The start time of an event or process. It represents the time at
     * which the event or process begins.
     * @property {int} vertex_id - The vertex_id property is an integer that represents the unique
     * identifier of a vertex in a graph.
     */
    struct third_entry {
        int src, dst;
        int start_time, end_time;
        int vertex_id;
    };

    #pragma pack(1)
    struct list_removal {
        int reach;
        robin_hood::unordered_flat_map<int, int> removal;
    };

    #pragma pack(1)
    struct count_removal {
        int reach = 0;
        robin_hood::unordered_flat_map<int, short> removal;
    };

    /**
     * The TBP_index struct represents an index with two vectors, L_out and L_in, each containing
     * vectors of algo::entry objects.
     * @property L_out - L_out is a vector of vectors of type algo::entry. It represents the outgoing
     * edges of a graph. Each element in the outer vector represents a vertex in the graph, and each
     * element in the inner vector represents an outgoing edge from that vertex.
     * @property L_in - L_in is a vector of vectors of algo::entry objects. It represents the incoming
     * edges of a graph. Each element in the outer vector represents a vertex in the graph, and each
     * element in the inner vector represents an incoming edge to that vertex.
     */
    struct TBP_index {
        TBP_index(size_t a, size_t b) : L_out(a), L_in(b) {};
        std::vector<std::vector<algo::entry>> L_out;
        std::vector<std::vector<algo::entry>> L_in;
    };

    class Core {
    public:
        std::vector<edge> edges;
        algo::Graph graph;
        
        Core(void);

        int main(int ac, char **av);
        void init_data(std::string const &src);
        void init_data3(std::string const &src);

        static std::vector<int> get_vertex_order(algo::Graph const &graph);
        static std::vector<algo::partition> get_partitions_of_node(std::vector<algo::edge> const &edges, algo::Node &node);

        algo::TBP_index TBP_build(std::vector<int> const &vertex_order, algo::Graph &graph);
        bool TBP_query(TBP_index const &index, Node const &u, Node const &w, int s, int e);
        bool edge_colliding(int a, int b);
        std::vector<std::vector<algo::entry>> inverted_in_label_set(std::vector<std::vector<algo::entry>> L_in);
        int SSRQ(algo::Node &node, std::vector<std::vector<algo::entry>> &L_out, std::vector<std::vector<algo::entry>> &invL_in);
        int maxR_baseline(std::vector<std::vector<algo::entry>> &L_out, std::vector<std::vector<algo::entry>> &invL_in);
        std::tuple<int, std::vector<int>> maxReach(algo::Graph &graph);
        algo::Graph remove_vertex(Node &node, Graph &graph);
        int choose_vertex(algo::Graph &graph);
        std::vector<int> greedy1(algo::Graph &graph, int k);

        // third index
        void sort_te(std::vector<algo::wedge> &list_w, algo::wedge &wedge);
        bool dominant_wedge(algo::wedge &wedge1, algo::wedge &wedge2);
        std::vector<algo::wedge> necessary_wedges(std::vector<algo::wedge> &list_w, algo::wedge &wedge);
        std::vector<std::map<algo::vertices_pair, std::vector<algo::wedge>, algo::pairCompare>> third_index_build(algo::Graph &graph, int const &start, int const &end);
        std::vector<algo::third_entry> sort_list_entries(std::vector<std::map<algo::vertices_pair, std::vector<algo::wedge>, algo::pairCompare>> &list_dict);
        int third_choose_vertex_v2(std::vector<std::map<algo::vertices_pair, std::vector<algo::wedge>, algo::pairCompare>> &list_dict, int &size);
        std::vector<int> greedy3_v2(algo::Graph &graph, int k, int const &start, int const &end);
        std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>> list_position(std::vector<algo::third_entry> &LEO);
        bool SPRQ_static2(std::vector<algo::third_entry> &LEO, std::vector<std::vector<int>> &src_pos, std::vector<std::vector<int>> &dst_pos, int &idu, int &idw);
        std::vector<int> SSRQ_static2(std::vector<algo::third_entry> &LEO, std::vector<std::vector<int>> &src_pos, int &idu);
        int max_reach_static(std::vector<algo::third_entry> &LEO);
        int max_reach_static2(std::vector<algo::third_entry> &LEO, std::vector<std::vector<int>> &src_pos);
        std::vector<algo::third_entry> vertex_removal(std::vector<algo::third_entry> &LEO, int vertex_id, bool layer);
        std::vector<algo::third_entry> edge_removal(std::vector<algo::third_entry> &LEO, algo::edge &edge);
        bool edge_colliding2(algo::edge const &e1, algo::edge const &e2);
        void entry_insertion(std::vector<algo::third_entry> &LEO, algo::third_entry &entry);
        std::vector<algo::third_entry> edge_addition(std::vector<algo::third_entry> &LEO, algo::edge &edge);

        int test_SPRQ_baseline(TBP_index const &index, std::vector<int> test, std::vector<int> test2);
        int test_SPRQ_WTB(std::vector<algo::third_entry> &LEO, std::vector<std::vector<int>> &src_pos, std::vector<std::vector<int>> &dst_pos, std::vector<int> test, std::vector<int> test2);
        int test_SSRQ_baseline(std::vector<std::vector<algo::entry>> &L_out, std::vector<std::vector<algo::entry>> &invL_in, std::vector<int> test);
        int test_SSRQ_WTB(std::vector<algo::third_entry> &LEO, std::vector<std::vector<int>> &src_pos, std::vector<int> test);
        int test_maxR_baseline(std::vector<std::vector<algo::entry>> &L_out, std::vector<std::vector<algo::entry>> &invL_in);
        int test_maxR_WTB(std::vector<algo::third_entry> &LEO, std::vector<std::vector<int>> &src_pos);
        int test_VRCQ_baseline(algo::Graph &graph);
        int test_VRCQ_WTB(std::vector<std::map<algo::vertices_pair, std::vector<algo::wedge>, algo::pairCompare>> &list_dict, int size);
    };
    
    int rndm(int a, int b, double g, int size);
    algo::partition const *find_partition(std::vector<algo::partition> const &partitions, int edge_id);
}

#endif