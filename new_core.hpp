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
     * We define a structure called "edge" with four integer fields and a print function.
     * @property {int} src - The source vertex of the edge.
     * @property {int} dst - The destination vertex of the edge.
     * @property {int} start_time - The start time of the edge, indicating when the edge becomes active
     * @property {int} end_time - The ending time of the edge, indicating when the edge stops being active
     * or ends.
     */
    struct edge {
        int src, dst;
        int start_time, end_time;

        void print(void) {
            std::cout << src << " " << dst << " " << start_time << " " << end_time << std::endl;
        }
    };

    typedef std::vector<int> partition;

    /* The Node class represents a node in a graph, with an ID, order, edges, and partitions.
    For a given graph, nodes are ranked by decreasing degree (when degree is equal, ranked by increasing ID),
    the node of order 0 is the one of highest degree.
    */
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

    /* The Graph class represents a bipartite graph with two sets of nodes, upper and lower. */
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

    struct entry {
        int node;
        int start_time;
        int end_time;
    };

    struct wedge {
        int start_time;
        int end_time;
    };

    struct wedge1 {
        int start_time;
        int end_time;
        int v;
    };

    struct vertices_pair {
        int src;
        int dst;
    };

    struct pairCompare {
        bool operator() (algo::vertices_pair const &u1, algo::vertices_pair const &u2) const
        {
            return (u1.src < u2.src || u1.src == u2.src && u1.dst < u2.dst);
        }
    };

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
    
    int rndm(int a, int b, double g);
    algo::partition const *find_partition(std::vector<algo::partition> const &partitions, int edge_id);
}

#endif