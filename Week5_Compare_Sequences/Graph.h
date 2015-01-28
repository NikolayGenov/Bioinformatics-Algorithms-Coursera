#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <list>
using namespace std;

struct Edge {
    string from, to;
    int weight;
    Edge():from("-1"),to("-1"),weight(-1){}
    Edge(string f, string t,int c) : from(f), to(t), weight(c) {
    }
    bool operator<(const Edge& other) const {
        return from <= other.from;
    }
    bool operator==(const Edge& other) const {
        return from == other.from && to == other.to;
    }
};

typedef string Vert;
typedef vector<string> Neighbours;
typedef vector<Edge> Edges;
typedef map<string, Neighbours> AdjList;
typedef map<string, Edges> AdjListEdges;
typedef map<string, bool> Visited;
typedef vector<string> Verts;
typedef map<string, int> Counts;
typedef vector<string> Path;


class Graph {
    public:
    Graph() {
    }
    Graph(string str, int k);
    Graph(vector<string> strs);

    friend ostream& operator<<(ostream& os, const Graph& graph);
    void read_graph_from_file(ifstream& f);
    void read_graph_from_file_by_edges(ifstream& f);
    void print_topo_order(ofstream& f);
       void longest_path(ostream &f);

    private:
    AdjList g;
    Verts verts;
    Counts counts = {};
    Visited visited;
    AdjListEdges preds = {};
    AdjListEdges ge = {};
    Vert source,sink;
    map<string, bool> has_incoming_edges = {};
    int k, d;

    void create_graph_by_string(string s);
    void create_graph_by_vector(vector<string> s);
    void make_counts();
    void check_for_incoming_edges();
    void check_for_incoming_edges_with_weights();


    Path recreate_path(map<string, Edge> &paths);
    Path topological_ordering();
    Path topological_ordering_by_edges();

};

#endif // GRAPH_H
