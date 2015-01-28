#ifndef DEBRUIJNGRAPH_H
#define DEBRUIJNGRAPH_H

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
using namespace std;

typedef vector<string> Neighbours;

struct Pair {
    pair<string, string> p;
    string front, back;
    Pair() {
    }
    Pair(string f, string b) : front(f), back(b) {
        p = make_pair(f, b);
    }

    bool operator<(const Pair& other) const {
        return (p < other.p);
    }
    bool operator==(const Pair& other) const {
        return front == other.front && back == other.back;
    }
};
typedef vector<Pair> PairNeighbours;

struct Edge {
    string from, to;
    Edge(string f, string t) : from(f), to(t) {
    }
    bool operator<(const Edge& other) const {
        return from <= other.from;
    }
    bool operator==(const Edge& other) const {
        return from == other.from && to == other.to;
    }
};

typedef map<string, Neighbours> AdjList;
typedef map<Pair, PairNeighbours> PairAdjList;
typedef map<string, bool> Visited;
typedef vector<string> Path;
typedef vector<Pair> PairedPath;
typedef vector<Edge> Edges;
typedef map<string, int> Counts;
typedef map<Pair, int> PCounts;

class DeBruijnGraph {
    public:
    DeBruijnGraph() {
    }
    DeBruijnGraph(string str, int k);
    DeBruijnGraph(vector<string> strs);
    DeBruijnGraph(vector<string> strs, int k, int d);

    friend ostream& operator<<(ostream& os, const DeBruijnGraph& graph);
    void read_graph_from_file(ifstream& f);
    void print_eulerian_cycle(ofstream& f);
    void print_eulerian_cycle_as_DNA(ofstream& f);
    void string_spelled_by_gapped_patterns(ofstream& f);
    void print_max_nonbranching_paths(ofstream& f);

    private:
    AdjList g;
    PairAdjList pg;
    Visited visited = {};
    Counts counts = {};
    PCounts pcounts = {};

    Edges edges;

    int k, d;

    void create_graph_by_string(string s);
    void create_graph_by_vector(vector<string> s);
    void create_graph_by_paired_vector(vector<string> s);

    void make_counts();
    void make_counts_paired();
    Path find_eulerian_cycle();
    PairedPath find_eulerian_cycle_paired();
    void find_circuit(string from, Path& circuit);
    void find_circuit_paired(Pair from, PairedPath& circuit);


    vector<Path> maximal_nonBranching_paths();

    Pair make_split_pair(string s);
    Pair find_unbalanced_nodes_and_add_edge();
    string string_spelled_by_genome_path_problem(Path data);
};

#endif // DEBRUIJNGRAPH_H
