#include "debruijngraph.h"

DeBruijnGraph::DeBruijnGraph(string str, int _k) {
    k = _k;
    create_graph_by_string(str);
}

DeBruijnGraph::DeBruijnGraph(vector<string> strs) {
    create_graph_by_vector(strs);
}
DeBruijnGraph::DeBruijnGraph(vector<string> strs, int _k, int _d) : k(_k), d(_d) {
    create_graph_by_paired_vector(strs);
}

// bool DeBruijnGraph::equal_k_prefix_and_suffix(const string& a, const string& b){

//    string s1 = a.substr(1, a.size() - 1) ;
//    string s2 = b.substr(0, a.size() - 1);
//    return  s1 == s2;
//}


void DeBruijnGraph::make_counts() {
    for (const auto n : g) {
        counts[n.first] += n.second.size();
        for (size_t i = 0; i < n.second.size(); ++i) {
            --counts[n.second[i]];
        }
    }
}
void DeBruijnGraph::make_counts_paired() {
    for (const auto n : pg) {
        pcounts[n.first] += n.second.size();
        for (size_t i = 0; i < n.second.size(); ++i) {
            --pcounts[n.second[i]];
        }
    }
}

Path DeBruijnGraph::find_eulerian_cycle() {
    make_counts();
    Path p;
    // string start =find_unbalanced_nodes_and_add_edge();
    //    find_circuit(start,p); //Must be reversed
    //    p.erase(p.begin()); //Delete for the path
    //    return p;
}

PairedPath DeBruijnGraph::find_eulerian_cycle_paired() {
    make_counts_paired();
    PairedPath p;
    Pair start = find_unbalanced_nodes_and_add_edge();
    find_circuit_paired(start, p); // Must be reversed
    //  p.erase(p.begin()); //Delete for the path
    std::reverse(p.begin(), p.end());
    return p;
}

void DeBruijnGraph::find_circuit(string from, Path& circuit) {
    if (g[from].size() == 0) {
        circuit.push_back(from);
    } else {
        while (g[from].size() > 0) {
            string to = g[from].back();
            g[from].pop_back();
            find_circuit(to, circuit);
        }
        circuit.push_back(from);
    }
}

void DeBruijnGraph::find_circuit_paired(Pair from, PairedPath& circuit) {
    if (pg[from].size() == 0) {
        circuit.push_back(from);
    } else {
        while (pg[from].size() > 0) {
            std::random_shuffle(pg[from].begin(), pg[from].end());
            Pair to = pg[from].back();
            pg[from].pop_back();
            find_circuit_paired(to, circuit);
        }
        circuit.push_back(from);
    }
}

Pair DeBruijnGraph::find_unbalanced_nodes_and_add_edge() {
    Pair first, second;
    for (auto pair : pcounts) {
        if (pair.second != 0) {
            if (pair.second > 0)
                first = pair.first;
            else
                second = pair.first;
        }
    }
    pg[second].push_back(first);
    return first;
}

void DeBruijnGraph::create_graph_by_string(string s) {
    for (size_t i = 0; i < s.size() - k + 1; ++i) {
        string a = s.substr(i, k);
        string vertex = s.substr(i, k - 1);

        for (size_t j = i; j < s.size() - k + 1; ++j) {
            string b = s.substr(j, k - 1);
            string neig = s.substr(j + 1, k - 1);
            if (vertex == b)
                g[vertex].push_back(neig);
        }
    }
}

void DeBruijnGraph::create_graph_by_vector(vector<string> s) {
    int size = s[0].size() - 1;
    k = size;
    string first, second;
    for (string str : s) {
        first = str.substr(0, size);
        second = str.substr(1, size);
        g[first].push_back(second);
    }
}

void DeBruijnGraph::create_graph_by_paired_vector(vector<string> s) {
    // int size = s[0].size() - 1;
    for (string str : s) {
        Pair p = make_split_pair(str);
        string prefix_p_f = p.front.substr(0, k - 1);
        string prefix_p_b = p.back.substr(0, k - 1);
        string suff_p_f = p.front.substr(1, k - 1);
        string suff_p_b = p.back.substr(1, k - 1);

        Pair prefix_pair(prefix_p_f, prefix_p_b);
        Pair suffix_pair(suff_p_f, suff_p_b);
        pg[prefix_pair].push_back(suffix_pair);
    }
}

void DeBruijnGraph::read_graph_from_file(ifstream& f) {
    string line;
    string first, rest;
    while (getline(f, line, '\n')) {
        stringstream linestream(line);
        getline(linestream, first, ' ');

        linestream.ignore(3);

        while (getline(linestream, rest, ',')) {
            g[first].push_back(rest);
        }
    }
}

void DeBruijnGraph::print_eulerian_cycle(ofstream& f) {
    Path p = find_eulerian_cycle();
    size_t i;
    for (i = p.size() - 1; i > 0; --i) {
        f << p[i] << "->";
    }
    f << p[i];
    f << endl;
}

void DeBruijnGraph::print_eulerian_cycle_as_DNA(ofstream& f) {
    Path p = find_eulerian_cycle();
    string DNA = p[p.size() - 1];
    size_t k = DNA.size() - 1;

    for (int i = p.size() - 2; i >= 0; --i) {
        DNA += p[i][k];
    }
    f << DNA << endl;
}


string DeBruijnGraph::string_spelled_by_genome_path_problem(Path data) {
    string DNA = data[0];
    size_t k = DNA.size() - 1;

    for (size_t i = 1; i < data.size(); ++i) {
        DNA += data[i][k];
    }

    return DNA;
}

void DeBruijnGraph::string_spelled_by_gapped_patterns(ofstream& f) {
    PairedPath path = find_eulerian_cycle_paired();
    vector<string> firsts, seconds;
    for (auto a : path) {
        firsts.push_back(a.front);
        seconds.push_back(a.back);
    }

    string prefix_string = string_spelled_by_genome_path_problem(firsts);
    string suffix_string = string_spelled_by_genome_path_problem(seconds);
    f << prefix_string << endl;
    f << "    " << suffix_string << endl;
    for (int i = k + d + 1; i < prefix_string.size(); ++i)
        if (prefix_string[i] != suffix_string[i - k - d])
            cout << "there is no string spelled by the gapped patterns" << endl;

    string result = prefix_string + suffix_string.substr(suffix_string.size() - k - d, k + d);
    f << result << endl;
}

void DeBruijnGraph::print_max_nonbranching_paths(ofstream& f) {
    auto paths = maximal_nonBranching_paths();
    for (auto path : paths) {
        f << string_spelled_by_genome_path_problem(path);
        f << endl;
    }
}

vector<Path> DeBruijnGraph::maximal_nonBranching_paths() {
    make_counts();
    vector<Path> paths;
    Path nonbranching_path;
    for (auto v : g) {
        if (g[v.first].size() != 1 || counts[v.first] != 0)
            if (v.second.size() > 0) {
                for (auto edge : v.second) {
                    nonbranching_path = {v.first};
                    visited[v.first] = true;
                    while (g[edge].size() == 1 && counts[edge] == 0) {

                        nonbranching_path.push_back(edge);
                        visited[edge] = true;
                        edge = g[edge].back();
                    }
                    nonbranching_path.push_back(edge);
                    paths.push_back(nonbranching_path);
                }
            }
    }

    for (auto v : g) {
        if (counts[v.first] == 0 && visited[v.first] == false) {
            nonbranching_path = {v.first};
            auto edge = v.second.back();
            while (edge != v.first) {
                nonbranching_path.push_back(edge);
                visited[edge] = true;
                edge = g[edge].back();
            }
            nonbranching_path.push_back(edge);
            paths.push_back(nonbranching_path);
        }
    }

    return paths;
}


Pair DeBruijnGraph::make_split_pair(string s) {
    string first = s.substr(0, s.find('|'));
    s.erase(0, s.find('|') + 1);
    string second = s;
    return Pair(first, second);
}


ostream& operator<<(ostream& os, const DeBruijnGraph& graph) {
    for (const auto n : graph.g) {
        if (n.second.size() != 0) {
            os << n.first << " -> ";
            size_t i = 0;
            if (n.second.size() > 1)
                for (i = 0; i < n.second.size() - 1; ++i) {
                    os << n.second[i] << ",";
                }
            os << n.second[i] << std::endl;
        }
    }
    return os;
}
