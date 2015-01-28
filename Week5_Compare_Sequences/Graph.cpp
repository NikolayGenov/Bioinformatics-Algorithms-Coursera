#include "Graph.h"

Graph::Graph(string str, int _k) {
    k = _k;
    create_graph_by_string(str);
}

Graph::Graph(vector<string> strs) {
    create_graph_by_vector(strs);
}

Path Graph::topological_ordering() {
    Path list;
    check_for_incoming_edges();
    Verts candidates;
    for (auto n : has_incoming_edges)
        if (!n.second)
            candidates.push_back(n.first);

    while (!candidates.empty()) {
        Vert a = candidates.back();
        list.push_back(a);
        candidates.pop_back();
        for (int i = g[a].size() - 1; i >= 0; --i) {
            Vert b = g[a].back();
            g[a].pop_back();
            check_for_incoming_edges();
            if (!has_incoming_edges[b])
                candidates.push_back(b);
        }
    }
    bool is_empty = true;
    for (auto n : g)
        is_empty = is_empty && n.second.empty();
    if (!is_empty) {
        cout << "Error!!!" << endl;
        exit(1);
    }
    return list;
}
Path Graph::topological_ordering_by_edges() {
    Path list;
    check_for_incoming_edges_with_weights();
    Verts candidates;
    for (auto n : has_incoming_edges)
        if (!n.second)
            candidates.push_back(n.first);

    while (!candidates.empty()) {
        Vert a = candidates.back();
        list.push_back(a);
        candidates.pop_back();
        for (int i = ge[a].size() - 1; i >= 0; --i) {
            Edge e = ge[a].back();
            Vert b = e.to;
            ge[a].pop_back();
            check_for_incoming_edges_with_weights();
            if (!has_incoming_edges[b])
                candidates.push_back(b);
        }
    }
    bool is_empty = true;
    for (auto n : g)
        is_empty = is_empty && n.second.empty();
    if (!is_empty) {
        cout << "Error!!!" << endl;
        exit(1);
    }
    return list;
}

Path Graph::recreate_path(map<string, Edge>& paths) {
    Path path;
    Vert target = sink;
    while (target != source) {
        path.push_back(target);
        target = paths[target].to;
    }
    path.push_back(source);
    std::reverse(path.begin(), path.end());
    return path;
}

void Graph::longest_path(ostream& f) {
    map<string, int> s;
    map<string, Edge> paths;
    for (auto v : verts)
        s[v] = -10000000;
    s[source] = 0;
    Path topo_order = topological_ordering_by_edges();
    auto iter_source = find(topo_order.begin(), topo_order.end(), source);
    for (auto b = iter_source; b != topo_order.end(); ++b) {
        if (*b != source) {
            int max_val = -1000000;
            for (auto edge : preds[*b]) {
                int val = s[edge.to] + edge.weight;
                if (val >= max_val) {
                    max_val = val;
                    paths[*b] = edge;
                }
            }
            s[*b] = max_val;
        }
    }
    f << s[sink] << endl;

    Path path = recreate_path(paths);
    for (auto p : path)
        f << p << "->";
    f << endl;
}

void Graph::print_topo_order(ofstream& f) {
    auto list = topological_ordering();
    for (auto s : list)
        f << s << ", ";
    f << endl;
}

void Graph::check_for_incoming_edges() {
    has_incoming_edges = {};
    for (auto vert : verts) {
        bool has_vert = false;
        for (auto adjlist : g) {
            if (has_vert)
                break;
            for (auto v : adjlist.second)
                if (vert == v) {
                    has_vert = true;
                    break;
                }
        }
        has_incoming_edges[vert] = has_vert;
    }
}
void Graph::check_for_incoming_edges_with_weights() {
    has_incoming_edges = {};
    for (auto vert : verts) {
        bool has_vert = false;
        for (auto adjlist : ge) {
            if (has_vert)
                break;
            for (auto v : adjlist.second)
                if (vert == v.to) {
                    has_vert = true;
                    break;
                }
        }
        has_incoming_edges[vert] = has_vert;
    }
}


void Graph::make_counts() {
    counts = {};
    for (const auto n : g) {
        counts[n.first] += n.second.size();
        for (size_t i = 0; i < n.second.size(); ++i) {
            --counts[n.second[i]];
        }
    }
}

void Graph::create_graph_by_string(string s) {
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

void Graph::create_graph_by_vector(vector<string> s) {
    int size = s[0].size() - 1;
    k = size;
    string first, second;
    for (string str : s) {
        first = str.substr(0, size);
        second = str.substr(1, size);
        g[first].push_back(second);
    }
}


void Graph::read_graph_from_file(ifstream& f) {
    string line;
    string first, rest;
    while (getline(f, line, '\n')) {
        stringstream linestream(line);
        getline(linestream, first, ' ');

        linestream.ignore(3);

        while (getline(linestream, rest, ',')) {
            g[first].push_back(rest);
            if (g[rest].empty())
                g[rest] = {};
        }
    }
    verts = {};
    for (auto n : g)
        verts.push_back(n.first);
}

void Graph::read_graph_from_file_by_edges(ifstream& f) {
    f >> source >> sink;
    string line;
    string first, second;
    int w;

    getline(f, line, '\n');
    while (getline(f, line, '\n')) {
        stringstream linestream(line);
        getline(linestream, first, '-');
        linestream.ignore(1);
        getline(linestream, second, ':');
        linestream >> w;
        ge[first].emplace_back(Edge(first, second, w));
        preds[second].emplace_back(Edge(second, first, w));
        if (ge[second].empty())
            ge[second] = {};
        if (preds[first].empty())
            preds[first] = {};
    }

    verts = {};
    for (auto n : ge)
        verts.push_back(n.first);
}

ostream& operator<<(ostream& os, const Graph& graph) {
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
