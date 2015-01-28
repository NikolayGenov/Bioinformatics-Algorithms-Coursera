#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <map>
using namespace std;

typedef vector<int> Permutation;
typedef vector<int> Chromosome;
typedef vector<Chromosome> Chromosomes;


struct Edge {
    int from, to;
    Edge() : from(-1), to(-1){};
    Edge(int a, int b) : from(a), to(b) {
    }
    bool operator==(const Edge& other) const {
        return (from == other.from && to == other.to) || (to == other.from && from == other.to);
    }
};
typedef vector<Edge> Edges;
typedef vector<Edges> Cycles;

ofstream a("answer.txt");

Permutation read_input(string filename) {
    ifstream f(filename.c_str());
    Permutation data;
    f.ignore(1);
    int tmp;
    while (f.peek() != ')') {
        f >> tmp;
        data.push_back(tmp);
    }
    return data;
}
Chromosomes read_input_multiple(ifstream& f) {
    f.ignore(1);
    Chromosomes result;
    do {
        Permutation data;
        int tmp;
        while (f.peek() != ')') {
            f >> tmp;
            data.push_back(tmp);
        }
        f.ignore(1);
        result.push_back(data);
        if (f.peek() == '\n') {
            f.ignore(1);

            break;
        }
        f.ignore(1);

    } while (true);
    return result;
}
Edges read_edges(FILE* f) {
    Edges result;
    while (true) {
        int from, to;
        fscanf(f, "(%d, %d)", &from, &to);
        Edge data(from, to);
        result.push_back(data);
        if (fgetc(f) == '\n')
            break;
        if (fgetc(f) == '\n')
            break;
    }
    return result;
}


string permutation_to_string(const Permutation& perm) {
    string output = "(";
    for (auto s : perm) {
        if (s > 0)
            output += "+";
        output += to_string(s) + " ";
    }
    output[output.size() - 1] = ')';
    return output;
}
string chromosomes_to_string(const Chromosomes& P) {
    string output;
    for (auto ch : P)
        output += permutation_to_string(ch);
    return output;
}

string permutation_without_signs_to_string(const Permutation& perm) {
    string output = "(";
    for (auto s : perm)
        output += to_string(s) + " ";

    output[output.size() - 1] = ')';
    return output;
}

string edges_to_string(const Edges& edges) {
    string text;
    for (Edge e : edges) {
        text += "(" + to_string(e.from) + ", " + to_string(e.to) + "), ";
    }
    text[text.size() - 2] = ' ';
    return text;
}

void reverse_between(Permutation& perm, int i, int j) {
    vector<int>::const_iterator first = perm.begin() + i;
    vector<int>::const_iterator last = perm.begin() + j + 1;
    Permutation subvec(first, last);
    reverse(subvec.begin(), subvec.end());
    for (int k = 0; i < j + 1; ++i, ++k)
        perm[i] = -subvec[k];
}

int greedy_sorting(Permutation& perm) {
    auto approx_reversal_distance = 0;
    int k = 1;
    for (size_t i = 0; i < perm.size(); ++i, ++k) {
        if (k != abs(perm[i])) {
            auto it = find(perm.begin(), perm.end(), k);
            if (it == perm.end())
                it = find(perm.begin(), perm.end(), -k);
            int j = distance(perm.begin(), it);

            reverse_between(perm, i, j);
            ++approx_reversal_distance;
            a << permutation_to_string(perm) << endl;
        }

        if (perm[i] == -k) {
            perm[i] = k;
            ++approx_reversal_distance;
            a << permutation_to_string(perm) << endl;
        }
    }
    return approx_reversal_distance;
}
int count_breakpoints(Permutation perm) {
    perm.push_back(perm.size() + 1);
    perm.insert(perm.begin(), 0);
    int breakpoints = 0;
    for (size_t i = 0; i < perm.size() - 1; ++i) {
        if ((perm[i + 1] - perm[i]) != 1)
            ++breakpoints;
    }
    return breakpoints;
}

Permutation chromosome_to_cycle(Permutation chromosome) {
    Permutation node;
    for (size_t j = 0; j < chromosome.size(); ++j) {
        int i = chromosome[j];
        if (i > 0) {
            node.push_back(2 * i - 1);
            node.push_back(2 * i);
        } else {
            node.push_back(-2 * i);
            node.push_back(-2 * i - 1);
        }
    }
    return node;
}

Permutation cycle_to_chromosome(Permutation nodes) {
    Permutation chromosome;
    for (size_t j = 0; j < nodes.size() / 2; ++j)
        if (nodes[2 * j] < nodes[2 * j + 1])
            chromosome.push_back(nodes[2 * j + 1] / 2);
        else
            chromosome.push_back(-nodes[2 * j] / 2);

    return chromosome;
}

Edges colored_edges(Chromosomes p) {
    Edges edges;

    for (Chromosome chromosome : p) {
        Permutation perm = chromosome_to_cycle(chromosome);
        int sz = perm.size();
        for (size_t j = 0; j < chromosome.size(); ++j)
            edges.push_back({perm[2 * j + 1 % sz], perm[2 * (j + 1) % sz]});
    }
    return edges;
}
Edges black_edges(Chromosomes p) {
    Edges edges;

    for (Chromosome chromosome : p) {
        Permutation perm = chromosome_to_cycle(chromosome);
        int sz = perm.size();
        for (size_t j = 0; j < chromosome.size(); ++j)
            edges.push_back({perm[2 * j + 1 % sz], perm[2 * (j + 1) % sz]});
    }
    return edges;
}

Edges all_edges(Chromosomes p) {
    Edges edges;

    for (Chromosome chromosome : p) {
        Permutation perm = chromosome_to_cycle(chromosome);
        size_t sz = perm.size();
        for (size_t j = 0; j < sz; ++j)
            edges.push_back({perm[j], perm[(j + 1) % sz]});
    }
    return edges;
}
Chromosomes graph_to_genome(Edges genome_graph) {
    Chromosomes P;
    Chromosome chrom;
    for (Edge e : genome_graph) {
        if (e.from < e.to) {
            chrom.push_back(e.from);
            chrom.push_back(e.to);
        }

        else {
            chrom.insert(chrom.begin(), e.to);
            chrom.push_back(e.from);

            Chromosome ch = cycle_to_chromosome(chrom);
            P.push_back(ch);
            chrom.clear();
        }
    }
    return P;
}
void add_unique(Chromosome& ch, int val) {
    auto it = find(ch.begin(), ch.end(), val);
    if (it == ch.end())
        ch.push_back(val);
}

Chromosomes graph_to_genome_given_cycles(Cycles& cycles) {
    Chromosomes P;
    Chromosome chrom;
    for (Edges egs : cycles) {
        for (Edge e : egs) {
            add_unique(chrom, e.from);
            add_unique(chrom, e.to);
        }

        Chromosome ch = cycle_to_chromosome(chrom);
        P.push_back(ch);
        chrom.clear();
    }

    return P;
}
Edges::iterator find_edge(Edges& lst, int value) {
    for (Edges::iterator it = lst.begin(); it != lst.end(); ++it)
        if (it->from == value || it->to == value)
            return it;
    return lst.end();
}
bool follow_cycle(Edges& search_lst, Edge& e) {
    auto it = find_edge(search_lst, e.to);
    if (it != search_lst.end()) {
        e = *it;
        search_lst.erase(it, it + 1);
        return true;
    } else {
        auto it = find_edge(search_lst, e.from);
        if (it != search_lst.end()) {
            e = *it;
            search_lst.erase(it, it + 1);
            return true;
        }
    }
    return false;
}
Cycles find_all_cycles(Edges& list) {
    Cycles cls;
    Edge e;
    while (!list.empty()) {
        e = list.back();
        Edges current;
        while (follow_cycle(list, e))
            current.push_back(e);

        cls.push_back(current);
    }
    return cls;
}

int count_cycles(Edges& first, Edges& second) {
    int number_of_cycles = 0;
    bool is_first = true, res = false;
    Edge e;
    while (first.size() != 0) {
        if (!res) {
            e = first.back();
            first.pop_back();
            number_of_cycles++;
        }
        if (is_first) {
            res = follow_cycle(second, e);
            is_first = false;
        } else {
            res = follow_cycle(first, e);
            is_first = true;
        }
    }
    return number_of_cycles;
}
void remove_edge(Edges& lst, Edge e) {
    auto it = find(lst.begin(), lst.end(), e);
    if (it != lst.end())
        lst.erase(it);
}

void two_break_on_genome_graph(Edges& edges, int i, int i1, int j, int j1) {
    remove_edge(edges, Edge(i, j));
    edges.push_back({i1, j});
    remove_edge(edges, Edge(i1, j1));
    edges.push_back({i, j1});
}

Chromosomes break_on_genome(Chromosomes chrom, int i, int i1, int j, int j1) {
    Edges edges = all_edges(chrom);
    two_break_on_genome_graph(edges, i, i1, j, j1);
    auto cls = find_all_cycles(edges);
    auto P = graph_to_genome_given_cycles(cls);
    return P;
}

Chromosomes break_on_genome_by_edges(Edges edges, int i, int i1, int j, int j1) {
    two_break_on_genome_graph(edges, i, i1, j, j1);
    auto cls = find_all_cycles(edges);
    auto P = graph_to_genome_given_cycles(cls);
    return P;
}

bool has_nontrivial_cycle(Edges edges) {
    auto cycles = find_all_cycles(edges);
    for (auto cycle : cycles)
        if (cycle.size() != 1)
            return false;
    return true;
}

Edge pick_nontrivial_edge(Edges& blues, Edges breakpoint_graph) {
    auto cycles = find_all_cycles(breakpoint_graph);
    for (auto cycle : cycles)
        if (cycle.size() > 1)
            for (auto it = cycle.begin(); it != cycle.end(); ++it)
                if (find(blues.begin(), blues.end(), *it) != blues.end())
                    return *it;
    return Edge();
}

Edge find_edge_ending_at(Edges& lst, int j) {
    for (auto it = lst.begin(); it != lst.end(); ++it)
        if (it->from == j || it->to == j)
            return *it;
    return Edge();
}

Edges union_edges(const Edges& lst, const Edges& lst2) {
    Edges result = lst;
    for (Edge e : lst2)
        if (find(result.begin(), result.end(), e) == result.end())
            result.push_back(e);
    return result;
}

void shortest_rearrangement_scenario(Chromosomes& P, Chromosomes& Q) {
    a << chromosomes_to_string(P) << endl;

    auto red_edges = colored_edges(P);
    auto blue_edges = colored_edges(Q);
    auto breakpoint_graph = union_edges(red_edges, blue_edges);

    while (!has_nontrivial_cycle(breakpoint_graph)) {

        Edge j_i1 = pick_nontrivial_edge(blue_edges, breakpoint_graph);

        int j = j_i1.from, i1 = j_i1.to, i, j1;

        Edge i_j = find_edge_ending_at(red_edges, j);
        Edge i1_j1 = find_edge_ending_at(red_edges, i1);

        i = j == i_j.to ? i_j.from : i_j.to;
        j1 = i1 == i1_j1.from ? i1_j1.to : i1_j1.from;

        remove_edge(red_edges, i_j);
        remove_edge(red_edges, i1_j1);
        red_edges.push_back({i1, j});
        red_edges.push_back({i, j1});

        P = break_on_genome(P, i, i1, j, j1);
        breakpoint_graph = union_edges(red_edges, blue_edges);
        a << chromosomes_to_string(P) << endl;
    }
}
string reverse_complement(const string& input) {
    string rev(input.size(), ' ');
    for (int i = input.size() - 1, j = 0; i >= 0; --i, ++j)
        switch (input[i]) {
        case 'A':
            rev[j] = 'T';
            break;
        case 'T':
            rev[j] = 'A';
            break;
        case 'C':
            rev[j] = 'G';
            break;
        case 'G':
            rev[j] = 'C';
            break;
        default:
            break;
        }
    return rev;
}

long long count_k_mers(string& f, string& s, int k) {
    map<string, int> m;
    long long count = 0;
    long long sf = f.size() - k + 1;
    long long ss = s.size() - k + 1;
    for (long long i = 0; i < sf; ++i)
        m[f.substr(i, k)]++;
    for (long long i = 0; i < ss; ++i) {
        const string sub = s.substr(i, k);
        count += m[sub];
        count += m[reverse_complement(sub)];
    }
    return count;
}

int main() {
    FILE* fl = fopen("data.txt", "r");
    ifstream f("data.txt");

    int k;
    string s1, s2;
    f >> k >> s1 >> s2;
    // auto i = count_k_mers(s1, s2, k); //slow for the huge file
    //  cout << i << endl;
    //    auto P = read_input_multiple(f);
    //    auto Q = read_input_multiple(f);
    //
    //    int i, j, i1, j1;
    //    fscanf(fl, "%d, %d, %d, %d", &i, &i1, &j, &j1);

    // auto perms2 = read_input_multiple(f);
    //    auto first = colored_edges(perms);
    //  auto second = colored_edges(perms2);
    //    int blocks = first.size();
    //    int cycles = count_cycles(first,second);

    //    shortest_rearrangement_scenario(P, Q);

    fclose(fl);
    return 0;
}
