#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <stdio.h>
#include <algorithm>
#include "debruijngraph.h"
#include <stdio.h>
#include <cmath>

using namespace std;


// FILE* a  = fopen("answer.txt", "w");
ofstream a("answer.txt");
vector<string> composition_k(string DNA, int k) {
    set<string> s;
    for (size_t i = 0; i < DNA.length() - k + 1; ++i) {
        s.insert(DNA.substr(i, k));
    }
    vector<string> result(s.begin(), s.end());
    return result;
}

string string_spelled_by_genome_path_problem(vector<string> data) {
    string DNA = data[0];
    size_t k = DNA.size() - 1;

    for (size_t i = 1; i < data.size(); ++i) {
        DNA += data[i][k];
    }

    return DNA;
}

bool equal_prefix_and_suffix(const string& a, const string& b) {
    return a.substr(1, a.size() - 1) == b.substr(0, a.size() - 1);
}

void overlap_graph_problem(vector<string>& data) {

    sort(data.begin(), data.end());

    for (auto seq : data) {
        for (auto seq2 : data)
            if (equal_prefix_and_suffix(seq, seq2)) {
                ;
            }
        // fprintf(a, "%s -> %s\n",seq.c_str(), seq2.c_str());
    }
}

void binary(int n, string& A, vector<string>& data) {
    if (n < 1)
        data.push_back(A);
    else {
        A[n - 1] = '0';
        binary(n - 1, A, data);
        A[n - 1] = '1';
        /*
        Feel free to copy but please acknowledge studyalgorithms.com
        */
        binary(n - 1, A, data);
    }
}

vector<pair<string, string>> generate_k_d_mer(string str, int k, int d) {
    vector<pair<string, string>> s;
    for (size_t i = 0; i < str.length() - 2 * k - d + 1; ++i) {
        string first = str.substr(i, k);
        string second = str.substr(i + k + d, k);
        s.push_back(make_pair(first, second));
    }
    sort(s.begin(), s.end());
    vector<pair<string, string>> result(s.begin(), s.end());
    return result;
}
pair<string, string> make_split_pair(string s) {
    string first = s.substr(0, s.find('|'));
    s.erase(0, s.find('|') + 1);
    string second = s;
    return make_pair(first, second);
}

string string_spelled_by_gapped_patterns(vector<string> gappedPatterns, int d) {

    vector<string> firsts, seconds;
    for (auto a : gappedPatterns) {
        auto p = make_split_pair(a);
        firsts.push_back(p.first);
        seconds.push_back(p.second);
    }

    int k = firsts[0].size();
    string prefix_string = string_spelled_by_genome_path_problem(firsts);
    string suffix_string = string_spelled_by_genome_path_problem(seconds);

    for (size_t i = k + d + 1; i < prefix_string.size(); ++i)
        if (prefix_string[i] != suffix_string[i - k - d])
            cout << "there is no string spelled by the gapped patterns" << endl;

    string result = prefix_string + suffix_string.substr(suffix_string.size() - k - d, k + d);
    return result;
}


int main() {
    ifstream f("data.txt");
    string str;
    vector<string> data;

    //    auto result =  generate_k_d_mer("TAATGCCATGGGATGTT",3,2);
    //    for(auto pair:result){
    //        cout<<"("<<pair.first<<"|"<<pair.second<<")"<<" ";

    //    }
    //    cout<<endl;
    //    int n = 8; //HAVE TO DELETE N - 2 from the final answer !!!
    //    string A(n,'0');
    //    binary(n,A,data);
    // int k, d;
    //  f >>k>> d;
    while (f >> str)
        data.push_back(str);

    // a<<string_spelled_by_gapped_patterns(data,d);

    DeBruijnGraph g(data);
    //  g.read_graph_from_file(f);
    g.print_max_nonbranching_paths(a);
    // g.string_spelled_by_gapped_patterns(a);
    //  cout << g;
    return 0;
}
