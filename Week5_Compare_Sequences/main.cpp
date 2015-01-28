#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "Graph.h"
#include <algorithm>
#include "BLOSUM62.h"
#include "PAM250.h"
#include <map>
using namespace std;

typedef vector<vector<int>> Matrix;

typedef vector<vector<string>> CharMatrix;
ofstream a("answer.txt");


void read_input(int& n, int& m, Matrix& down, Matrix& right, ifstream& f) {
    f >> n >> m;
    int number;
    char sep;
    down.resize(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m + 1; ++j) {
            f >> number;
            down[i].push_back(number);
        }
    f >> sep;
    right.resize(n + 1);
    for (int i = 0; i < n + 1; ++i)
        for (int j = 0; j < m; ++j) {
            f >> number;
            right[i].push_back(number);
        }
}
int manhattan_tourist_problem(int n, int m, const Matrix& down, const Matrix& right) {
    Matrix dp;
    size_t n_size = n + 1, m_size = m + 1;
    dp.resize(n_size);

    for (size_t i = 0; i < n_size; ++i)
        dp[i].resize(m_size);
    dp[0][0] = 0;
    for (size_t i = 1; i < n_size; ++i)
        dp[i][0] = dp[i - 1][0] + down[i - 1][0];

    for (size_t i = 1; i < m_size; ++i)
        dp[0][i] = dp[0][i - 1] + right[0][i - 1];

    for (size_t i = 1; i < n_size; ++i)
        for (size_t j = 1; j < m_size; ++j)
            dp[i][j] = max(dp[i - 1][j] + down[i - 1][j - 1], dp[i][j - 1] + right[i - 1][j - 1]);

    return dp[n][m];
}

CharMatrix LCS_backtrack(string& v, string& w, int& max_value, int indel_value) {
    CharMatrix backtrack;
    Matrix s;
    string sym;
    size_t v_size = v.size() + 1, w_size = w.size() + 1;
    int indel = indel_value; // TODO FIXME
    backtrack.resize(v_size);
    s.resize(v_size);
    for (size_t i = 0; i < v_size; ++i) {
        backtrack[i].resize(w_size);
        s[i].resize(w_size);
    }

    for (size_t i = 1; i < v_size; ++i)
        s[i][0] = -indel * i;

    for (size_t i = 1; i < w_size; ++i)
        s[0][i] = -indel * i;

    for (size_t i = 1; i < v_size; ++i)
        for (size_t j = 1; j < w_size; ++j) {
            char c1 = v[i - 1], c2 = w[j - 1];
            int match = s[i - 1][j - 1];        //+ BLOSUM62[c1][c2];
            int mismatch = s[i - 1][j - 1] - 1; // + BLOSUM62[c1][c2];
            int down = s[i - 1][j] - indel;
            int right = s[i][j - 1] - indel;

            if (c1 == c2)
                s[i][j] = match;
            else
                s[i][j] = max(max(down, right), mismatch);

            if (s[i][j] == match)
                sym = "match";
            if (s[i][j] == down)
                sym = "down";
            if (s[i][j] == right)
                sym = "right";
            if (s[i][j] == mismatch)
                sym = "mismatch";

            backtrack[i][j] = sym;
        }
    max_value = s[v.size()][w.size()];
    return backtrack;
}

// Not working right
/*
void output_LCS(CharMatrix& backtrack, string& v, string& w, int i, int j) {
    if (i < 0 || j < 0) {
        a << endl;
        return;
    }
    if (backtrack[i][j] == "↓") {
        output_LCS(backtrack, v, w, i - 1, j);
        a << "-";
    } else if (backtrack[i][j] == "→") {
        output_LCS(backtrack, v, w, i, j - 1);
        a << "-";
    } else if (backtrack[i][j] == "X") {
        output_LCS(backtrack, v, w, i - 1, j - 1);
        a << w[j];
    } else if (backtrack[i][j] == "↘") {
        output_LCS(backtrack, v, w, i - 1, j - 1);
        a << v[i];
    }
}
*/
string insert_indel(string str, int pos) {
    str.insert(str.begin() + pos, '-');
    return str;
}

void output_global_alignment(CharMatrix& backtrack, string& v, string& w) {
    string v_aligned = v, w_aligned = w;

    int i = v.size(), j = w.size();
    while (i >= 0 && j >= 0) {
        if (backtrack[i][j] == "down") {
            --i;
            w_aligned = insert_indel(w_aligned, j);
        } else if (backtrack[i][j] == "right") {
            --j;
            v_aligned = insert_indel(v_aligned, i);
        } else {
            --i;
            --j;
        }
    }
    for (int k = 0; k <= i; ++k)
        w_aligned = insert_indel(w_aligned, 0);
    for (int k = 0; k <= j; ++k)
        v_aligned = insert_indel(v_aligned, 0);

    a << v_aligned << endl;
    a << w_aligned << endl;
}

void find_max_position(const Matrix& s, int& pos_i, int& pos_j) {

    int max = 0, cur = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        auto it = max_element(s[i].begin(), s[i].end());
        cur = distance(s[i].begin(), it);
        if (*it > max) {
            max = *it;
            pos_i = i;
            pos_j = cur;
        }
    }
}
void find_max_position_given_param(const Matrix& s, int& pos_i, int j) {
    int max = s[0][j];
    for (size_t i = 0; i < s.size(); ++i) {
        auto it = s[i][j];
        if (it > max) {
            max = it;
            pos_i = i;
        }
    }
}

void find_max_position_last_row_col(const Matrix& s, int& pos_i, int& pos_j) {
    int max = -1000000;
    for (size_t i = 0; i < s.size(); ++i)
        for (size_t j = 0; j < s[0].size(); ++j)
            if (i == s.size() - 1 || j == s[0].size() - 1)
                if (s[i][j] > max) {
                    max = s[i][j];
                    pos_i = i;
                    pos_j = j;
                }
}
void local_alignment(string& v, string& w, int indel_value) {
    CharMatrix backtrack;
    Matrix s;
    string sym;
    size_t v_size = v.size() + 1, w_size = w.size() + 1;
    int indel = indel_value; // TODO FIXME
    backtrack.resize(v_size);
    s.resize(v_size);
    for (size_t i = 0; i < v_size; ++i) {
        backtrack[i].resize(w_size);
        s[i].resize(w_size);
    }

    for (size_t i = 1; i < v_size; ++i)
        s[i][0] = 0; //-indel * i;

    for (size_t i = 1; i < w_size; ++i)
        s[0][i] = 0; //-indel * i;

    for (size_t i = 1; i < v_size; ++i)
        for (size_t j = 1; j < w_size; ++j) {
            char c1 = v[i - 1], c2 = w[j - 1];
            int match = s[i - 1][j - 1] + PAM250[c1][c2];
            int mismatch = s[i - 1][j - 1] + PAM250[c1][c2];
            int down = s[i - 1][j] - indel;
            int right = s[i][j - 1] - indel;
            int free_jump = 0;

            if (c1 == c2)
                s[i][j] = match;
            else
                s[i][j] = max(max(down, right), max(free_jump, mismatch));

            if (s[i][j] == match)
                sym = "match";
            if (s[i][j] == down)
                sym = "down";
            if (s[i][j] == right)
                sym = "right";
            if (s[i][j] == mismatch)
                sym = "mismatch";
            if (s[i][j] == free_jump)
                sym = "free jump";

            backtrack[i][j] = sym;
        }

    // Get the position of the maximal cell
    int i, j;
    find_max_position(s, i, j);

    int max_value = s[i][j];
    a << max_value << endl;
    string v_aligned = v.substr(0, i), w_aligned = w.substr(0, j);


    while (backtrack[i][j] != "free jump" && i >= 0 && j >= 0) {
        if (backtrack[i][j] == "down") {
            --i;
            w_aligned = insert_indel(w_aligned, j);
        } else if (backtrack[i][j] == "right") {
            --j;
            v_aligned = insert_indel(v_aligned, i);
        } else {
            --i;
            --j;
        }
    }
    v_aligned = v_aligned.substr(i), w_aligned = w_aligned.substr(j);


    a << v_aligned << endl;
    a << w_aligned << endl;
}

void fitting_alignment(string& v, string& w, int indel_value) {
    CharMatrix backtrack;
    Matrix s;
    string sym;
    size_t v_size = v.size() + 1, w_size = w.size() + 1;
    int indel = indel_value; // TODO FIXME
    backtrack.resize(v_size);
    s.resize(v_size);
    for (size_t i = 0; i < v_size; ++i) {
        backtrack[i].resize(w_size);
        s[i].resize(w_size);
    }

    for (size_t i = 1; i < v_size; ++i)
        s[i][0] = 0; //-indel * i;

    for (size_t i = 1; i < w_size; ++i)
        s[0][i] = 0; //-indel * i;

    for (size_t i = 1; i < v_size; ++i)
        for (size_t j = 1; j < w_size; ++j) {
            char c1 = v[i - 1], c2 = w[j - 1];
            int match = s[i - 1][j - 1] + 1;        // PAM250[c1][c2];
            int mismatch = s[i - 1][j - 1] - indel; // PAM250[c1][c2];
            int down = s[i - 1][j] - indel;
            int right = s[i][j - 1] - indel;

            if (c1 == c2)
                s[i][j] = match;
            else
                s[i][j] = max(max(down, right), mismatch);

            if (s[i][j] == match)
                sym = "match";
            if (s[i][j] == down)
                sym = "down";
            if (s[i][j] == right)
                sym = "right";
            if (s[i][j] == mismatch)
                sym = "mismatch";

            backtrack[i][j] = sym;
        }

    int i, j = w.size();
    find_max_position_given_param(s, i, j);

    int max_value = s[i][j];
    a << max_value << endl;

    string v_aligned = v.substr(0, i), w_aligned = w.substr(0, j);

    while (i > 0 && j > 0) {
        if (backtrack[i][j] == "down") {
            --i;
            w_aligned = insert_indel(w_aligned, j);
        } else if (backtrack[i][j] == "right") {
            --j;
            v_aligned = insert_indel(v_aligned, i);
        } else {
            --i;
            --j;
        }
    }
    v_aligned = v_aligned.substr(i);

    a << v_aligned << endl;
    a << w_aligned << endl;
}
void overlap_alignment(string& v, string& w, int indel_value) {
    CharMatrix backtrack;
    Matrix s;
    string sym;
    size_t v_size = v.size() + 1, w_size = w.size() + 1;
    int indel = indel_value; // TODO FIXME
    backtrack.resize(v_size);
    s.resize(v_size);
    for (size_t i = 0; i < v_size; ++i) {
        backtrack[i].resize(w_size);
        s[i].resize(w_size);
    }

    for (size_t i = 1; i < v_size; ++i)
        s[i][0] = 0; //-indel * i;

    for (size_t i = 1; i < w_size; ++i)
        s[0][i] = 0; //-indel * i;

    for (size_t i = 1; i < v_size; ++i)
        for (size_t j = 1; j < w_size; ++j) {
            char c1 = v[i - 1], c2 = w[j - 1];
            int match = s[i - 1][j - 1] + 1;        // PAM250[c1][c2];
            int mismatch = s[i - 1][j - 1] - indel; // PAM250[c1][c2];
            int down = s[i - 1][j] - indel;
            int right = s[i][j - 1] - indel;

            if (c1 == c2)
                s[i][j] = match;
            else
                s[i][j] = max(max(down, right), mismatch);

            if (s[i][j] == match)
                sym = "match";
            if (s[i][j] == down)
                sym = "down";
            if (s[i][j] == right)
                sym = "right";
            if (s[i][j] == mismatch)
                sym = "mismatch";

            backtrack[i][j] = sym;
        }

    int i, j;
    find_max_position_last_row_col(s, i, j);
    int max_value = s[i][j];
    a << max_value << endl;

    string v_aligned = v.substr(0, i), w_aligned = w.substr(0, j);

    while (i > 0 && j > 0) {
        if (backtrack[i][j] == "down") {
            --i;
            w_aligned = insert_indel(w_aligned, j);
        } else if (backtrack[i][j] == "right") {
            --j;
            v_aligned = insert_indel(v_aligned, i);
        } else {
            --i;
            --j;
        }
    }
    v_aligned = v_aligned.substr(i);
    w_aligned = w_aligned.substr(j);

    a << v_aligned << endl;
    a << w_aligned << endl;
}

void global_alignment_affine_gap_penalty(string& v, string& w, int indel_value, int eps) {
    CharMatrix backtrack[3]; // 0 - low, 1 - mid, 2 - upper
    Matrix s[3];
    string sym;
    size_t v_size = v.size() + 1, w_size = w.size() + 1;
    int indel = indel_value; // TODO FIXME
    for (int i = 0; i < 3; ++i) {
        backtrack[i].resize(v_size);
        s[i].resize(v_size);

        for (size_t j = 0; j < v_size; ++j) {
            backtrack[i][j].resize(w_size);
            s[i][j].resize(w_size);
        }
    }
    for (size_t i = 1; i < v_size; ++i) {
        s[0][i][0] = -indel - (i - 1) * eps;
        s[1][i][0] = -indel - (i - 1) * eps;
        s[2][i][0] = -10 * indel;
    }

    for (size_t j = 1; j < w_size; ++j) {
        s[2][0][j] = -indel - (j - 1) * eps;
        s[1][0][j] = -indel - (j - 1) * eps;
        s[0][0][j] = -10 * indel;
    }

    for (size_t i = 1; i < v_size; ++i)
        for (size_t j = 1; j < w_size; ++j) {

            char c1 = v[i - 1], c2 = w[j - 1];
            int lower_scores = max(s[0][i - 1][j] - eps, s[1][i - 1][j] - indel);
            s[0][i][j] = lower_scores;
            if (s[0][i][j] == s[0][i - 1][j] - eps)
                sym = "down";
            if (s[0][i][j] == s[1][i - 1][j] - indel)
                sym = "insert";

            backtrack[0][i][j] = sym;


            int upper_scores = max(s[2][i][j - 1] - eps, s[1][i][j - 1] - indel);
            s[2][i][j] = upper_scores;
            if (s[2][i][j] == s[2][i][j - 1] - eps)
                sym = "right";
            if (s[2][i][j] == s[1][i][j - 1] - indel)
                sym = "delete";

            backtrack[2][i][j] = sym;


            int middle_scores = max(max(s[0][i][j], s[1][i - 1][j - 1] + BLOSUM62[c1][c2]), s[2][i][j]);
            s[1][i][j] = middle_scores;

            if (s[1][i][j] == s[0][i][j])
                sym = "down";
            if (s[1][i][j] == s[1][i - 1][j - 1] + BLOSUM62[c1][c2])
                sym = "match";
            if (s[1][i][j] == s[2][i][j])
                sym = "right";

            backtrack[1][i][j] = sym;
        }


    int i = v.size(), j = w.size();
    vector<int> max_values = {s[0][i][j], s[1][i][j], s[2][i][j]};
    auto it = max_element(max_values.begin(), max_values.end());
    int backtrack_matrix = distance(max_values.begin(), it);
    int max_value = *it;


    string v_aligned = v.substr(0, i), w_aligned = w.substr(0, j);

    while (i > 0 && j > 0) {
        if (backtrack_matrix == 0) {
            if (backtrack[0][i][j] == "insert")
                backtrack_matrix = 1;
            --i;
            w_aligned = insert_indel(w_aligned, j);
        } else if (backtrack_matrix == 1) {

            if (backtrack[1][i][j] == "down")
                backtrack_matrix = 0;
            else if (backtrack[1][i][j] == "right")
                backtrack_matrix = 2;
            else {
                --i;
                --j;
            }

        } else {
            if (backtrack[2][i][j] == "delete")
                backtrack_matrix = 1;
            --j;
            v_aligned = insert_indel(v_aligned, i);
        }
    }
    for (int k = 0; k < i; ++k)
        w_aligned = insert_indel(w_aligned, 0);

    for (int k = 0; k < j; ++k)
        v_aligned = insert_indel(v_aligned, 0);

    a << max_value << endl;
    a << v_aligned << endl;
    a << w_aligned << endl;
}

vector<int> middle_column_score(string& v, string& w, int sigma, vector<int>& backtrack) {
    vector<vector<int>> s;
    s.resize(v.size() + 1);
    for (size_t i = 0; i < v.size() + 1; ++i) {
        s[i] = {0, 0};
        for (ssize_t j = -1; j < 1; ++j)
            s[i][j] = i * j * sigma;
    }
    s[0][1] = -sigma;
    backtrack.resize(v.size() + 1);
    for (size_t i = 0; i < backtrack.size(); ++i)
        backtrack[i] = 0;


    for (size_t j = 1; j < w.size() / 2 + 1; ++j) {
        for (size_t i = 0; i < v.size() + 1; ++i) {
            if (i == 0)
                s[i][1] = -j * sigma;
            else {
                char c1 = v[i - 1], c2 = w[j - 1];

                auto scores = {s[i - 1][0] + BLOSUM62[c1][c2], s[i][0] - sigma, s[i - 1][1] - sigma};
                auto el_it = max_element(scores.begin(), scores.end());
                s[i][1] = *el_it;
                backtrack[i] = distance(scores.begin(), el_it);
            }
        }
        if (j != w.size() / 2) {
            for (size_t i = 0; i < s.size(); ++i)
                s[i][0] = s[i][1];
        }
    }
    vector<int> mid_column;
    mid_column.resize(s.size());
    for (size_t i = 0; i < s.size(); ++i)
        mid_column[i] = s[i][1];


        for (size_t i = 0; i < s.size(); ++i)
            s.pop_back();
    s.clear();

    return mid_column;
}

void middle_edge_in_linear_space_problem(string& v, string& w, int indel_value) {
    vector<int> backtrack;
    vector<int> source_to_middle(middle_column_score(v, w, indel_value, backtrack));

    string v_rev(v.rbegin(), v.rend());
    string w_rev(w.rbegin(), w.rend());

    auto middle_to_sink = middle_column_score(v_rev, w_rev, indel_value, backtrack);

    int max = -1000000;
    size_t pos = 0;
    for (size_t i = 0; i < source_to_middle.size(); ++i) {
        int val = source_to_middle[i] + middle_to_sink[i];
        if (val > max) {
            max = val;
            pos = i;
        }
    }
    cout << "(" << pos << "," << (w.size() / 2) << endl;

    if (pos == source_to_middle.size() - 1)
        cout << "(" << pos << "," << (w.size() / 2) + 1 << endl;
    else {
        if (backtrack[pos] == 0)
            cout << "(" << pos + 1 << "," << (w.size() / 2) + 1 << endl;
        else if (backtrack[pos] == 1)
            cout << "(" << pos << "," << (w.size() / 2) + 1 << endl;
        else
            cout << "(" << pos << "," << (w.size() / 2) << endl;
    }
}

int main() {
    ifstream f("data.txt");

    // Graph g;

    // g.read_graph_from_file_by_edges(f);
    // g.longest_path(a);
    // cout<<g;
    // g.print_topo_order(a);

    string first, second;
    f >> first >> second;
    // string first = "PLEASANTLY", second = "MEANLY";
    // cout << first << endl;
    //        for (auto row : backtrack) {
    //            for (auto i : row)
    //                cout << i << " ";
    //            cout << endl;
    //}

    int indel = 5;
    int eps = 1;
    int arr[2];
    arr[-1] = 2;

    middle_edge_in_linear_space_problem(first, second, indel);
    // local_alignment(first, second, indel);

    //    auto backtrack =  LCS_backtrack(first,second,max_val,indel);
    //    cout << abs(max_val) << endl;
    //   output_global_alignment(backtrack, first, second);

    // int n, m;
    //    Matrix down, right;
    //    read_input(n, m, down, right, f);
    // cout << manhattan_tourist_problem(n, m, down, right) << endl;

    return 0;
}
