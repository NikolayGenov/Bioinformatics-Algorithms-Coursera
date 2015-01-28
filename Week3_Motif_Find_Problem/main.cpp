#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <set>
#include <iterator>
#include <stdlib.h> /* srand, rand */
#include <time.h>
#include <math.h>
#include <map>
#include <random>
#include <algorithm>
using namespace std;

static map<char, int> letters = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};
static map<int, char> lettersInverse = {{0, 'A'}, {1, 'C'}, {2, 'G'}, {3, 'T'}};

ofstream answer("answer.txt");
static vector<double> powmap;

static void calculatePows(int n) {
    powmap.resize(n);
    for (int i = 0; i < n; ++i)
        powmap[i] = pow(4, i);
}
template <class T> void transpose(std::vector<std::vector<T>>& a, std::vector<std::vector<T>>& b, int width, int height) {
    for (int i = 0; i < width; i++)
        for (int j = 0; j < height; j++)
            b[j][i] = a[i][j];
}

int hammingDistance(const string& A, const string& B) {
    int distance = 0;
    for (size_t i = 0; i < A.length(); ++i) {
        if (A[i] != B[i])
            distance++;
    }
    return distance;
}

int PatternToNumber(string pattern) {
    if (pattern.empty())
        return 0;
    int letc = letters[pattern.back()];
    pattern.pop_back();
    return 4 * PatternToNumber(pattern) + letc;
}

long long patternToNumber(const string& pattern) {
    long long sum = 0;
    for (int i = pattern.length() - 1, j = 0; i >= 0; --i, ++j) {
        sum += letters[pattern[j]] * powmap[i];
    }
    return sum;
}


string NumberToPattern(int index, int k) {
    if (k == 1)
        return string(1, lettersInverse[index]);
    int prefixIndex = index / 4;
    int r = index % 4;
    return NumberToPattern(prefixIndex, k - 1) + lettersInverse[r];
}

string numbersToPattern(int input, int len) {
    int sum = input, digit = 0, powof4;

    string pattern(len, ' ');
    for (int i = len - 1; i >= 0; --i) {
        powof4 = pow(4, i);
        digit = 0;
        for (int d = 0; d <= 3; ++d) {
            if (d * powof4 > sum)
                break;
            else
                digit = d;
        }
        pattern[len - 1 - i] = lettersInverse[digit];
        sum -= digit * powof4;
    }
    return pattern;
}

void computingFrequencies(vector<int>& freqs, const string& input, int length) {
    int lok = 0;
    for (size_t i = 0; i < input.size() - length + 1; ++i) {
        lok = patternToNumber(input.substr(i, length));
        freqs[lok]++;
    }
}

set<string> neighbors(string pattern, int d) {
    if (d == 0)
        return {pattern};
    if (pattern.length() == 1)
        return {"A", "C", "G", "T"};
    set<string> Neighborhood;
    string suffix = pattern.substr(1);
    set<string> SuffixNeighbors = neighbors(suffix, d);
    for (string text : SuffixNeighbors) {
        if (hammingDistance(suffix, text) < d)
            for (int i = 0; i < 4; ++i)
                Neighborhood.insert(lettersInverse[i] + text);
        else
            Neighborhood.insert(pattern[0] + text);
    }
    return Neighborhood;
}


bool can_find_pattern_in_DNAs(vector<string>& DNAs, const string& pattern, int k, int d) {
    bool can_be_found_in_one_seq = true;
    size_t i;
    for (i = 0; i < DNAs.size() && can_be_found_in_one_seq; ++i) {
        can_be_found_in_one_seq = false;
        for (size_t j = 0; j < DNAs[i].size() - k + 1; ++j)
            if (hammingDistance(pattern, DNAs[i].substr(j, k)) <= d)
                can_be_found_in_one_seq = true;
    }
    if (i == DNAs.size() && can_be_found_in_one_seq)
        return true;

    return false;
}

set<string> motif_enumeration(vector<string>& DNAs, int k, int d) {
    set<string> patterns = {};
    string pattern;
    set<string> Neighborhood;
    int four_power_k = pow(4, k);
    vector<bool> close(four_power_k, false);

    for (size_t i = 0; i < DNAs.size(); ++i) {
        string seq = DNAs[i];
        for (size_t i = 0; i < seq.size() - k + 1; ++i) {
            Neighborhood = neighbors(seq.substr(i, k), d);
            for (string pattern : Neighborhood)
                close[patternToNumber(pattern)] = true;
        }
        for (int i = 0; i < four_power_k; ++i) {
            if (close[i]) {
                pattern = numbersToPattern(i, k);
                if (can_find_pattern_in_DNAs(DNAs, pattern, k, d))
                    patterns.insert(pattern);
            }
        }
    }

    return patterns;
}

int distance_between_pattern_and_strings(string pattern, const vector<string>& DNA) {
    size_t k = pattern.size();
    int distance = 0;
    for (auto text : DNA) {
        int hamming_distance = 1 << 20;

        for (size_t j = 0; j < text.size() - k + 1; ++j) {
            int dist = hammingDistance(pattern, text.substr(j, k));
            if (hamming_distance > dist)
                hamming_distance = dist;
        }

        distance += hamming_distance;
    }
    return distance;
}

string median_string(const vector<string>& DNA, int k) {
    int distance = 1 << 20;
    string median;
    for (int i = 0; i < powmap[k]; ++i) {
        string pattern = NumberToPattern(i, k);
        int current_dist = distance_between_pattern_and_strings(pattern, DNA);
        if (distance > current_dist) {
            distance = current_dist;
            // Question 6 and add >= up there ^^^
            // cout<<pattern<< "   "<< distance<<endl;

            median = pattern;
        }
    }
    return median;
}

double compute_profile_probability(const string& k_mer, const vector<vector<double>>& matrix) {
    double probability = 1;
    for (size_t i = 0; i < k_mer.size(); ++i) {
        probability *= matrix[letters[k_mer[i]]][i];
        if (probability == 0)
            break;
    }
    return probability;
}

string profile_most_probable_k_mer(string DNA, int k, const vector<vector<double>>& matrix) {
    double probability = 0;
    string most_probable_k_mer;

    for (size_t i = 0; i < DNA.size() - k + 1; ++i) {
        string k_mer = DNA.substr(i, k);
        double prob = compute_profile_probability(k_mer, matrix);
        if (prob > probability) {
            probability = prob;
            most_probable_k_mer = k_mer;
        }
    }
    if (probability == 0)
        return DNA.substr(0, k);
    return most_probable_k_mer;
}

vector<vector<double>> create_profile_matrix(const vector<string>& motifs, int k) {
    double addent = 1.0 / (2 * motifs.size()); // WITH pseudocounts
    vector<vector<double>> matrix(letters.size(), vector<double>(k, addent));

    // vector<vector<double>> matrix(letters.size(), vector<double>(k, 0));
    // double addent = 1.0 / motifs.size(); //For WITHOUT pseudocounts
    for (size_t i = 0; i < motifs.size(); ++i) {
        for (int j = 0; j < k; ++j) {
            matrix[letters[motifs[i][j]]][j] += addent;
        }
    }
    return matrix;
}

string find_consensus(const vector<string>& motifs, int k) {
    string consensus = "";
    vector<vector<double>> matrix = create_profile_matrix(motifs, k);
    vector<vector<double>> trans_matrix(k, vector<double>(letters.size(), 0));

    transpose(matrix, trans_matrix, letters.size(), k);

    for (int i = 0; i < k; ++i) {
        size_t index = std::distance(trans_matrix[i].begin(), max_element(trans_matrix[i].begin(), trans_matrix[i].end()));
        consensus += lettersInverse[index];
    }
    return consensus;
}

int score_motifs(const vector<string>& motifs, int k) {
    string consensus = find_consensus(motifs, k);
    int score = 0;
    for (string motif : motifs)
        score += hammingDistance(consensus, motif);
    return score;
}

vector<string> greedy_motif_search(const vector<string>& DNA, int k, int t) {

    vector<string> best_motifs;
    for (string strand : DNA)
        best_motifs.push_back(strand.substr(0, k));

    string base_strand = DNA[0];

    for (size_t i = 0; i < base_strand.size() - k + 1; ++i) {
        string k_mer = base_strand.substr(i, k), next_motif;

        vector<string> motifs = {k_mer};
        for (int j = 1; j < t; ++j) {
            vector<vector<double>> matrix = create_profile_matrix(motifs, k);
            next_motif = profile_most_probable_k_mer(DNA[j], k, matrix);
            motifs.push_back(next_motif);
        }
        int score_mot = score_motifs(motifs, k);
        int score_best_mot = score_motifs(best_motifs, k);
        if (score_mot < score_best_mot)
            best_motifs = motifs;
    }
    return best_motifs;
}

// Random algos

vector<string> motifs_for_random(const vector<string>& DNA, const vector<vector<double>>& matrix, int k) {

    vector<string> motifs;
    for (string strand : DNA) {
        motifs.push_back(profile_most_probable_k_mer(strand, k, matrix));
    }

    return motifs;
}

pair<vector<string>, int> randomized_motif_search(const vector<string>& DNA, int k) {

    int max_pos = DNA[0].size() - k + 1;
    vector<string> best_motifs, motifs;
    for (string strand : DNA) {
        int rand_pos = rand() % max_pos;
        best_motifs.push_back(strand.substr(rand_pos, k));
    }

    motifs = best_motifs;

    while (true) {
        vector<vector<double>> profile = create_profile_matrix(motifs, k);

        motifs = motifs_for_random(DNA, profile, k);

        int score_mot = score_motifs(motifs, k);
        int score_best_mot = score_motifs(best_motifs, k);
        if (score_mot < score_best_mot)
            best_motifs = motifs;
        else
            return make_pair(best_motifs, score_best_mot);
    }
}


string Profile_randomly_generated_k_mer(string motif_i, const vector<vector<double>>& matrix, int k) {
    std::random_device rd;
    std::mt19937 gen(rd());
    vector<double> prob;
    for (size_t i = 0; i < motif_i.size() - k + 1; ++i) {
        string k_mer = motif_i.substr(i, k);
        prob.push_back(compute_profile_probability(k_mer, matrix));
    }

    std::discrete_distribution<> random(prob.begin(), prob.end());
    int pos = random(gen);
    return motif_i.substr(pos, k);
}

pair<vector<string>, int> gibbs_sampler_motif_search(const vector<string>& DNA, int k, int t, int N) {
    int max_pos = DNA[0].size() - k + 1;
    int score_best_mot;
    vector<string> best_motifs, motifs;
    for (string strand : DNA) {
        int rand_pos = rand() % max_pos;
        best_motifs.push_back(strand.substr(rand_pos, k));
    }
    motifs = best_motifs;

    for (int j = 0; j < N; ++j) {

        int i = rand() % t;

        string row_i = DNA[i];
        motifs.erase(motifs.begin() + i);
        vector<vector<double>> profile = create_profile_matrix(motifs, k);

        string motif_i = Profile_randomly_generated_k_mer(row_i, profile, k);
        motifs.insert(motifs.begin() + i, motif_i);

        int score_mot = score_motifs(motifs, k);
        score_best_mot = score_motifs(best_motifs, k);
        if (score_mot < score_best_mot)
            best_motifs = motifs;
    }
    return make_pair(best_motifs, score_best_mot);
}

vector<string> run_random_motif_search_times(const vector<string>& DNA, int k, int t, int N, int times) {
    pair<vector<string>, int> motif_pair, best_motif_pair;
    best_motif_pair = randomized_motif_search(DNA, k);

    for (int i = 0; i < times; ++i) {

        motif_pair = gibbs_sampler_motif_search(DNA, k, t, N);

        if (best_motif_pair.second > motif_pair.second)
            best_motif_pair = motif_pair;
    }

    cout << "Best rusult is : " << best_motif_pair.second << endl;

    return best_motif_pair.first;
}


/*
vector<vector<double>> entr = {
{0.5, 0, 0, 0.5},
{0.25, 0.25, 0.25, 0.25},
{0, 0, 0, 1},
{0.25, 0, 0.5, 0.25}
        };

double term(double f, double s) {
    if (f == 0)
        return 0;
    return f * s;
}
void entropy (){
    double sum = 0;
    for (size_t i = 0; i < entr.size(); ++i) {
        sum = 0;
        for (size_t j = 0; j < entr[0].size(); ++j) {
            sum += term(entr[i][j], log2(entr[i][j]));
        }
        sum = -sum;
        cout<<sum<<endl;
    }

}
*/

int main() {
    //// Precalculations
    calculatePows(200);
    srand(time(NULL));

    // Input data
    const char* name = "data.txt";
    ifstream f(name);
    // int k, t,N;
    vector<string> DNA;
    // f >> k >> t>>N;
    //    // k = 15,d = 4;

    string s; //,pat;
    //    // f>>pat;
    while (f >> s)
        DNA.push_back(s);

    int k = 15;


    auto set_patterns = run_random_motif_search_times(DNA, k, DNA.size(), 2000, 30);

    cout << find_consensus(set_patterns, k) << endl;
    // Print Answer

    //        for (auto s : set_patterns)
    //            answer << s << endl;
    //        answer<<endl;

    //// Real calculations
    // auto set_patterns = motif_enumeration(DNAs, k, d);

    //// Print Answer
    // for (auto s : set_patterns)
    // answer << s << " ";


    //    auto best_motifs = greedy_motif_search(DNA, k, t);

    //    // Print Answer
    //    for (auto s : best_motifs)
    //        answer << s << endl;
    //    answer << endl;

    // answer<<profile_most_probable_k_mer(DNA, k, matrix);
    //  answer<<median_string(DNA, k)<<endl;

    //    //Task 3 - entropy
    //    entropy();

    return 0;
}
