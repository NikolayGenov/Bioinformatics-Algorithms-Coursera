#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <set>
#include <math.h>
#include <queue>
#include "protein_table.h"
using namespace std;

ofstream answer("answer");
const int CODON_LENGTH = 3;

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

string translate_DNA_to_RNA(const string& DNA) {
    string RNA(DNA);
    for (size_t i = 0; i < DNA.length(); ++i)
        if (DNA[i] == 'T')
            RNA[i] = 'U';
    return RNA;
}

string translate_RNA_to_Protein(const string& str) {
    string protein(str.length() / CODON_LENGTH, ' ');
    for (size_t i = 0, j = 0; i < str.length(); i += CODON_LENGTH, ++j)
        protein[j] = protein_map[str.substr(i, CODON_LENGTH)];
    return protein;
}

vector<string> substring_of_genome_and_complement_matching_aminos(const string& DNA, const string& aminos) {
    size_t amino_to_DNA_length = aminos.length() * CODON_LENGTH;
    vector<string> valid_translations_of_DNA_to_the_amino;
    for (size_t i = 0; i < DNA.length(); ++i) {
        string substrDNA = DNA.substr(i, amino_to_DNA_length);
        string revDNA = reverse_complement(substrDNA);
        string trans_aminos = translate_RNA_to_Protein(translate_DNA_to_RNA(substrDNA));
        string reverse_trans_aminos = translate_RNA_to_Protein(translate_DNA_to_RNA(revDNA));
        if (aminos == trans_aminos || aminos == reverse_trans_aminos)
            valid_translations_of_DNA_to_the_amino.push_back(substrDNA);
    }
    return valid_translations_of_DNA_to_the_amino;
}

vector<int> cyclic_spectrum(string peptide) {
    vector<int> prefix_mass(peptide.length() + 1, 0);
    for (size_t i = 1; i <= peptide.length(); ++i) {
        prefix_mass[i] = prefix_mass[i - 1] + integer_mass_table[peptide[i - 1]];
    }

    int peptide_mass = prefix_mass[peptide.length()];
    vector<int> cyclic_spectrum_list(1, 0);
    for (size_t i = 0; i < peptide.length(); ++i)
        for (size_t j = i + 1; j <= peptide.length(); ++j) {
            int i_j = prefix_mass[j] - prefix_mass[i];
            int cyc_i_j = peptide_mass - i_j;
            cyclic_spectrum_list.push_back(i_j);
            if (i != 0 && j != peptide.length())
                cyclic_spectrum_list.push_back(cyc_i_j);
        }

    std::sort(cyclic_spectrum_list.begin(), cyclic_spectrum_list.end());
    return cyclic_spectrum_list;
}

vector<int> linear_spectrum(string peptide) {
    vector<int> prefix_mass(peptide.length() + 1, 0);
    for (size_t i = 1; i <= peptide.length(); ++i) {
        prefix_mass[i] = prefix_mass[i - 1] + integer_mass_table[peptide[i - 1]];
    }

    vector<int> linear_spectrum_list(1, 0);
    for (size_t i = 0; i < peptide.length(); ++i)
        for (size_t j = i + 1; j <= peptide.length(); ++j)
            linear_spectrum_list.push_back(prefix_mass[j] - prefix_mass[i]);

    std::sort(linear_spectrum_list.begin(), linear_spectrum_list.end());
    return linear_spectrum_list;
}

vector<int> cyclic_spectrum_with_ints(vector<int>& peptide) {
    vector<int> prefix_mass(peptide.size() + 1, 0);

    for (size_t i = 1; i <= peptide.size(); ++i) {
        prefix_mass[i] = prefix_mass[i - 1] + peptide[i - 1];
    }

    int peptide_mass = prefix_mass[peptide.size()];
    vector<int> cyclic_spectrum_list(1, 0);
    for (size_t i = 0; i < peptide.size(); ++i)
        for (size_t j = i + 1; j <= peptide.size(); ++j) {
            int i_j = prefix_mass[j] - prefix_mass[i];
            int cyc_i_j = peptide_mass - i_j;
            cyclic_spectrum_list.push_back(i_j);
            if (i != 0 && j != peptide.size()) {
                cyclic_spectrum_list.push_back(cyc_i_j);
            }
        }

    std::sort(cyclic_spectrum_list.begin(), cyclic_spectrum_list.end());
    return cyclic_spectrum_list;
}

vector<int> linear_spectrum_with_ints(vector<int>& peptide) {
    vector<int> prefix_mass(peptide.size() + 1, 0);
    for (size_t i = 1; i <= peptide.size(); ++i) {
        prefix_mass[i] = prefix_mass[i - 1] + peptide[i - 1];
    }

    vector<int> linear_spectrum_list(1, 0);
    for (size_t i = 0; i < peptide.size(); ++i)
        for (size_t j = i + 1; j <= peptide.size(); ++j)
            linear_spectrum_list.push_back(prefix_mass[j] - prefix_mass[i]);

    std::sort(linear_spectrum_list.begin(), linear_spectrum_list.end());
    return linear_spectrum_list;
}

void expand_peptides(vector<vector<int>>& peptides, const vector<int>& amino_acid_masses) {
    vector<vector<int>> new_peptides;
    for (auto peptide = peptides.begin(); peptide != peptides.end(); ++peptide)
        for (size_t i = 0; i < amino_acid_masses.size(); ++i) {
            vector<int> copy_vec(*peptide);
            copy_vec.push_back(amino_acid_masses[i]);
            new_peptides.push_back(copy_vec);
        }
    peptides = new_peptides;
}

bool is_consistent_with_spectrum(const vector<int>& peptide, map<int, int>& spectrum_map) {
    bool is_consistent = true;

    std::map<int, int> peptide_map;
    for (size_t i = 0; i < peptide.size(); ++i)
        ++peptide_map[peptide[i]];

    for (auto it = peptide_map.begin(); it != peptide_map.end(); ++it) {
        if (it->second > spectrum_map[it->first]) {
            is_consistent = false;
            break;
        }
    }
    return is_consistent;
}

bool is_consistent_with_spectrum_given_string(string& peptide, map<int, int>& spectrum_map) {
    bool is_consistent = true;

    std::map<int, int> peptide_map;
    for (size_t i = 0; i < peptide.size(); ++i)
        ++peptide_map[integer_mass_table[peptide[i]]];

    for (auto it = peptide_map.begin(); it != peptide_map.end(); ++it) {
        if (it->second > spectrum_map[it->first]) {
            is_consistent = false;
            break;
        }
    }
    return is_consistent;
}

vector<vector<int>> cyclo_peptide_squencing(vector<int> spectrum) {
    int parent_sum_spectrum = spectrum[spectrum.size() - 1];
    std::map<int, int> spectrum_map;

    spectrum.erase(spectrum.begin()); // because we add zero in the linear_spectrum
    vector<int> linear_spec_spec = linear_spectrum_with_ints(spectrum);
    for (size_t i = 0; i < linear_spec_spec.size(); ++i) {
        ++spectrum_map[linear_spec_spec[i]];
    }
    spectrum.insert(spectrum.begin(), 0); // Readd it again

    bool iter = true;
    vector<vector<int>> peptides = {{}}, result; // Must be empty and not zero

    while (!peptides.empty()) {
        expand_peptides(peptides, integer_table);

        for (auto peptide = peptides.begin(); peptide != peptides.end();) {
            int sum_peptide = std::accumulate(peptide->begin(), peptide->end(), 0);

            if (parent_sum_spectrum == sum_peptide) {

                vector<int> cyclic_spec = cyclic_spectrum_with_ints(*peptide);
                if (cyclic_spec == spectrum)
                    result.push_back(*peptide);
                peptide = peptides.erase(peptide);
                iter = false;
            } else if (!(is_consistent_with_spectrum(linear_spectrum_with_ints(*peptide), spectrum_map))) {
                peptide = peptides.erase(peptide);
                iter = false;
            }
            if (iter)
                ++peptide;
            iter = true;
        }

        set<vector<int>> st(peptides.begin(), peptides.end());
        peptides.clear();
        std::copy(st.begin(), st.end(), std::back_inserter(peptides));
    }

    return result;
}

int scoring_problem(const vector<int>& spectrum, const vector<int> experimental) {
    bool end = true;
    int count_equal = 0;
    for (size_t i = 0, j = 0; end;) {
        if (spectrum[i] == experimental[j]) {
            ++count_equal, ++i, ++j;
        } else if (spectrum[i] < experimental[j]) {
            ++i;
        } else {
            ++j;
        }
        if (i == spectrum.size() || j == experimental.size())
            break;
    }
    return count_equal;
}
int cyclopeptide_scoring_problem_given_string(const string& peptide, const vector<int> experimental) {
    return scoring_problem(cyclic_spectrum(peptide), experimental);
}
int linear_scoring_problem_given_string(const string& peptide, const vector<int> experimental) {
    return scoring_problem(linear_spectrum(peptide), experimental);
}
int linear_scoring_problem_given_ints(vector<int>& peptide, const vector<int> experimental) {
    return scoring_problem(linear_spectrum_with_ints(peptide), experimental);
}
bool compare_linear_peptides_by_int(pair<int, vector<int>>& a, pair<int, vector<int>>& b) {
    return a.first < b.first;
}
void trim(vector<vector<int>>& leaderboard, const vector<int>& spectrum, size_t n) {

    priority_queue<pair<int, vector<int>>, vector<pair<int, vector<int>>>, std::function<bool(pair<int, vector<int>>&, pair<int, vector<int>>&)>> pq(
    compare_linear_peptides_by_int);

    for (auto peptide = leaderboard.begin(); peptide != leaderboard.end(); ++peptide) {
        pq.push(make_pair(linear_scoring_problem_given_ints(*peptide, spectrum), *peptide));
    }
    size_t size = leaderboard.size();
    leaderboard.clear();

    pair<int, vector<int>> n_th_peptide;
    for (size_t i = 0; i < size; ++i) {
        pair<int, vector<int>> p = pq.top();
        if (i < n) {
            if (i == n - 1) {
                n_th_peptide = p;
            }
        } else if (p.first < n_th_peptide.first) {
            break;
        }
        leaderboard.push_back(p.second);
        pq.pop();
    }
}
bool compare_linear_peptides_by_int_given_string(pair<int, string>& a, pair<int, string>& b) {
    return a.first < b.first;
}
void trim_given_string(vector<string>& leaderboard, const vector<int>& spectrum, int n) {

    priority_queue<pair<int, string>, vector<pair<int, string>>, std::function<bool(pair<int, string>&, pair<int, string>&)>> pq(compare_linear_peptides_by_int_given_string);

    for (auto peptide = leaderboard.begin(); peptide != leaderboard.end(); ++peptide) {
        pq.push(make_pair(linear_scoring_problem_given_string(*peptide, spectrum), *peptide));
    }
    int size = leaderboard.size();
    leaderboard.clear();

    pair<int, string> n_th_peptide;
    for (int i = 0; i < size; ++i) {
        pair<int, string> p = pq.top();
        if (i < n) {
            if (i == n - 1) {
                n_th_peptide = p;
            }
        } else if (p.first < n_th_peptide.first) {
            break;
        }
        leaderboard.push_back(p.second);
        pq.pop();
    }
}
string format_answer(const vector<int>& peptide) {
    string pep_str;
    for (size_t i = 0; i < peptide.size(); ++i) {
        pep_str += to_string(peptide[i]);
        pep_str += "-";
    }
    pep_str.erase(pep_str.end() - 1);
    return pep_str;
}
string format_answer_to_chars(vector<int>& peptide) {
    string pep_str;
    for (size_t i = 0; i < peptide.size(); ++i) {
        pep_str += mass_to_char_table[peptide[i]];
    }
    return pep_str;
}

vector<int> leaderboard_cyclopeptide_squencing(vector<int> spectrum, int n, const vector<int>& amino_masses) {
    int parent_sum_spectrum = spectrum[spectrum.size() - 1];

    bool iter = true;
    vector<vector<int>> leaderboard = {{}}, old_leaderboard;
    vector<int> leader_peptide = {};
    int score_leader_peptide = 0;
    //   int count = 0;
    while (!leaderboard.empty()) {
        old_leaderboard.assign(leaderboard.begin(), leaderboard.end());
        expand_peptides(leaderboard, amino_masses);

        for (auto peptide = leaderboard.begin(); peptide != leaderboard.end();) {
            int sum_peptide = std::accumulate(peptide->begin(), peptide->end(), 0);

            if (parent_sum_spectrum == sum_peptide) {

                vector<int> cyclic_spec = cyclic_spectrum_with_ints(*peptide);
                int score_peptide = scoring_problem(cyclic_spec, spectrum);
                //                if (score_peptide == 82){
                //                    answer<<format_answer(*peptide)<<" ";
                //                }
                if (score_peptide > score_leader_peptide) {
                    // cout<<format_answer(*peptide)<<" ";
                    // cout<<score_peptide<<endl;
                    leader_peptide = *peptide;
                    score_leader_peptide = score_peptide;
                }
            } else if (sum_peptide > parent_sum_spectrum) {
                peptide = leaderboard.erase(peptide);
                iter = false;
            }

            if (iter)
                ++peptide;
            iter = true;
        }

        set<vector<int>> st(leaderboard.begin(), leaderboard.end());
        leaderboard.clear();
        std::copy(st.begin(), st.end(), std::back_inserter(leaderboard));

        trim(leaderboard, spectrum, n);
    }

    return leader_peptide;
}

vector<int> trim_spectral_convolution(const map<int, int>& convolution_map, int n) {
    auto cmp = [](std::pair<int, int> const& a, std::pair<int, int> const& b) { return a.second > b.second; };

    vector<pair<int, int>> v;
    copy(convolution_map.begin(), convolution_map.end(), back_inserter(v));
    std::sort(v.begin(), v.end(), cmp);

    vector<int> result_convolution;
    int i = 0;
    pair<int, int> n_th_value;
    for (auto a : v) {
        pair<int, int> p = a;
        if (i < n) {
            if (i == n - 1) {
                n_th_value = p;
            }
        } else if (p.second < n_th_value.second) {
            return result_convolution;
        }
        if (a.first >= 57 && a.first <= 200) {
            for (int i = 0; i < a.second; ++i) {
                result_convolution.push_back(a.first);
            }
            ++i;
        }
    }
    return result_convolution;
}

vector<int> spectral_convolution_problem(vector<int>& spectrum, int N) {
    map<int, int> convolution_map;
    sort(spectrum.begin(), spectrum.end());

    for (size_t i = 0; i < spectrum.size(); ++i)
        for (size_t j = i + 1; j < spectrum.size(); ++j)
            convolution_map[spectrum[j] - spectrum[i]]++;
    N = N == 0 ? convolution_map.size() : N;
    return trim_spectral_convolution(convolution_map, N);
}

vector<int> convolution_cyclopeptide_sequencing(vector<int> spectrum, int M, int N) {
    vector<int> amino_masses_remaining = spectral_convolution_problem(spectrum, M);
    set<int> set_amino_masses(amino_masses_remaining.begin(), amino_masses_remaining.end());
    amino_masses_remaining.clear();
    std::copy(set_amino_masses.begin(), set_amino_masses.end(), std::back_inserter(amino_masses_remaining));

    return leaderboard_cyclopeptide_squencing(spectrum, N, amino_masses_remaining);
}

void print_to_file(const vector<string>& vec) {
    for (auto str : vec)
        answer << str << " ";
    answer << endl;
}

void print_to_file_int(const vector<int>& vec) {
    for (auto str : vec)
        answer << str << " ";
    answer << endl;
}

int main() {
    ifstream data_file("data.txt");
    string DNA, aminos, peptide;
    // data_file>>DNA>>aminos;
    int num = 0;
    int N, M;
    vector<int> input;
    data_file >> M >> N;
    while (data_file >> num)
        input.push_back(num);

    // ConvolutionCyclopeptideSequencing
    //     vector<int> result = convolution_cyclopeptide_sequencing(input,M,N);
    //   answer<<format_answer(result)<<endl;

    /* Leaderboard
    // vector<int> leader =
    //    leaderboard_cyclopeptide_squencing(input,1000);


    //answer<<format_answer(leader)<<endl;

    // vector<int> pep = { 0, 71, 87 ,101 ,113 ,158 ,184 ,188 ,259 ,271 ,372};

    //    cout<<linear_scoring_problem_given_ints(pep,input)<<endl;
    // cout<<linear_scoring_problem_given_ints(leader,input)<<endl;
    */

    // Question 5
    /*
       vector<vector<int>> res = cyclo_peptide_squencing({0,71 ,101,113 ,131 ,184
       ,202 ,214 ,232 ,285 ,303 ,315 ,345 ,416});
       for (vector<int> vec : res)
       cout << format_answer_to_chars(vec) << " ";
       */

    /* //Question 6
       int result = cyclopeptide_scoring_problem_given_string("MAMA",{0 ,57 ,71
       ,71
       ,71 ,104 ,131 ,202 ,202 ,202 ,256 ,333 ,333 ,403, 404});
       cout<<result<<endl;
       */
    /* //Question 7
       int result = linear_scoring_problem_given_string("PEEP",{0 ,97 ,97 ,129
       ,129
       ,194 ,203 ,226 ,226 ,258 ,323 ,323 ,323 ,355 ,403 ,452});
       cout<<result<<endl;
       */

    // Question 8

    /*
       vector<int> spectrum = {0 ,71 ,99 ,101 ,103 ,128 ,129 ,199 ,200 ,204 ,227
       ,230 ,231 ,298 ,303 ,328 ,330 ,332 ,333};
       std::map<int, int> spectrum_map;

       spectrum.erase(spectrum.begin());       // because we add zero in the
       linear_spectrum
       vector<int> linear_spec_spec = linear_spectrum_with_ints(spectrum);
       for (size_t i = 0; i < linear_spec_spec.size(); ++i){
       ++spectrum_map[linear_spec_spec[i]];
       }

       string str = "TCE";
       vector<int> str_spec = linear_spectrum(str);
       cout<<is_consistent_with_spectrum (str_spec, spectrum_map)<<endl;

    // Question 9

    vector<int> res = {0 ,86 ,160 ,234 ,308 ,320 ,382};
    print_to_file_int(spectral_convolution_problem(res,0));
    */

    /*
     * //Trim given string
    //    string str;
    //     vector<string> strs;
    //    while ((data_file.peek() != '\n') && (data_file>>str))
    strs.push_back(str);
    //        data_file.ignore();
    //     while ((data_file.peek()!='\n') && (data_file>>num))
    input.push_back(num);
    //     data_file>>n;
    //   trim_given_string(strs,input,n);

    //    for(auto a: strs)
    //        answer<<a<<" ";
    //    answer<<endl;

  */

    cout << "Done!" << endl;
    return 0;
}
