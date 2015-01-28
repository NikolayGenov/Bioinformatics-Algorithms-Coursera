#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <algorithm>
#include <set>
#include <cmath>

using namespace std;

static ofstream answer("answer.txt");
string readData(const char* filename){
    ifstream f(filename);
    string input,tmp;
    while(f>>tmp)
        input+=tmp;
    return input;
}

static map<char, int> letters = {{'A',0},{'C',1}, {'G',2}, {'T',3} };
static map<int, char> lettersInverse = {{ 0,'A'},{ 1,'C'}, { 2,'G'}, { 3,'T'} };

vector<int> minSkewCounter(string input){
    int skew = 0,minSkew = 0;
    vector<int> minPos;
    for (size_t i = 0; i < input.length(); ++i) {
        if (input[i] == 'G') skew++;
        if (input[i] == 'C') skew--;
        if(skew == minSkew)  minPos.push_back(i+1);
        if(skew < minSkew){
            minPos.clear();
            minPos.push_back(i+1);
            minSkew = skew;
        }
    }
    return minPos;
}

string makeReverse(string input){
    string rev(input.size(),' ');
    for(int i = input.size()-1,j=0 ;i>=0; --i,++j)
        switch (input[i]) {
        case 'A':rev[j] = 'T'; break;
        case 'T':rev[j] = 'A'; break;
        case 'C':rev[j] = 'G'; break;
        case 'G':rev[j] = 'C'; break;
        default:
            break;
        }
    return rev;
}

long long  PatternToNumber(string pattern){
    if(pattern.empty()) return 0;
    int letc = letters[pattern.back()];
    pattern.pop_back();
    return 4 * PatternToNumber(pattern) + letc;
}

string NumberToPattern(int index,int  k){
    if (k == 1) return string(1,lettersInverse[index]);
    int prefixIndex = index /  4;
    int r = index % 4;
    return NumberToPattern(prefixIndex, k - 1) +  lettersInverse[r] ;

}

void computingFrequencies(vector<int>& freqs, const string& input,int length){
    int lok = 0;
    for(size_t i = 0 ; i< input.size() - length + 1; ++i){
        lok = PatternToNumber(input.substr(i,length));
        freqs[lok]++;
    }
}

int hammingDistance(const string& A,const string& B){
    int distance = 0;
    for (size_t i = 0; i < A.length(); ++i)
        if(A[i] != B[i])
            distance++;
    return distance;
}

int ApproximatePatternCount(string input,string pattern, int d){

    int hammingCount = 0;
    string text;
    int len = pattern.length();
    for(size_t i = 0 ; i< input.size() - len + 1; ++i){
        text = input.substr(i,len);
        if(hammingDistance(pattern,text) <= d)
            hammingCount++;
    }
    return hammingCount;
}

set<string> neighbors(string pattern, int d){
    if (d == 0)
        return {pattern};
    if (pattern.length() == 1)
        return {"A", "C", "G", "T"};
    set<string> Neighborhood;
    string suffix = pattern.substr(1);
    set<string> SuffixNeighbors = neighbors(suffix, d);
    for(string text : SuffixNeighbors){
        if (hammingDistance(suffix,text) < d)
            for(int i =0; i< 4;++i)
                Neighborhood.insert( lettersInverse[i] + text);
        else
            Neighborhood.insert(pattern[0] + text);
    }
    return Neighborhood;
}


set<string>  frequentWordsWithMismatchesAndReverse(string text, int k,int d){
    string pattern,reverse_pattern;
    set<string> frequentPatterns, Neighborhood;
    int four_power_k = pow(4,k);
    vector<bool> close(four_power_k,false);
    vector<int> frequencyArray(four_power_k,0);

    for(size_t i = 0; i < text.size() - k + 1; ++i){
        Neighborhood = neighbors(text.substr(i, k), d);
        for(string pattern: Neighborhood)
            close[PatternToNumber(pattern)] = true;
    }
    for (int i = 0; i < four_power_k; ++i) {
        if(close[i]){
            pattern = NumberToPattern(i, k);
            reverse_pattern = makeReverse(pattern);
            frequencyArray[i] = ApproximatePatternCount(text, pattern, d) + ApproximatePatternCount(text, reverse_pattern, d);
        }
    }
    int maxCount = *max_element(frequencyArray.begin(),frequencyArray.end());
    answer<<"Max count: "<<maxCount<<endl;
    for (int i = 0; i < four_power_k; ++i) {
        if(frequencyArray[i] == maxCount){
            pattern = NumberToPattern(i, k);
            frequentPatterns.insert(pattern);
        }
    }
    return frequentPatterns;
}

void saveResults(set<string>& result,int d, int k, int minPos,int clumpLen){
    answer<<"For d = "<<d<<" k = "<<k<<" starting position: "<<minPos<<"/ len: "<<clumpLen<<endl
         <<"The best answer is:"<<endl;
    for(string s: result)
        answer<<s<<endl;
}

int main()
{
    string text = readData("data.txt");
    vector<int> minPos = minSkewCounter(text);
    int minPosStart = *minPos.begin();
    int length_clump = 500;
    string clump = text.substr(minPosStart,length_clump);
    int d = 2, k = 10;
    set<string> freq_words = frequentWordsWithMismatchesAndReverse(clump,k,d);
    saveResults(freq_words,d,k,minPosStart,length_clump);

    return 0;
}
