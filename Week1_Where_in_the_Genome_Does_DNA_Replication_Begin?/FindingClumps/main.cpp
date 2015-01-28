#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <fstream>
#include <set>
#include <cmath>
#include <cstring>
#include <array>

using namespace std;
//For task 1.4.4

void L_t_clumps (char* input,int size,int k, int L, int t); //version 1.0

void ClumpFinding(string input,int k, int L, int t); //version 2.0

void BetterClumpFinding (string input,int k, int L, int t); //version 2.0
//bool myFunction(std::pair<char[10],int> fi, std::pair<char[10],int> se){
//    return fi.second < se.second;
//}

struct cmp_str
{
    bool operator()(char const *a, char const *b)
    {
        return std::strcmp(a, b) < 0;
    }
};

static vector<double> powmap;
static void calculatePows(int n){
    powmap.resize(n);
    for (int i = 0; i < n; ++i)
        powmap[i] = pow(4,i);

}

void printResult(set<string>& res){
    ofstream f("answer");
    for(auto s: res)
        f<<s<<" ";

    f.close();
}

void readAndExec(const char * filename){
    ifstream f(filename);
    string input;
    // int k,L,t;
    // f>>input>>k>>L>>t;
    f>>input;

    ClumpFinding("GCACAAGGCCGACAATAGGACGTAGCCTTGAAGACGACGTAGCGTGGTCGCATAAGTACAGTAGATAGTACCTCCCCCGCGCATCCTATTATTAAGTTAATT",4,30,3);
    f.close();
}
/* //Problem
set<string> k_mersFoundT_times(char* input,int start, int end,int length, int times){
    map<char* ,int,cmp_str> stringCounts;
    int size = end -  start;

    for(size_t i = 0 ; i< size - length + 1; ++i){
     //  char* buffer = new  char[10]; //problem here - go back to str ?

        strncpy(buffer,input+start+i,length);
        stringCounts[buffer]++;

}
    set<string> words;
    for(const auto it: stringCounts)
        if(it.second >= times){
            cout<<it.first<<endl;
            words.insert(it.first);
    }
            return words;
}

void L_t_clumps (char* input,int size,int k, int L, int t){

    set<string> result, words;
    cout<<size<<endl;
    for(size_t i = 0 ; i< size - L + 1; ++i){
        words =  k_mersFoundT_times(input,i,i+L,k,t);

        result.insert(words.begin(), words.end());
        if (i %1000 == 0&& i >0 )
            printf("%d\n",i);
    }
    printResult(result);
}
*/
//Version 2.0

static map<char, int> letters = {{'A',0},{'C',1}, {'G',2}, {'T',3} };
static map<int, char> lettersInverse = {{ 0,'A'},{ 1,'C'}, { 2,'G'}, { 3,'T'} };

int patternToNumber(const string& pattern){
    int sum = 0;
    for (int i = pattern.length() - 1,j = 0 ; i >=0; --i,++j){
        sum+= letters[pattern[j]]*powmap[i];
    }
    return sum;
}

string numbersToPattern(int input, int len){
    int sum = input, digit = 0, powof4;

    string pattern(len,' ');
    for (int i = len-1; i >=0; --i) {
        powof4 = pow(4,i);
        digit =  0;
        for(int d = 0; d <= 3; ++d){
            if(d * powof4 > sum)
                break;
            else
                digit = d;
        }
        pattern[len - 1 - i] = lettersInverse[digit];
        sum-= digit * powof4;
    }
    return pattern;
}

void computingFrequencies(vector<int>& freqs, const string& input,int length){
    int lok = 0;
    for(size_t i = 0 ; i< input.size() - length + 1; ++i){
        lok = patternToNumber(input.substr(i,length));
        freqs[lok]++;
    }
}

void ClumpFinding(string input,int k, int L, int t){

    set<string> freqPatterns;
    int po4 = pow(4,k);
    static vector<bool> clump(po4,false);
    vector<int> freqs(po4, 0);
    cout<<po4<<endl;
    for(size_t i = 0 ; i< input.size() - L + 1; ++i){

        computingFrequencies(freqs,input.substr(i,L),k);
        for (int j = 0; j < po4; ++j) {
            if(freqs[j] >= t)
                clump[j] = true;
        }
        if (i % 1000 == 0)
            printf("%d\n",i);
        std::fill(freqs.begin(), freqs.end(), 0);
    }

    for (int i = 0; i < po4; ++i)
        if( clump[i])
            freqPatterns.emplace(numbersToPattern(i,k));

    printResult(freqPatterns);
}

//Version 3.0 - Fastest possible

void BetterClumpFinding(string input,int k, int L, int t){

    int po4 = pow(4,k),j;
    vector<bool> clump(po4,false);
    vector<int> freqs(po4, 0);

    string text = input.substr(0,L);
    computingFrequencies(freqs,text,k);
    for (int j = 0; j < po4; ++j)
        if(freqs[j] >= t)
            clump[j] = true;

    for(size_t i = 1 ; i< input.size() - L + 1; ++i){
        j = patternToNumber(input.substr(i-1,k));
        --freqs[j];

        j = patternToNumber(input.substr(i + L -k,k));
        ++freqs[j];

        if(freqs[j] >= t)
            clump[j] = true;
/* //Speed test
        if (i % 1000 == 0)
            printf("%d\n",i);
*/
    }

    set<string> freqPatterns;
    for (int i = 0; i < po4; ++i)
        if( clump[i])
            freqPatterns.emplace(numbersToPattern(i,k));

   // printResult(freqPatterns);
    cout<<"The answer for E-coli is: "<< freqPatterns.size()<<endl; //The answer
}


int main()
{

    calculatePows(100);
   readAndExec("data.txt");
   //   string str = "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA";
   // BetterClumpFinding("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA",5,50,4);

   return 0;
}

