#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <algorithm>
#include <set>
#include <cmath>

using namespace std;

static ofstream a("answer.txt");
//Using in both of them
int hammingDistance(const string& A,const string& B){
    int distance = 0;
    for (size_t i = 0; i < A.length(); ++i) {
        if(A[i] != B[i])
            distance++;
    }

    return distance;
}

int ApproximatePatternCount(string input,string pattern, int d){

    int hammingCount=0;
    string text;
    int len = pattern.length();
    for(size_t i = 0 ; i< input.size() - len + 1; ++i){
        text = input.substr(i,len);
        if(hammingDistance(pattern,text) <= d)
            hammingCount++;
    }

    return hammingCount;

}
static map<char, int> letters = {{'A',0},{'C',1}, {'G',2}, {'T',3} };
static map<int, char> lettersInverse = {{ 0,'A'},{ 1,'C'}, { 2,'G'}, { 3,'T'} };


bool myFunction(std::pair<string,int> fi, std::pair<string,int> se){
    return fi.second < se.second;
}


//Version 1.0

void generateAllK_words(map<string,int>& m, string seed, int pos,int k){
    if(pos == k)
        return;

    for(int i =0; i< 4;++i){
        seed[pos] = lettersInverse[i];
        m[seed] = 0;
        generateAllK_words(m,seed,pos+1,k);
    }
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

set<string> frequentWords_mismatch_and_reverse(string input, int k, int d){

    map<string,int> m;
    generateAllK_words(m,string(k,'A'),0,k);

    for(size_t i = 0 ; i< input.size() - k + 1; ++i)
        for(auto it:m)
            if(hammingDistance(it.first,input.substr(i,k)) <= d)
                m[it.first]++;


    map<string,int> pattern_only;
    typedef std::map<string, int>::iterator it_type;

    for(it_type it = m.begin(); it != m.end(); it++) {
        string rev = makeReverse(it->first);
        if(m[it->first] > 0 ){
            pattern_only[it->first] = m[it->first] + m[rev];
            m[rev] = 0;
        }
    }

    int max = (*max_element(pattern_only.begin(),pattern_only.end(),&myFunction)).second;
    set<string> result;
    for(const auto it:pattern_only)
        if(it.second == max){
            result.insert(it.first);
            result.insert(makeReverse(it.first));
        }
    return result;

}
/* ==================================================================================================================
   ==================================================================================================================
   ==================================================================================================================*/
//Version 2.0

// Generating the Neighborhood of a String

static vector<double> powmap;
static void calculatePows(int n){
    powmap.resize(n);
    for (int i = 0; i < n; ++i)
        powmap[i] = pow(4,i);
}

__int128  PatternToNumber(string pattern){
    if(pattern.empty()) return 0;
    int letc = letters[pattern.back()];
    pattern.pop_back();
    return 4 * PatternToNumber(pattern) + letc;
}

long long patternToNumber(const string& pattern){
    long long sum = 0;
    for (int i = pattern.length() - 1,j = 0 ; i >=0; --i,++j){
        sum+= letters[pattern[j]]*powmap[i];
    }
    return sum;
}


string NumberToPattern(int index,int  k){
    if (k == 1) return string(1,lettersInverse[index]);
    int prefixIndex = index /  4;
    int r = index % 4;
    return NumberToPattern(prefixIndex, k - 1) +  lettersInverse[r] ;

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
            close[patternToNumber(pattern)] = true;
    }
    for (int i = 0; i < four_power_k; ++i) {
        if(close[i]){
            pattern = numbersToPattern(i, k);
            reverse_pattern = makeReverse(pattern);
            frequencyArray[i] = ApproximatePatternCount(text, pattern, d) + ApproximatePatternCount(text, reverse_pattern, d);
        }
    }
    int maxCount = *max_element(frequencyArray.begin(),frequencyArray.end());
    cout<<"Max count: "<<maxCount<<endl;
    for (int i = 0; i < four_power_k; ++i) {
        if(frequencyArray[i] == maxCount){
            pattern = numbersToPattern(i, k);
            frequentPatterns.insert(pattern);
        }
    }
    return frequentPatterns;
}

void printResult(const set<string>& s){
    for(auto i:s)
        cout<<i<<" ";
    cout<<endl;
}
//Black magic
std::ostream&
operator<<( std::ostream& dest, __int128_t value )
{
    std::ostream::sentry s( dest );
    if ( s ) {
        __uint128_t tmp = value < 0 ? -value : value;
        char buffer[ 128 ];
        char* d = std::end( buffer );
        do
        {
            -- d;
            *d = "0123456789"[ tmp % 10 ];
            tmp /= 10;
        } while ( tmp != 0 );
        if ( value < 0 ) {
            -- d;
            *d = '-';
        }
        int len = std::end( buffer ) - d;
        if ( dest.rdbuf()->sputn( d, len ) != len ) {
            dest.setstate( std::ios_base::badbit );
        }
    }
    return dest;
}

void test(){
    ifstream f("data.txt");
    string text("AAATCAACGAGTCGACTAAGGGATATGTGAAGAGAGCAGGTATGTCGGACGGTA");
    //    f>>text;
    //    int number = 5895;
    //    int k = 9;
    __int128 answer = PatternToNumber(text);
    a   <<answer;
    f.close();
}

int main()
{

    calculatePows(200);

    //	ifstream f("ecoli.txt");
    //	string input;
    //	int d = 1,k  =9 ;
    //	f>>input;
    //	int position_to_start=  3923820;
    //	int length = 400;
    //	string text = input.substr(position_to_start,length);


    //	printResult(frequentWordsWithMismatchesAndReverse(text,k,d));
//cout<<ApproximatePatternCount("CGTGACAGTGTATGGGCATCTTT","TGT",1);
    // cout<<hammingDistance("CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA", "CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG");
   int count = 0;
    for(auto b: neighbors("CCCC",3))
       count++;
    cout<<count;

    return 0;
}

