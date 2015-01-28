#include <iostream>
#include <string.h>
#include <string>
#include <map>
#include <algorithm>
#include <vector>
#include <utility>
#include <fstream>
using namespace std;


static ofstream a("answer.txt");
bool myFunction(std::pair<string,int> fi, std::pair<string,int> se){
    return fi.second < se.second;
}


//Task 1.2.6
int countPattern(string input, string pattern){
    unsigned long long where =  input.find(pattern);
    unsigned long long  counter = 0;
    while (where!=std::string::npos){
        counter++;
        where = input.find(pattern,where + 1);
    }
    return counter;
}
//EndTask

//Task 1.2.9
void frequentWords(string input,int length){

    map<string,int> stringCounts;

    for(size_t i = 0 ; i< input.size() - length + 1; ++i)
        stringCounts[input.substr(i,length)]++;

    int max = (*max_element(stringCounts.begin(),stringCounts.end(),&myFunction)).second;

    //  for( map<string,int>::const_iterator it = stringCounts.begin(); it != stringCounts.end(); ++it)
    for(const auto it: stringCounts)
        if(it.second == max)
            a<<it.first<<" ";
}

//EndTask

void findPositions(string input, string genome){

    unsigned long long where =  input.find(genome);
    ofstream f("text.txt");

    while (where!=std::string::npos){
        f<<where<<" ";
        where = input.find(genome,where + 1);
    }

    f.close();

}

/*
 * //Defacto - slower version of task 1.2.9 - using faster stuff ...
int patternToNumber(string pattern){
    int sum = 0, val= 0;
    int len = pattern.length();
    for (int i = 0; i < len; ++i) {
        switch (pattern[len -1 -  i]) {
        case 'A': val = 0;  break;
        case 'C': val = 1; break;
        case 'G': val = 2; break;
        case 'T': val = 3; break;
        }
        sum+= val*pow(4,i);
    }
    return sum;
}

string numbersToPattern(int input, int len){
    int sum = input, digit = 0, powof4;
    map<int, char> letters = {{ 0,'A'},{ 1,'C'}, { 2,'G'}, { 3,'T'} };
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
        pattern[len - 1 - i] = letters[digit];
        sum-= digit * powof4;
    }
    return pattern;
}

vector<int> computingFrequencies(string input,int length){
    string pattern;
    vector<int> freqs(pow(4,length), 0);

    for(size_t i = 0 ; i< input.size() - length + 1; ++i){
        pattern = input.substr(i,length);
        int lok = patternToNumber(pattern);
        freqs[lok]++;
    }
    return freqs;
}


void fasterFrequentWords(string input,int length){
    vector<string> words;
    vector<int> freqs =  computingFrequencies(input, length);
    int maxCount =  *max_element(freqs.begin(),freqs.end());
    int pow4 = pow(4,length);
    for (int i = 0; i < pow4 ; ++i) {
        if(freqs[i] == maxCount)
            words.push_back(numbersToPattern(i,length));
    }
    for(auto word:words)
        a<<word<<" ";

}

*/



int main()
{

//    string input;
//    int length = 0;
//    ifstream f("data.txt");
//    f>>input>>length;
//    fasterFrequentWords(input,length);
//    f.close();


    //    findPositions(input,"ATGATCAAG");
    //cout<<endl;

    return 0;
}

