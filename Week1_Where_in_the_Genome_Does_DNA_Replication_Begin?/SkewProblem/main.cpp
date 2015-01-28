#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using namespace std;


static ofstream a("answer.txt");

void minSkewCounter(string input){
    int skew = 0;
    int minSkew = 0;
    vector<int> minPos;
    for (int i = 0; i < input.length(); ++i) {
        if (input[i] == 'G')
            skew++;
        if (input[i] == 'C')
            skew--;
        if(skew == minSkew){
            minPos.push_back(i+1);
        }
        if(skew < minSkew){
            minPos.clear();
            minPos.push_back(i+1);
            minSkew = skew;
        }

    }

    for(int i:minPos){
        cout<<i<<" ";
    }
}


int main()
{
    ifstream f("data.txt");
    string input;
    f>>input;
    minSkewCounter("GATACACTTCCCAGTAGGTACTG");
    return 0;
}

