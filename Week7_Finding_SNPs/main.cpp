#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>

using namespace std;
typedef char Symbol;
typedef vector<string> Patterns;

class Trie {
    private:
    struct Node {
        Node() : value(numberNodes++) {
        }
        map<Symbol, Node*> childrens;

        bool has_edge_with_symbol(Symbol s) {
            return childrens[s] != nullptr;
        }

        Node* get_node_with_symbol(Symbol s) {
            return childrens[s];
        }

        Node* add_node(Symbol s) {
            childrens[s] = new Node();
            return childrens[s];
        }

        ~Node() {
            for (auto pair : childrens)
                delete pair.second;
        }

        int value;
    };

    public:
    Trie(Patterns patterns) {
        root = new Node();
        for (auto str : patterns) {
            Node* currentNode = root;
            for (auto sym : str)
                if (currentNode->has_edge_with_symbol(sym))
                    currentNode = currentNode->get_node_with_symbol(sym);
                else
                    currentNode = currentNode->add_node(sym);
        }
    }

    void print_trie(ostream& os) {
        rec_print(root, os);
    }

    ~Trie() {
        delete root;
    }

    private:
    Node* root;
    static int numberNodes;

    void rec_print(Node* start, ostream& os) {
        Node* node = start;
        for (auto child : node->childrens) {
            os << node->value << "->" << child.second->value << ":" << child.first << endl;
            rec_print(child.second, os);
        }
    }
};
int Trie::numberNodes = 0;

void read_file(ifstream& file, Patterns& p) {
    string str;
    while (file >> str)
        p.push_back(str);
}

int main() {
    ofstream answer("answer.txt");
    ifstream file("data.txt");
    Patterns p;
    read_file(file, p);
    Trie trie(p);
    trie.print_trie(answer);

    return 0;
}
