#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>

using namespace std;
typedef char Symbol;
typedef string Symbols;
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

        bool is_leaf() {
            return childrens.size() == 0;
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

    vector<int> matching(const Symbols& text) {
        size_t pos = 0;
        vector<int> positions;
        while (pos < text.size()) {
            if (prefix_matching(text.substr(pos)))
                positions.push_back(pos);
            ++pos;
        }
        return positions;
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

    bool prefix_matching(const Symbols& text) {
        int count = 0;
        Symbol sym = text[count];
        Node* v = root;
        while (true) {
            if (v->is_leaf()) {
                return true;
            } else if (v->has_edge_with_symbol(sym)) {
                ++count;
                v = v->get_node_with_symbol(sym);
                sym = text[count];
            } else
                return false;
        }
    }
};
int Trie::numberNodes = 0;

void read_file(ifstream& file, Symbols& s, Patterns& p) {
    string str;
    file >> s;
    while (file >> str)
        p.push_back(str);
}
void print_vector(ofstream& os, vector<int>& vec) {
    for (auto i : vec)
        os << i << " ";
    os << endl;
}

int main() {
    ofstream answer("answer.txt");
    ifstream file("data.txt");
    Symbols syms;
    Patterns p;
    read_file(file, syms, p);
    Trie trie(p);
    auto vec = trie.matching(syms);
    print_vector(answer, vec);
    return 0;
}
