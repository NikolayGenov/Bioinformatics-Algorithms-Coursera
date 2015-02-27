#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>

using namespace std;
typedef char Symbol;
typedef string Symbols;
typedef vector<string> Patterns;
const Symbol SPECIAL_SYMBOLS[] = {'#', '$'};

class Trie {
    private:
    struct Node {
        Node() : value(numberNodes++) {
        }
        map<Symbol, Node*> children;

        bool has_edge_with_symbol(Symbol s) {
            return children[s] != nullptr;
        }

        Node* get_node_with_symbol(Symbol s) {
            return children[s];
        }

        Node* add_node(Symbol s) {
            children[s] = new Node();
            return children[s];
        }

        bool is_leaf() {
            return children.size() == 0;
        }

        ~Node() {
            for (auto pair : children)
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
        for (auto child : node->children) {
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

class SuffixTree {
    struct Node; // Forward declaration

    private:
    string text;
    size_t size;
    static Node* root;
    Node* lastNewNode;

    Node* activeNode;
    int activeEdge;
    int activeLength;

    int leafEnd;
    int* rootEnd;
    int* splitEnd;
    int remainingSuffixCount;

    struct Node {
        map<char, Node*> children;
        Node* suffixLink;

        int start;
        int* end;
        bool ownEnd;
        int suffixIndex;

        Node(int s, int* e, bool o) : children({}), suffixLink(root), start(s), end(e), ownEnd(o), suffixIndex(-1) {
        }
        ~Node() {
            for (auto it : children)
                delete it.second;

            if (ownEnd)
                delete end;
        }
        int get_edge_length() {
            return *end - start + 1;
        }
    };
    bool try_walk_down(Node& currNode) {
        int length = currNode.get_edge_length();
        if (activeLength >= length) {
            activeEdge += length;
            activeLength -= length;
            activeNode = &currNode;
            return true;
        }
        return false;
    }

    void extend_suffix_tree(int pos) {
        leafEnd = pos;
        remainingSuffixCount++;
        lastNewNode = nullptr;

        while (remainingSuffixCount > 0) {
            if (activeLength == 0)
                activeEdge = pos;

            if (activeNode->children[text[activeEdge]] == nullptr) {
                activeNode->children[text[activeEdge]] = new Node(pos, &leafEnd, false);

                if (lastNewNode != nullptr) {
                    lastNewNode->suffixLink = activeNode;
                    lastNewNode = nullptr;
                }
            } else {
                Node& next = *activeNode->children[text[activeEdge]];
                if (try_walk_down(next))
                    continue;

                if (text[next.start + activeLength] == text[pos]) {
                    if (lastNewNode != nullptr && activeNode != root) {
                        lastNewNode->suffixLink = activeNode;
                        lastNewNode = nullptr;
                    }

                    activeLength++;
                    break;
                }

                // New internal node
                splitEnd = new int(next.start + activeLength - 1);
                auto split = new Node(next.start, splitEnd, true);
                activeNode->children[text[activeEdge]] = split;

                // New leaf coming out of new internal node
                split->children[text[pos]] = new Node(pos, &leafEnd, false);
                next.start += activeLength;
                split->children[text[next.start]] = &next;

                if (lastNewNode != nullptr)
                    lastNewNode->suffixLink = split;

                lastNewNode = split;
            }

            remainingSuffixCount--;
            if (activeNode == root && activeLength > 0) {
                activeLength--;
                activeEdge = pos - remainingSuffixCount + 1;
            } else if (activeNode != root)
                activeNode = activeNode->suffixLink;
        }
    }

    void set_suffix_index_by_DFS(Node* n, int labelHeight) {
        if (n == nullptr)
            return;

        bool leaf = true;

        for (auto child : n->children) {
            leaf = false;
            set_suffix_index_by_DFS(child.second, labelHeight + child.second->get_edge_length());
        }
        if (leaf == true)
            n->suffixIndex = size - labelHeight;
    }

    void build_suffix_tree() {
        size = text.size();

        rootEnd = new int(-1);

        root = new Node(-1, rootEnd, true);
        root->suffixLink = nullptr;
        activeNode = root;

        for (size_t i = 0; i < size; i++)
            extend_suffix_tree(i);
    }

    void print(int i, int j, ofstream& os) {
        for (int k = i; k <= j; k++)
            os << text[k];
        os << endl;
    }

    void print_tree_by_DFS(Node* n, ofstream& os) {
        if (n == nullptr)
            return;

        if (n->start != -1)
            print(n->start, *(n->end), os);

        for (auto child : n->children)
            print_tree_by_DFS(child.second, os);
    }

    enum NodeType { LEAF = 0, INTERNAL = -1, X_NODE = -2, Y_NODE = -3, XY_NODE = -4 };

    void find_deepest_interal_node(Node& n, int labelHeight, int& maxHeight, int& substringStartIndex) {
        if (&n == nullptr)
            return;
        if (n.suffixIndex == INTERNAL)
            for (auto child : n.children)
                find_deepest_interal_node(*child.second, labelHeight + child.second->get_edge_length(), maxHeight, substringStartIndex);
        else if (n.suffixIndex > INTERNAL && maxHeight < labelHeight - n.get_edge_length()) {
            maxHeight = labelHeight - n.get_edge_length();
            substringStartIndex = n.suffixIndex;
        }
    }

    int find_deepest_common_interal_node(Node& n, int labelHeight, int& maxHeight, int& substringStartIndex, int sizeFirstString) {
        if (&n == nullptr)
            return LEAF;
        int nodeType;
        if (n.suffixIndex < LEAF)
            for (auto child : n.children) {
                nodeType = find_deepest_common_interal_node(*child.second,
                                                            labelHeight + child.second->get_edge_length(),
                                                            maxHeight,
                                                            substringStartIndex,
                                                            sizeFirstString);

                if (n.suffixIndex == INTERNAL)
                    n.suffixIndex = nodeType;
                else if ((n.suffixIndex == X_NODE && nodeType == Y_NODE) ||
                         (n.suffixIndex == Y_NODE && nodeType == X_NODE) || n.suffixIndex == XY_NODE) {
                    n.suffixIndex = XY_NODE;
                    if (maxHeight < labelHeight) {
                        maxHeight = labelHeight;
                        substringStartIndex = *n.end - labelHeight + 1;
                    }
                }
            }
        else if (n.suffixIndex > INTERNAL && n.suffixIndex < sizeFirstString)
            return X_NODE;
        else if (n.suffixIndex >= sizeFirstString)
            return Y_NODE;
        return n.suffixIndex;
    }

    void find_highest_non_shared_interal_node(Node& n, int labelHeight, int& minHeight, int& substringStartIndex, int sizeFirstString) {
        if (&n == nullptr)
            return;
        if (n.suffixIndex == INTERNAL) {
            for (auto child : n.children)
                if (child.first != SPECIAL_SYMBOLS[0] && child.first != SPECIAL_SYMBOLS[1])
                    find_highest_non_shared_interal_node(*child.second,
                                                         labelHeight + child.second->get_edge_length(),
                                                         minHeight,
                                                         substringStartIndex,
                                                         sizeFirstString);
        } else if (minHeight > labelHeight - n.get_edge_length()) {
            minHeight = labelHeight - n.get_edge_length() + 1;
            substringStartIndex = n.suffixIndex;
        }
    }

    public:
    SuffixTree(string str)
    : text(std::move(str)), activeNode(nullptr), activeEdge(-1), activeLength(0), leafEnd(-1), rootEnd(nullptr),
      splitEnd(nullptr), remainingSuffixCount(0) {
        build_suffix_tree();
        int labelHeight = 0;
        set_suffix_index_by_DFS(root, labelHeight);
    }

    void print_tree(ofstream& os) {
        print_tree_by_DFS(root, os);
    }

    string longest_repeated_substring() {
        int maxHeight = 0;
        int startIndex = 0;
        find_deepest_interal_node(*root, 0, maxHeight, startIndex);
        if (maxHeight == 0)
            return "";
        return text.substr(startIndex, maxHeight);
    }

    string longest_common_substring(int sizeFirstString) {
        int maxHeight = 0;
        int startIndex = 0;
        find_deepest_common_interal_node(*root, 0, maxHeight, startIndex, sizeFirstString);
        if (maxHeight == 0)
            return "";
        return text.substr(startIndex, maxHeight);
    }

    string shortest_non_shared_substring(int sizeFirstString) {
        int minHeight = sizeFirstString;
        int startIndex = 0;
        find_highest_non_shared_interal_node(*root, 0, minHeight, startIndex, sizeFirstString);
        if (minHeight == 0)
            return "";
        return text.substr(startIndex, minHeight);
    }


    ~SuffixTree() {
        delete root;
    }
};
SuffixTree::Node* SuffixTree::root = nullptr;


void read_file(ifstream& file, Symbols& s, int& sizeFirstString) {
    string str;
    file >> s;
    sizeFirstString = s.size();
    file >> str;
    s += SPECIAL_SYMBOLS[0] + str + SPECIAL_SYMBOLS[1];
    //    while (file >> str)
    //       p.push_back(str);
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
    int sizeFirstString = 0;
    read_file(file, syms, sizeFirstString);
    SuffixTree tree(syms);
    answer << tree.shortest_non_shared_substring(sizeFirstString) << endl;
    return 0;
}
