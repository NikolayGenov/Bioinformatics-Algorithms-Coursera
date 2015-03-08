#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <algorithm>
#include <cstring>

using namespace std;
typedef char Symbol;
typedef string Symbols;
typedef string Pattern;
typedef vector<string> Patterns;
const Symbol SPECIAL_SYMBOLS[] = {'$', '#'};

void sort_indices_by_lex_order(vector<int>& indices, const string& text) {
    size_t length = text.size();
    sort(indices.begin(), indices.end(), [length, text](int a, int b) {
        for (size_t i = 0; i < length; ++i)
            if (text[(a + i) % length] != text[(b + i) % length])
                return text[(a + i) % length] < text[(b + i) % length];
        return text[a] < text[b];
    });
}

struct OccurPair {
    int index;
    char ch;
    OccurPair(char c, int i) : index(i), ch(c) {
    }
    bool operator==(const OccurPair& other) const {
        return ch == other.ch && index == other.index;
    }
    bool operator<(const OccurPair& other) const {
        return ch < other.ch || (ch == other.ch && index < other.index);
    }
};

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
        if (n->start != -1)
            print(n->start, *(n->end), os);

        for (auto child : n->children)
            print_tree_by_DFS(child.second, os);
    }

    enum NodeType { LEAF = 0, INTERNAL = -1, X_NODE = -2, Y_NODE = -3, XY_NODE = -4 };

    void find_deepest_interal_node(Node& n, int labelHeight, int& maxHeight, int& substringStartIndex) {
        if (n.suffixIndex == INTERNAL)
            for (auto child : n.children)
                find_deepest_interal_node(*child.second, labelHeight + child.second->get_edge_length(), maxHeight, substringStartIndex);
        else if (n.suffixIndex > INTERNAL && maxHeight < labelHeight - n.get_edge_length()) {
            maxHeight = labelHeight - n.get_edge_length();
            substringStartIndex = n.suffixIndex;
        }
    }

    int find_deepest_common_interal_node(Node& n, int labelHeight, int& maxHeight, int& substringStartIndex, int sizeFirstString) {
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

// Taken and a bit modified from the version here:
// http://www.geeksforgeeks.org/suffix-array-set-2-a-nlognlogn-algorithm/
class SuffixArray {

    public:
    SuffixArray(const string& s) {
        buildSuffixArray(s);
    }

    vector<int> get_suffix_array() {
        return array;
    }

    vector<pair<int, int>> get_k_suffix_array(int k) {
        vector<pair<int, int>> k_array;
        for (size_t i = 0; i < array.size(); ++i)
            if (array[i] % k == 0)
                k_array.push_back(make_pair(i, array[i]));
        return k_array;
    }

    private:
    vector<int> array;

    struct suffix {
        int index;
        int rank[2];
    };

    bool cmp(struct suffix a, struct suffix b) const {
        return (a.rank[0] == b.rank[0]) ? (a.rank[1] < b.rank[1] ? true : false) : (a.rank[0] < b.rank[0] ? true : false);
    }

    void buildSuffixArray(const string& str) {
        size_t length = str.size();
        vector<suffix> suffixes(length);
        const char FIRST_LETTER = 'A';

        for (size_t i = 0; i < length; i++) {
            suffixes[i].index = i;
            suffixes[i].rank[0] = str[i] - FIRST_LETTER;
            suffixes[i].rank[1] = ((i + 1) < length) ? (str[i + 1] - FIRST_LETTER) : -1;
        }

        sort(suffixes.begin(), suffixes.end(), [this](suffix a, suffix b) { return cmp(a, b); });

        vector<int> ind(length);
        for (size_t k = 4; k < 2 * length; k = k * 2) {

            int rank = 0;
            int prev_rank = suffixes[0].rank[0];
            suffixes[0].rank[0] = rank;
            ind[suffixes[0].index] = 0;

            for (size_t i = 1; i < length; i++) {
                if (suffixes[i].rank[0] == prev_rank && suffixes[i].rank[1] == suffixes[i - 1].rank[1]) {
                    prev_rank = suffixes[i].rank[0];
                    suffixes[i].rank[0] = rank;
                } else {
                    prev_rank = suffixes[i].rank[0];
                    suffixes[i].rank[0] = ++rank;
                }
                ind[suffixes[i].index] = i;
            }

            for (size_t i = 0; i < length; i++) {
                size_t nextindex = suffixes[i].index + k / 2;
                suffixes[i].rank[1] = (nextindex < length) ? suffixes[ind[nextindex]].rank[0] : -1;
            }

            sort(suffixes.begin(), suffixes.end(), [this](suffix a, suffix b) { return cmp(a, b); });
        }

        array.resize(length);
        for (size_t i = 0; i < length; i++)
            array[i] = suffixes[i].index;
    }
};

class BurrowsWheelerTransformation {

    public:
    BurrowsWheelerTransformation(string s) : text(std::move(s)) {
        length = text.length();
        for (size_t i = 0; i < length; ++i)
            indices.push_back(i);

        sort_indices_by_lex_order(indices, text);

        for (size_t i = 0; i < length; ++i)
            transformed.push_back(text[(indices[i] + length - 1) % length]);
    }

    string get_transform() {
        return transformed;
    }

    static string get_inverse_tranformation(const string& s) {
        return inverse_transform(s);
    }

    int count_runs(int min_length) {
        int runs = 0;
        int counter = 1;
        bool new_word = false;
        for (size_t i = 0; i < length - 1; ++i)
            if (transformed[i] == transformed[i + 1]) {
                ++counter;
                if (counter >= min_length && new_word) {
                    ++runs;
                    new_word = false;
                }
            } else {
                new_word = true;
                counter = 1;
            }
        return runs;
    }

    private:
    size_t length;
    string text;
    vector<int> indices;
    string transformed = "";

    static string inverse_transform(const string& transformed) {
        string lexOrdered = transformed;
        string original;
        map<char, int> symbolsCount;
        vector<OccurPair> transformedPairs, lexOrderedPairs;

        sort(lexOrdered.begin(), lexOrdered.end());

        for (auto& elem : transformed)
            transformedPairs.emplace_back(elem, symbolsCount[elem]++);

        symbolsCount.clear();

        for (auto& elem : lexOrdered)
            lexOrderedPairs.emplace_back(elem, symbolsCount[elem]++);

        OccurPair p{SPECIAL_SYMBOLS[0], 0};
        for (auto& elem : transformed) {
            auto pos = find(transformedPairs.begin(), transformedPairs.end(), p) - transformedPairs.begin();
            p = lexOrderedPairs[pos];
            original.push_back(p.ch);
        }
        return original;
    }
};


vector<int> build_last_first_index(const string& first, const string& last) {
    map<char, int> symbolsCount;
    vector<OccurPair> firstPairs, lastPairs;

    for (auto& elem : last)
        lastPairs.emplace_back(elem, symbolsCount[elem]++);

    symbolsCount.clear();
    for (auto& elem : first)
        firstPairs.emplace_back(elem, symbolsCount[elem]++);

    vector<int> index;
    for (auto& lastPair : lastPairs) {
        auto pos = find(firstPairs.begin(), firstPairs.end(), lastPair) - firstPairs.begin();
        index.push_back(pos);
    }

    return index;
}

map<Symbol, int> build_first_occurrence_index(string last) {
    sort(last.begin(), last.end());
    map<Symbol, int> index;
    for (size_t i = 0; i < last.size(); ++i)
        if (index[last[i]] == 0)
            index[last[i]] = i;

    return index;
}

map<Symbol, vector<int>> build_checkout_symbols_index(const string& last, const vector<Symbol>& syms, int checkpointSize) {
    map<Symbol, vector<int>> index;
    map<Symbol, int> lastRow;

    for (Symbol sym : syms)
        lastRow[sym] = 0;


    for (size_t i = 0; i < last.size(); ++i)
        for (auto sym : syms) {
            int& value = lastRow[sym];

            if (i % checkpointSize == 0)
                index[sym].push_back(value);

            if (sym == last[i])
                value++;
        }
    return index;
}

map<Symbol, vector<int>> build_count_symbols_matrix(const string& last, const vector<Symbol>& syms) {
    map<Symbol, vector<int>> index;

    for (Symbol sym : syms)
        index[sym] = {0};

    for (size_t i = 0; i < last.size(); ++i)
        for (auto sym : syms) {
            int value = index[sym][i];
            if (sym == last[i])
                value++;
            index[sym].push_back(value);
        }
    return index;
}

vector<int> BWMatching(string lastColumn, Patterns patterns) {
    string firstColumn = lastColumn;
    sort(firstColumn.begin(), firstColumn.end());

    vector<int> indices = build_last_first_index(firstColumn, lastColumn);
    vector<int> numberOccurrences;

    for (auto pattern : patterns) {
        size_t top = 0;
        auto bottom = lastColumn.size() - 1;
        auto revIt = pattern.rbegin();
        auto revEndIt = pattern.rend();
        while (top <= bottom) {
            if (revIt != revEndIt) {
                Symbol sym = *revIt++;
                auto begin = lastColumn.begin();
                auto beginTop = begin + top;
                auto endBottom = begin + bottom + 1;
                auto first = find(beginTop, endBottom, sym);
                auto last = find_end(beginTop, endBottom, &sym, &sym + 1);
                if (first != endBottom) {
                    int topIndex = first - begin;
                    int bottomIndex = last - begin;
                    top = indices[topIndex];
                    bottom = indices[bottomIndex];
                } else {
                    numberOccurrences.push_back(0);
                    break;
                }
            } else {
                numberOccurrences.push_back(bottom - top + 1);
                break;
            }
        }
    }
    return numberOccurrences;
}


vector<int> better_BW_matching(string last, Patterns patterns) {

    auto firstOccurrence = build_first_occurrence_index(last);
    vector<Symbol> syms;
    for (auto pair : firstOccurrence)
        syms.push_back(pair.first);
    auto countSymbolsMatrix = build_count_symbols_matrix(last, syms);
    auto startColumn = last.begin();

    vector<int> numberOccurrences;
    for (auto pattern : patterns) {
        auto top = (size_t)0;
        auto bottom = last.size() - 1;
        auto rBegIt = pattern.rbegin();
        auto rEndIt = pattern.rend();
        while (top <= bottom) {
            if (rBegIt != rEndIt) {
                Symbol s = *rBegIt++;
                auto firstOc = startColumn + top;
                auto endOc = startColumn + bottom + 1;
                if (find(firstOc, endOc, s) != endOc) {
                    top = firstOccurrence[s] + countSymbolsMatrix[s][top];
                    bottom = firstOccurrence[s] + countSymbolsMatrix[s][bottom + 1] - 1;
                } else {
                    numberOccurrences.push_back(0);
                    break;
                }
            } else {
                numberOccurrences.push_back(bottom - top + 1);
                break;
            }
        }
    }
    return numberOccurrences;
}

int get_position_checkpoint(const string& last, Symbol s, map<Symbol, vector<int>>& index, int top, int checkpointSize) {
    int nearestCheckpoint = (top - 1) / checkpointSize;
    int start = nearestCheckpoint * checkpointSize;
    int count = 0;
    for (int i = start; i < top; ++i)
        if (last[i] == s)
            count++;
    return count + index[s][nearestCheckpoint];
}

void add_indices_with_partial_array(vector<int>& indices,
                                    const string last,
                                    map<Symbol, vector<int>>& index,
                                    const vector<pair<int, int>>& partialArray,
                                    map<Symbol, int>& first,
                                    int i,
                                    int checkpointSize) {

    int count = 0;
    int pos = i;
    auto it = find_if(partialArray.begin(), partialArray.end(), [pos](pair<int, int> k) { return k.first == pos; });
    while (it == partialArray.end()) {
        Symbol s = last[pos];
        pos = first[s] + get_position_checkpoint(last, s, index, pos, checkpointSize);
        count++;
        it = find_if(partialArray.begin(), partialArray.end(), [pos](pair<int, int> k) { return k.first == pos; });
    }
    int index_pos = it->second + count;
    indices.push_back(index_pos);
}

vector<int> efficient_index_BW_matching(const string& last,
                                        const Pattern& pattern,
                                        map<Symbol, int>& firstIndex,
                                        map<Symbol, vector<int>>& checkIndex,
                                        const vector<pair<int, int>>& suffArray,
                                        int checkpointSize,
                                        size_t& top,
                                        size_t& bottom) {

    auto startColumn = last.begin();
    vector<int> indices;

    auto rBegIt = pattern.rbegin();
    auto rEndIt = pattern.rend();

    while (top <= bottom) {
        if (rBegIt != rEndIt) {
            Symbol s = *rBegIt++;
            auto firstOc = startColumn + top;
            auto endOc = startColumn + bottom + 1;
            if (find(firstOc, endOc, s) != endOc) {
                top = firstIndex[s] + get_position_checkpoint(last, s, checkIndex, top, checkpointSize);
                bottom = firstIndex[s] + get_position_checkpoint(last, s, checkIndex, bottom + 1, checkpointSize) - 1;
            } else {
                break;
            }
        } else {
            for (size_t i = top; i < bottom + 1; ++i)
                add_indices_with_partial_array(indices, last, checkIndex, suffArray, firstIndex, i, checkpointSize);
            break;
        }
    }

    return indices;
}

vector<int> find_substring_pos_with_BWT(const string& str, const Patterns& patterns, int checkpointSize) {
    BurrowsWheelerTransformation trans(str);
    SuffixArray array(str);
    auto BWT = trans.get_transform();
    auto partialSuffixArray = array.get_k_suffix_array(checkpointSize);
    auto firstIndex = build_first_occurrence_index(BWT);
    vector<Symbol> syms;
    for (auto pair : firstIndex)
        syms.push_back(pair.first);

    auto checkIndex = build_checkout_symbols_index(BWT, syms, checkpointSize);

    vector<int> allIndices = {};
    for (auto pattern : patterns) {
        size_t top = 0;
        size_t bottom = BWT.size() - 1;
        auto indices =
        efficient_index_BW_matching(BWT, pattern, firstIndex, checkIndex, partialSuffixArray, checkpointSize, top, bottom);
        allIndices.insert(allIndices.end(), indices.begin(), indices.end());
    }

    sort(allIndices.begin(), allIndices.end());
    return allIndices;
}

vector<int> approximate_BWT_search(const string& last,
                                   const Pattern& pattern,
                                   map<Symbol, int>& firstIndex,
                                   map<Symbol, vector<int>>& checkIndex,
                                   const vector<pair<int, int>>& suffArray,
                                   int checkpointSize,
                                   int maxMismatches) {


    vector<int> indices;

    auto top = (size_t)0;
    auto bottom = last.size() - 1;
    auto rBegIt = pattern.rbegin();
    /*
        // Code for exact matching for the first
        // Need lots of changes to make it work - for example start
        // from the first minExacMatches + 1 symbols and find the rest
        // and not only the last one
        auto minExactMatches = pattern.size() / (maxMismatches + 1);
        auto remaining = pattern.size() - minExactMatches;
        string exactMatch;
        for (size_t i = 0; i < minExactMatches; ++i) {
            exactMatch.insert(exactMatch.begin(), *rBegIt++);
        }

        auto pos = efficient_index_BW_matching(last, exactMatch, firstIndex, checkIndex, suffArray, checkpointSize, top,
       bottom);
        if (pos.size() == 0)
            return indices;
    */
    vector<int> oldErrorCount(last.size(), 0);

    for (size_t step = 0; step < pattern.size(); ++step) {

        vector<int> errorCount(last.size(), -1);
        int bestTop = (int)last.size() + 1;
        int bestBottom = -1;
        for (size_t i = top; i < bottom + 1; ++i) {
            if (oldErrorCount[i] >= 0) {
                Symbol s = last[i];

                int localTop = firstIndex[s] + get_position_checkpoint(last, s, checkIndex, i, checkpointSize);
                int localBottom = firstIndex[s] + get_position_checkpoint(last, s, checkIndex, i + 1, checkpointSize) - 1;
                if (s != *rBegIt) {
                    errorCount[localTop] = oldErrorCount[i] + 1;
                    if (errorCount[localTop] > maxMismatches)
                        errorCount[localTop] = -1;
                } else {
                    errorCount[localTop] = oldErrorCount[i];
                }

                bestTop = min(localTop, bestTop);
                bestBottom = max(localBottom, bestBottom);
            }
        }
        if (bestTop == (int)last.size() + 1)
            return indices;
        oldErrorCount = errorCount;
        top = bestTop;
        bottom = bestBottom;
        rBegIt++;
    }

    for (size_t i = top; i < bottom + 1; ++i)
        if (oldErrorCount[i] >= 0)
            add_indices_with_partial_array(indices, last, checkIndex, suffArray, firstIndex, i, checkpointSize);

    return indices;
}

vector<int> multiple_approximate_BWT_search(const string& str, const Patterns& patterns, int checkpointSize, int maxMismatches) {
    BurrowsWheelerTransformation trans(str);
    SuffixArray array(str);
    auto BWT = trans.get_transform();
    auto partialSuffixArray = array.get_k_suffix_array(checkpointSize);
    auto firstIndex = build_first_occurrence_index(BWT);
    vector<Symbol> syms;
    for (auto pair : firstIndex)
        syms.push_back(pair.first);

    auto checkIndex = build_checkout_symbols_index(BWT, syms, checkpointSize);

    vector<int> allIndices = {};
    for (auto pattern : patterns) {
        auto indices =
        approximate_BWT_search(BWT, pattern, firstIndex, checkIndex, partialSuffixArray, checkpointSize, maxMismatches);
        allIndices.insert(allIndices.end(), indices.begin(), indices.end());
    }

    sort(allIndices.begin(), allIndices.end());
    return allIndices;
}

void read_file(ifstream& file, Symbols& s, Patterns& p) {
    string str;
    file >> s;
    s += "$";
    while (file >> str)
        p.push_back(str);
}
void print_vector(ofstream& os, const vector<int>& vec) {
    for (auto i : vec)
        os << i << " ";
    os << endl;
}

int main() {
    ofstream answer("answer.txt");
    ifstream file("data.txt");
    string str;
    Patterns p;
    read_file(file, str, p);
    int maxMismatches = stoi(p[p.size() - 1]);
    p.pop_back();

    const int CHECKPOINT_ARRAY_SIZE_CONST = 5;
    auto vec = multiple_approximate_BWT_search(str, p, CHECKPOINT_ARRAY_SIZE_CONST, maxMismatches);
    print_vector(answer, vec);
    return 0;
}
