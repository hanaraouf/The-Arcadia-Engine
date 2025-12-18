// ArcadiaEngine.cpp - STUDENT TEMPLATE
// TODO: Implement all the functions below according to the assignment requirements

#include "ArcadiaEngine.h"
#include <algorithm>
#include <queue>
#include <numeric>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <set>

using namespace std;

// =========================================================
// PART A: DATA STRUCTURES (Concrete Implementations)
// =========================================================

// --- 1. PlayerTable (Double Hashing) ---

class ConcretePlayerTable : public PlayerTable {
private:
    // TODO: Define your data structures here
    // Hint: You'll need a hash table with double hashing collision resolution

    static const int TABLE_SIZE = 101;
    vector<pair<int, string>> table;
    vector<bool> occupied;
    
    int hashFunction1(int key) {
        // Mid-Square method
        long long squared = (long long)key * key;
        int middleDigits = (squared / 100) % 10000;
        return middleDigits % TABLE_SIZE;
    }
    
    // for double hashing
    int hashFunction2(int key) {
        return 1 + (key % (TABLE_SIZE - 1));
    }

public:
    ConcretePlayerTable() {
        // TODO: Initialize your hash table
        table.resize(TABLE_SIZE, {-1, ""});
        occupied.resize(TABLE_SIZE, false);
    }

    void insert(int playerID, string name) override {
        // TODO: Implement double hashing insert
        // Remember to handle collisions using h1(key) + i * h2(key)
        int h1 = hashFunction1(playerID);
        int h2 = hashFunction2(playerID);
        
        for (int i = 0; i < TABLE_SIZE; i++) {
            int index = (h1 + i * h2) % TABLE_SIZE;
            if (!occupied[index]) {
                table[index] = {playerID, name};
                occupied[index] = true;
                return;
            }
        }
        cout << "Error: Table is full" << endl;
    }

    string search(int playerID) override {
        // TODO: Implement double hashing search
        // Return "" if player not found
        int h1 = hashFunction1(playerID);
        int h2 = hashFunction2(playerID);
        
        for (int i = 0; i < TABLE_SIZE; i++) {
            int index = (h1 + i * h2) % TABLE_SIZE;
            if (!occupied[index]){
                break;
            }
            if (table[index].first == playerID) {
                return table[index].second;
            }
        }
        return "";
    }
};

// --- 2. Leaderboard (Skip List) ---

class ConcreteLeaderboard : public Leaderboard {
private:
    // TODO: Define your skip list node structure and necessary variables
    // Hint: You'll need nodes with multiple forward pointers
    struct SkipNode {
        int playerID;
        int score;
        SkipNode** forward;
        int level;
        
        SkipNode(int id, int sc, int lvl) : playerID(id), score(sc), level(lvl) {
            forward = new SkipNode*[lvl + 1];
            for (int i = 0; i <= lvl; i++) {
                forward[i] = nullptr;
            }
        }
        
        ~SkipNode() {
            delete[] forward;
        }
    };
    
    static const int MAX_LEVEL = 16;
    SkipNode* header;
    int currentLevel;
    
    int randomLevel() {
        int level = 0;
        while (level < MAX_LEVEL && (rand() % 2 == 0)) {
            level++;
        }
        return level;
    }
    
    bool shouldComeAfter(int score1, int id1, int score2, int id2) {
        if (score1 < score2) {
            return true;
        }
        if (score1 > score2) {
            return false;
        }
        return id1 > id2;
    }

public:
    ConcreteLeaderboard() {
        // TODO: Initialize your skip list
        header = new SkipNode(-1, INT_MAX, MAX_LEVEL);
        currentLevel = 0;
    }

    void addScore(int playerID, int score) override {
        // TODO: Implement skip list insertion
        // Remember to maintain descending order by score
        SkipNode* update[MAX_LEVEL + 1];
        for (int i = 0; i <= MAX_LEVEL; i++) {
            update[i] = nullptr;
        }
        
        SkipNode* current = header;
        
        for (int i = currentLevel; i >= 0; i--) {
            while (current->forward[i] != nullptr &&
                   !shouldComeAfter(current->forward[i]->score, current->forward[i]->playerID, score, playerID)) {
                current = current->forward[i];
            }
            update[i] = current;
        }
        
        int newLevel = randomLevel();
        if (newLevel > currentLevel) {
            for (int i = currentLevel + 1; i <= newLevel; i++) {
                update[i] = header;
            }
            currentLevel = newLevel;
        }
        
        SkipNode* newNode = new SkipNode(playerID, score, newLevel);
        for (int i = 0; i <= newLevel; i++) {
            newNode->forward[i] = update[i]->forward[i];
            update[i]->forward[i] = newNode;
        }
    }

    void removePlayer(int playerID) override {
        // TODO: Implement skip list deletion
        SkipNode* update[MAX_LEVEL + 1];
        for (int i = 0; i <= MAX_LEVEL; i++) {
            update[i] = nullptr;
        }
        
        SkipNode* current = header;
        
        for (int i = currentLevel; i >= 0; i--) {
            while (current->forward[i] != nullptr &&
                   current->forward[i]->playerID != playerID) {
                current = current->forward[i];
            }
            update[i] = current;
        }
        
        current = current->forward[0];
        
        if (current != nullptr && current->playerID == playerID) {
            for (int i = 0; i <= currentLevel; i++) {
                if (update[i]->forward[i] != current) break;
                update[i]->forward[i] = current->forward[i];
            }
            delete current;
            
            while (currentLevel > 0 && header->forward[currentLevel] == nullptr) {
                currentLevel--;
            }
        }
    }

    vector<int> getTopN(int n) override {
        // TODO: Return top N player IDs in descending score order
        vector<int> result;
        SkipNode* current = header->forward[0];
        
        while (current != nullptr && result.size() < (size_t)n) {
            result.push_back(current->playerID);
            current = current->forward[0];
        }
        
        return result;
    }
};

// --- 3. AuctionTree (Red-Black Tree) ---

class ConcreteAuctionTree : public AuctionTree {
private:
    // TODO: Define your Red-Black Tree node structure
    // Hint: Each node needs: id, price, color, left, right, parent pointers
    
    static const int RED = 0;
    static const int BLACK = 1;
    
    // Node Structure
    struct RBNode {
        int itemID;
        int price;
        int color;
        RBNode* left;
        RBNode* right;
        RBNode* parent;
        
        RBNode(int id, int p) : itemID(id), price(p), color(RED),
                                 left(nullptr), right(nullptr), parent(nullptr) {}
    };
    
    RBNode* root;
    
// Left Rotation: Rotate node x down to the left, make its right child y parent
    void rotateLeft(RBNode* x) {
        RBNode* y = x->right;
        x->right = y->left;
        
        // Update parent pointer of y's left child
        if (y->left != nullptr) {
            y->left->parent = x;
        }
        
        y->parent = x->parent;
        
        if (x->parent == nullptr) {
            root = y;
        } else if (x == x->parent->left) {
            x->parent->left = y;
        } else {
            x->parent->right = y;
        }
        
        y->left = x;
        x->parent = y;
    }
    
// Right Rotation: Rotate node x down to the right, make its left child y parent
    void rotateRight(RBNode* x) {
        RBNode* y = x->left;
        x->left = y->right;
        
        if (y->right != nullptr) {
            y->right->parent = x;
        }
        
        y->parent = x->parent;
        
        // Update x's parent to point to y instead of x
        if (x->parent == nullptr) {
            root = y;
        } else if (x == x->parent->right) {
            x->parent->right = y;
        } else {
            x->parent->left = y;
        }
        
        y->right = x;
        x->parent = y;
    }
    
    
//     Fix violations after insertion
//     Root must be BLACK
//     RED node must have BLACK children

//     Cases :
//     case 1. Uncle is RED: Recolor parent, uncle, grandparent
//     case 2. Uncle is BLACK, node is right child: Rotate to convert to case 3
//     case 3. Uncle is BLACK, node is left child: Rotate and recolor
    
    void fixInsert(RBNode* x) {
        while (x->parent != nullptr && x->parent->color == RED) {
            
// Case: P is LEFT child of GP
            if (x->parent == x->parent->parent->left) {
                // y is uncle
                RBNode* y = x->parent->parent->right;
                
                // Case 1
                if (y != nullptr && y->color == RED) {
                    x->parent->color = BLACK;           // Recolor P to BLACK
                    y->color = BLACK;                    // Recolor U to BLACK
                    x->parent->parent->color = RED;     // Recolor GP to RED
                    x = x->parent->parent;              // Move x up to GP
                }
                else {
                    // Case 2
                    if (x == x->parent->right) {
                        x = x->parent;
                        rotateLeft(x);
                    }
                    // Case 3
                    x->parent->color = BLACK;           // Recolor P to BLACK
                    x->parent->parent->color = RED;     // Recolor GP to RED
                    rotateRight(x->parent->parent);     // Rotate right at GP
                }
            }
            
            
// Case: P is RIGHT child of GP
            else {
                // y is uncle
                RBNode* y = x->parent->parent->left;
                
                // Case 1
                if (y != nullptr && y->color == RED) {
                    x->parent->color = BLACK;           // Recolor P to BLACK
                    y->color = BLACK;                    // Recolor U to BLACK
                    x->parent->parent->color = RED;     // Recolor GP to RED
                    x = x->parent->parent;              // Move x up to GP
                }
                else {
                    // Case 2
                    if (x == x->parent->left) {
                        x = x->parent;
                        rotateRight(x);
                    }
                    // Case 3
                    x->parent->color = BLACK;           // Recolor P to BLACK
                    x->parent->parent->color = RED;     // Recolor GP to RED
                    rotateLeft(x->parent->parent);      // Rotate left at GP
                }
            }
        }
        root->color = BLACK;
    }
    

// Replace subtree at root u with subtree at root v
    void moveSubtree(RBNode* u, RBNode* v) {
        if (u->parent == nullptr) {
            root = v;
        } else if (u == u->parent->left) {
            u->parent->left = v;
        } else {
            u->parent->right = v;
        }
        if (v != nullptr) {
            v->parent = u->parent;
        }
    }
    
//minimum
    RBNode* minimum(RBNode* node) {
        while (node->left != nullptr) {
            node = node->left;
        }
        return node;
    }
    

    void fixDelete(RBNode* x, RBNode* xParent) {
        // Continue while x is not root and is BLACK (double-black)
        while (x != root && (x == nullptr || x->color == BLACK)) {
            
            // Case x is LEFT child
            if (x == xParent->left) {
                RBNode* w = xParent->right;
                
                if (w->color == RED) {
                    w->color = BLACK;               // Recolor w to BLACK
                    xParent->color = RED;           // Recolor P to RED
                    rotateLeft(xParent);            // Rotate left at P
                    w = xParent->right;             // Update sibling
                }
                
                if ((w->left == nullptr || w->left->color == BLACK) &&
                    (w->right == nullptr || w->right->color == BLACK)) {
                    w->color = RED;                 // Recolor w to RED
                    x = xParent;
                    xParent = x->parent;
                }
                else {
                    if (w->right == nullptr || w->right->color == BLACK) {
                        if (w->left != nullptr) w->left->color = BLACK;
                        w->color = RED;
                        rotateRight(w);
                        w = xParent->right;
                    }
                    w->color = xParent->color;                          // w takes P color
                    xParent->color = BLACK;                             // P becomes BLACK
                    if (w->right != nullptr) w->right->color = BLACK;
                    rotateLeft(xParent);                                // Rotate left at P
                    x = root;                                           //stop loop
                }
            }
            // Case: x is RIGHT child
            else {
                RBNode* w = xParent->left;   // w is sibling
                
                if (w->color == RED) {
                    w->color = BLACK;               // Recolor w to BLACK
                    xParent->color = RED;           // Recolor P to RED
                    rotateRight(xParent);           // Rotate right at P
                    w = xParent->left;              // Update sibling
                }
                
                if ((w->right == nullptr || w->right->color == BLACK) &&
                    (w->left == nullptr || w->left->color == BLACK)) {
                    w->color = RED;
                    x = xParent;
                    xParent = x->parent;
                }
                else {
                    if (w->left == nullptr || w->left->color == BLACK) {
                        if (w->right != nullptr) w->right->color = BLACK;
                        w->color = RED;
                        rotateLeft(w);
                        w = xParent->left;
                    }
                    w->color = xParent->color;                      // w takes P color
                    xParent->color = BLACK;                         // P becomes BLACK
                    if (w->left != nullptr) w->left->color = BLACK;
                    rotateRight(xParent);                           // Rotate right at P
                    x = root;                                       // stop loop
                }
            }
        }
        if (x != nullptr) x->color = BLACK;
    }
    
//Find node by itemID
    RBNode* findNode(int itemID, RBNode* node) {
        if (node == nullptr) return nullptr;
        if (node->itemID == itemID) return node;
        
        RBNode* found = findNode(itemID, node->left);
        if (found != nullptr) return found;
        
        return findNode(itemID, node->right);
    }

public:
    ConcreteAuctionTree() {
        // TODO: Initialize your Red-Black Tree
        root = nullptr;  // Tree starts empty
    }


    void insertItem(int itemID, int price) override {
        // TODO: Implement Red-Black Tree insertion
        // Remember to maintain RB-Tree properties with rotations and recoloring
        
        // Create new RED node
        RBNode* z = new RBNode(itemID, price);
        RBNode* y = nullptr;
        RBNode* x = root;
        
        //  find position
        while (x != nullptr) {
            y = x;
            
            if (z->price < x->price || (z->price == x->price && z->itemID < x->itemID)) {
                x = x->left;
            } else {
                x = x->right;
            }
        }
        
        z->parent = y;
        
        // Tree was empty, z becomes root
        if (y == nullptr) {
            root = z;
        }
        // z goes to left
        else if (z->price < y->price || (z->price == y->price && z->itemID < y->itemID)) {
            y->left = z;
        }
        // z goes to right
        else {
            y->right = z;
        }
        
        fixInsert(z);
    }


    void deleteItem(int itemID) override {
        // TODO: Implement Red-Black Tree deletion
        // This is complex - handle all cases carefully
        
        // Step 1: Find the node to delete
        RBNode* z = findNode(itemID, root);
        if (z == nullptr) return;
        
        // Track the node that will be removed and its color
        RBNode* y = z;
        RBNode* x;
        RBNode* xParent;
        int yOriginalColor = y->color;
        
        // Case 1: z has no left child
        if (z->left == nullptr) {
            x = z->right;
            xParent = z->parent;
            moveSubtree(z, z->right);  // Replace z with right child
        }
        // Case 2: z has no right child
        else if (z->right == nullptr) {
            x = z->left;
            xParent = z->parent;
            moveSubtree(z, z->left);   // Replace z with left child
        }
        // Case 3: z has two children
        else {
            // Find successor
            y = minimum(z->right);
            yOriginalColor = y->color;
            x = y->right;  // y's right child will replace y
            
            // Check if y is direct child of z
            if (y->parent == z) {
                xParent = y;  // x's parent will be y
            } else {
                xParent = y->parent;
                moveSubtree(y, y->right);  // Replace y with its right child
                y->right = z->right;      // Give y's right subtree to z
                y->right->parent = y;
            }
            
            // Replace z with y
            moveSubtree(z, y);
            y->left = z->left;
            y->left->parent = y;
            y->color = z->color;  // y takes z's color
        }
        
        delete z;
        
        if (yOriginalColor == BLACK) {
            fixDelete(x, xParent);
        }
    }
};

// =========================================================
// PART B: INVENTORY SYSTEM (Dynamic Programming)
// =========================================================

int InventorySystem::optimizeLootSplit(int n, vector<int>& coins) {
    // TODO: Implement partition problem using DP
    // Goal: Minimize |sum(subset1) - sum(subset2)|
    // Hint: Use subset sum DP to find closest sum to total/2
    return 0;
}

int InventorySystem::maximizeCarryValue(int capacity, vector<pair<int, int>>& items) {
    // TODO: Implement 0/1 Knapsack using DP
    // items = {weight, value} pairs
    // Return maximum value achievable within capacity
    return 0;
}

long long InventorySystem::countStringPossibilities(string s) {
    // TODO: Implement string decoding DP
    // Rules: "uu" can be decoded as "w" or "uu"
    //        "nn" can be decoded as "m" or "nn"
    // Count total possible decodings
    return 0;
}

// =========================================================
// PART C: WORLD NAVIGATOR (Graphs)
// =========================================================

bool WorldNavigator::pathExists(int n, vector<vector<int>>& edges, int source, int dest) {
    if (source == dest) return true;
    
    // Build adjacency list for undirected graph
    vector<vector<int>> adj(n);
    for (const auto& edge : edges) {
        int u = edge[0], v = edge[1];
        adj[u].push_back(v);
        adj[v].push_back(u);  // Bidirectional
    }
    
    // BFS for undirected graph
    vector<bool> visited(n, false);
    queue<int> q;
    
    visited[source] = true;
    q.push(source);
    
    while (!q.empty()) {
        int curr = q.front();
        q.pop();
        
        if (curr == dest) return true;
        
        for (int neighbor : adj[curr]) {
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                q.push(neighbor);
            }
        }
    }
    
    return false;
}

// Union-Find (Disjoint Set Union) for Kruskal's MST
struct DSU {
    vector<int> parent, rank;
    
    DSU(int n) {
        parent.resize(n);
        rank.resize(n, 0);
        for (int i = 0; i < n; i++) parent[i] = i;
    }
    
    int find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]);  // Path compression
        }
        return parent[x];
    }
    
    bool unite(int x, int y) {
        int rx = find(x);
        int ry = find(y);
        
        if (rx == ry) return false;
        
        // Union by rank
        if (rank[rx] < rank[ry]) {
            parent[rx] = ry;
        } else if (rank[rx] > rank[ry]) {
            parent[ry] = rx;
        } else {
            parent[ry] = rx;
            rank[rx]++;
        }
        return true;
    }
};

long long WorldNavigator::minBribeCost(int n, int m, long long goldRate, long long silverRate,
                                       vector<vector<int>>& roadData) {
    // Build edge list with computed costs (undirected graph)
    vector<vector<long long>> edges;
    for (const auto& road : roadData) {
        int u = road[0], v = road[1];
        long long gold = road[2], silver = road[3];
        long long cost = gold * goldRate + silver * silverRate;
        edges.push_back({cost, (long long)u, (long long)v});
    }
    
    // Kruskal's MST algorithm for undirected graph
    sort(edges.begin(), edges.end());
    
    DSU dsu(n);
    long long total = 0;
    int edgesUsed = 0;
    
    for (const auto& e : edges) {
        long long cost = e[0];
        int u = e[1], v = e[2];
        
        if (dsu.unite(u, v)) {
            total += cost;
            edgesUsed++;
            if (edgesUsed == n - 1) break;  // MST complete
        }
    }
    
    // Check if all cities can be connected
    return (edgesUsed == n - 1) ? total : -1;
}

string WorldNavigator::sumMinDistancesBinary(int n, vector<vector<int>>& roads) {
    const long long INF = 1e18;
    vector<vector<long long>> dist(n, vector<long long>(n, INF));
    
    // Initialize distances
    for (int i = 0; i < n; i++) dist[i][i] = 0;
    
    // Build directed graph (teleporter network)
    for (const auto& road : roads) {
        int u = road[0], v = road[1];
        long long w = road[2];
        dist[u][v] = min(dist[u][v], w);  // Directed edge u->v
    }
    
    // Floyd-Warshall for all-pairs shortest paths in directed graph
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            if (dist[i][k] == INF) continue;  // Optimization
            for (int j = 0; j < n; j++) {
                if (dist[k][j] == INF) continue;  // Optimization
                dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
            }
        }
    }
    
    // Sum distances for all pairs i < j (directed distance i→j)
    long long sum = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (dist[i][j] < INF) {
                sum += dist[i][j];
            }
        }
    }
    
    // Convert to binary string
    if (sum == 0) return "0";
    
    string binary = "";
    while (sum > 0) {
        binary = to_string(sum % 2) + binary;
        sum /= 2;
    }
    
    return binary;
}
// =========================================================
// PART D: SERVER KERNEL (Greedy)
// =========================================================

int ServerKernel::minIntervals(vector<char>& tasks, int n) {
    // TODO: Implement task scheduler with cooling time
    // Same task must wait 'n' intervals before running again
    // Return minimum total intervals needed (including idle time)
    // Hint: Use greedy approach with frequency counting
    
    // Base cases 
    if(tasks.empty()) return 0;
    if(tasks.size() == 1) return 1;
    if(n == 0) return tasks.size(); // No cooling time

    vector<char> schedule; // Vector to store the task or idle  

    // Sort tasks in ascending order
    sort(tasks.begin(), tasks.end());  // O(n log n)

    // Frequency counter map O(n log k)
    map<char, int> counter;
    for( char task : tasks){
        counter[task]++;
    }

    // Build priority queue for unique tasks based on frequency and task  
    priority_queue<pair<int, char>> pq;
    for(pair<char, int> task:counter){
        pq.push({task.second, task.first}); // O( k log k)  
    }
   
    // Initialize cooling queue and total intervals
    queue <pair<int,pair<int ,char>>> cooling; // Pair< available_remaining, <task,frequency>
    int totalIntervals = 0;
    int remaining = tasks.size();
    // Process tasks
    while(remaining > 0){ // n times
        totalIntervals++;
        // If queue not empty or cooling task is available to reinsert into pq 
        while(!cooling.empty() && cooling.front().first <= totalIntervals){
            pq.push(cooling.front().second);// Reinsert into pq with remaining frequency
            // Remove previous cooling entry
            cooling.pop();
        }

        if(!pq.empty()){
            pair<int, char> current = pq.top(); //O(1)
            int freq = current.first;
            char task = current.second;
            pq.pop();// O(log k)
            schedule.push_back(task);
            // Execute task
            freq--;
            remaining--;
            if(freq > 0){
                // Put in cooling {next available time for current task, remaining frequency}
                cooling.push({totalIntervals + n + 1, {freq, task}});// O(1) 
            }
        }else{
            schedule.push_back('_'); // Idle
        }

    }
    cout << "Schedule: ";
    for (char c : schedule) {
    if (c == '_') cout << "idle ";
    else cout << c << " ";
}
    return totalIntervals;
    
    // Total Complexity = O(n log k) + O(k log k) + O(n log k) = O(n log k)
    // Since k <=26 constant makes O(n)
    // Time = O(n) +  O(n log n) = O(n log n)
}

// =========================================================
// FACTORY FUNCTIONS (Required for Testing)
// =========================================================

extern "C" {
    PlayerTable* createPlayerTable() { 
        return new ConcretePlayerTable(); 
    }

    Leaderboard* createLeaderboard() { 
        return new ConcreteLeaderboard(); 
    }

    AuctionTree* createAuctionTree() { 
        return new ConcreteAuctionTree(); 
    }
}
