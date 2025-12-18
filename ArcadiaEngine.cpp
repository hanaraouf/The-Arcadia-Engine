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

public:
    ConcretePlayerTable() {
        // TODO: Initialize your hash table
    }

    void insert(int playerID, string name) override {
        // TODO: Implement double hashing insert
        // Remember to handle collisions using h1(key) + i * h2(key)
    }

    string search(int playerID) override {
        // TODO: Implement double hashing search
        // Return "" if player not found
        return "";
    }
};

// --- 2. Leaderboard (Skip List) ---

class ConcreteLeaderboard : public Leaderboard {
private:
    // TODO: Define your skip list node structure and necessary variables
    // Hint: You'll need nodes with multiple forward pointers

public:
    ConcreteLeaderboard() {
        // TODO: Initialize your skip list
    }

    void addScore(int playerID, int score) override {
        // TODO: Implement skip list insertion
        // Remember to maintain descending order by score
    }

    void removePlayer(int playerID) override {
        // TODO: Implement skip list deletion
    }

    vector<int> getTopN(int n) override {
        // TODO: Return top N player IDs in descending score order
        return {};
    }
};

// --- 3. AuctionTree (Red-Black Tree) ---

class ConcreteAuctionTree : public AuctionTree {
private:
    // TODO: Define your Red-Black Tree node structure
    // Hint: Each node needs: id, price, color, left, right, parent pointers

public:
    ConcreteAuctionTree() {
        // TODO: Initialize your Red-Black Tree
    }

    void insertItem(int itemID, int price) override {
        // TODO: Implement Red-Black Tree insertion
        // Remember to maintain RB-Tree properties with rotations and recoloring
    }

    void deleteItem(int itemID) override {
        // TODO: Implement Red-Black Tree deletion
        // This is complex - handle all cases carefully
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
