#include <vector>

using namespace std;

struct Edge {
    int u;
    int v;
};

struct Service {
    int s;
    int d;
    int S;
    int L;
    int R;
    int V;
    vector<Edge> edges;
};

struct Input {
    int N;
    int M;
    vector<int> P;
    vector<Edge> edges;
    int J;
    vector<Service> services;
};

