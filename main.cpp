#include <vector>
#include <iostream>

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

Input input;


void read_input() {
    cin >> input.N >> input.M;
    input.P.resize(input.N);
    input.edges.resize(input.M);

    for (int& Pi : input.P) {
        cin >> Pi;
    }
    for (Edge& edge : input.edges) {
        int ui, vi;
        cin >> ui >> vi;
        edge.u = ui - 1;
        edge.v = vi - 1;
    }

    cin >> input.J;
    input.services.resize(input.J);

    for (int i = 0; i < input.J; i++) {
        Service& service = input.services[i];
        int s, d, L, R;
        cin >> s >> d >> service.S >> L >> R >> service.V;
        service.s = s - 1;
        service.d = d - 1;
        service.L = L - 1;
        service.R = R - 1;

        service.edges.resize(service.S);
        for (Edge& edge : service.edges) {
            int ej;
            cin >> ej;
            edge = input.edges[ej-1];
        }
    }
}


namespace Bottleneck {

    vector<vector<int>> failures;

    void print_bottleneck_results() {
        cout << failures.size() << endl;
        for (vector<int> failure : failures) {
            cout << failure.size() << endl;

            for (size_t i = 0; i < failure.size(); ++i) {
                cout << failure[i];
                if (i < failure.size() - 1) {
                    cout << " ";
                }
            }
        }
        fflush(stdout); // TODO: REVISAR !!!

    }

    void Bottleneck() {
        
        print_bottleneck_results();
    }
}


namespace Replanning {

    void Replanning() {
        
    }
}


int main() {
    read_input();
    
    Bottleneck::Bottleneck();
    
    Replanning::Replanning();

    return 0;
}
