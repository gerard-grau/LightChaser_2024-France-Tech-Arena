#include <vector>
#include <iostream>
#include <utility>
#include <tuple>

using namespace std;

using Node = int;

struct Edge {
    int idx;
    Node u;
    Node v;
};

struct Service {
    int id;
    Node s;
    Node d;
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


namespace ReadInput {
    void read_input() {
        cin >> input.N >> input.M;
        input.P.resize(input.N);
        input.edges.resize(input.M);

        for (int& Pi : input.P) {
            cin >> Pi;
        }

        for (int i = 0; i < input.edges.size(); ++i) {
            input.edges[i].idx = i;
        
            Node ui, vi;
            cin >> ui >> vi;
            input.edges[i].u = ui - 1;
            input.edges[i].v = vi - 1;
        }


        cin >> input.J;
        input.services.resize(input.J);

        for (int i = 0; i < input.J; i++) {
            Service& service = input.services[i];
            service.id = i;

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
}


namespace Bottleneck {

    using Scenario = vector<int>;
    vector<Scenario> total_scenarios;


    void print_bottleneck_results() {
        
        cout << total_scenarios.size() << endl;

        for (auto& scenario : total_scenarios) {
            cout << scenario.size() << endl;

            for (size_t i = 0; i < scenario.size(); ++i) {
                cout << scenario[i];
                if (i < scenario.size() - 1) {
                    cout << " ";
                }
            }
        }

    }


    void bottleneck() {
        
        print_bottleneck_results();
    }
}


namespace Replanning {
    
    int T;
    using EdgeWithChannel = tuple<Edge, int, int>;
    vector<pair<Service, vector<EdgeWithChannel>>> replanned_services;

    void replan_services(const Edge& e_failed) {

    }

    void print_replanned_services() {
        cout << replanned_services.size() << endl;
        for (auto& [service, serv_edges] : replanned_services) {
            cout << service.id << ' ' << serv_edges.size() << endl;

            for (auto& [edge, e_l, e_r] : serv_edges) {
                cout << edge.idx << e_l << e_r;
            }
            cout << endl;
        }
        fflush(stdout);
    }

    /**
     * @brief Handles the replanning process.
     * 
     * Reads the number of test cases (T) and for each test case, reads the number of failed services (e_failed),
     * calls the replanned_services function, and prints the result.
     */
    void replanning() {
        cin >> T;

        for (int i = 0; i < T; i++) {
            int e_failed_idx;
            cin >> e_failed_idx;
            if (e_failed_idx == -1) return;

            const Edge& e_failed = input.edges[e_failed_idx-1];
            
            replan_services(e_failed);
            print_replanned_services();
        }
    }

    
}


int main() {
    ReadInput::read_input();
    
    Bottleneck::bottleneck();
    
    Replanning::replanning();

    return 0;
}
