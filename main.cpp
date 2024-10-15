#include <algorithm>
#include <iostream>
#include <limits>
#include <queue>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;

const float INF = numeric_limits<float>::infinity();

using Node = int;

struct HyperNode {
    Node node;
    int wavelength;

    bool operator==(const HyperNode& other) const {
        return node == other.node && wavelength == other.wavelength;
    }
};

struct Service;

struct Edge {
    int idx;
    Node u;
    Node v;
    vector<int> channels = vector<int>(40, -1);
    unordered_set<int> services_set;

    const bool has_capacity(int wavelength, int bandwidth) {
        for (int i = wavelength; i < wavelength + bandwidth; i++) {
            if (channels[i] != -1) return false;
        }
        return true;
    }
};

struct Service {
    int id;
    Node source;
    Node dest;
    int Seq_length; // no és molt útlil ja que és igual a edges.size()
    int Left; // channel start 
    int Right; // channel end
    int Value;
    vector<Edge> edges;

    /**
     * @return Right - Left + 1
     */
    const int bandwidth() const {
        return Right - Left + 1;
    }
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

        for (int &Pi : input.P) {
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
            Service &service = input.services[i];
            service.id = i;

            int s, d, L, R;
            cin >> s >> d >> service.Seq_length >> L >> R >> service.Value;
            service.source = s - 1;
            service.dest = d - 1;
            service.Left = L - 1;
            service.Right = R - 1;

            service.edges.resize(service.Seq_length);
            for (Edge &edge : service.edges) {
                int ej;
                cin >> ej;
                edge = input.edges[ej - 1];
            }
        }
    }
}

namespace Bottleneck {

    using Scenario = vector<Edge>;
    vector<Scenario> total_scenarios;

    double jaccard_similarity(const std::vector<Edge>& set1, const std::vector<Edge>& set2) {
        if (set1.empty() && set2.empty()) return 1.0;
        if (set1.empty() || set2.empty()) return 0.0;

        size_t intersection_size = std::count_if(set1.begin(), set1.end(),
            [&set2](const Edge& elem) {
                return find(set2.begin(), set2.end(), elem) != set2.end();
            });

        size_t union_size = set1.size() + set2.size() - intersection_size;

        return static_cast<double>(intersection_size) / union_size;
    }


    void print_bottleneck_results()
    {

        cout << total_scenarios.size() << endl;

        for (auto &scenario : total_scenarios)
        {
            cout << scenario.size() << endl;

            for (size_t i = 0; i < scenario.size(); ++i)
            {
                cout << scenario[i];
                if (i < scenario.size() - 1)
                {
                    cout << " ";
                }
            }
        }
    }

    void random_bottleneck() {

    }

    void bottleneck() {
        
        print_bottleneck_results();
    }
}

namespace Replanning
{
    
    class Graph {
        
        private:
            vector<vector<pair<Node, double>>> adj; // TODO: Store the edge to track its capacity, not only to the nodes it connects.

        public:

            Graph (const vector<Edge>& edges) : adj(input.N) {                
                for (const auto& edge : edges) {
                    adj[edge.u].push_back(pair(edge.v, 1));
                }
            }

            /**
             * @brief Calculates the cost of changing a service at a given node.
             *
             * This function computes the cost associated with changing a service at a node `u` 
             * between the start and end nodes. If the node `u` is either the start or end node, 
             * the cost is zero. Otherwise, the cost is calculated based on the inverse of the 
             * product of the input parameter at node `u` and the service value.
             *
             * @param start The starting node.
             * @param end The ending node.
             * @param u The node at which the service change cost is being calculated.
             * @param Service The service for which the change cost is being calculated.
             * @return The cost of changing the service at node `u`.
             */
            float change_cost(Node start, Node end, Node u, Service Service) {
                if (u == start || u == end) return 0;
                return 100.0f / (input.P[u] * Service.Value); // TODO: review and modify heuristic
            }

            /**
             * @brief Constructs a path from the parent nodes.
             *
             * This function generates a path from the given parent nodes starting from the end node and tracing back to the start node.
             * It constructs the path by following the parent nodes and considering the wavelength of the edges.
             *
             * @param parent A 2D vector of HyperNode objects representing the parent nodes.
             * @param start The starting node of the path.
             * @param end The ending node of the path.
             * @return A vector of HyperNode objects representing the path from the start node to the end node.
             */
            vector<HyperNode> get_path_from_parents(vector<vector<HyperNode>> parent, Node start, Node end) {

                vector<HyperNode> path;
                Node v = end;
                int wavelength = parent[end][0].wavelength; // Use parent's wavelength if different, if equal this doesn't matter

                while (v != start) {
                    path.push_back({v, wavelength});
                    auto [v, wavelength] = parent[end][0];
                }
                path.push_back({start, wavelength});
                reverse(path.begin(), path.end());

                return path; // TODO: update edges' channels & Graph's and update Pi's
            }

            vector<HyperNode> dijkstra(Node start, Node end, Service& serv) {
                int k = 40 - serv.bandwidth() + 1;

                vector<vector<bool>> visited (input.N, vector<bool>(k, false));
                vector<vector<float>> dist (input.N, vector<float>(k, INF));
                vector<vector<HyperNode>> parent (input.N, vector<HyperNode>(k, {-1, -1}));
                
                auto cmp = [](const pair<HyperNode, float>& a, const pair<HyperNode, float>& b) {
                    return a.second > b.second;
                };
                priority_queue<pair<HyperNode, float>, vector<pair<HyperNode, float>>, decltype(cmp)> pqueue(cmp);

                dist[start][0] = 0;
                pqueue.push({{start, 0}, 0.0f});

                while (!pqueue.empty()) {
                    auto [hnode, cost] = pqueue.top();
                    auto [u, wavelength] = hnode;
                    pqueue.pop();
                    
                    if (visited[u][wavelength]) continue;
                    visited[u][wavelength] = true;
                    
                    if (hnode == HyperNode{end, 0}) break;
                    
                    for (const auto& [v, c] : adj[u]) {
                        // TODO: check if the edge [u, v] has capacity for this wavelength
                        // if (input.edges[u][v][wavelength:wavelength+r] != [-1..-1]) continue
                        float new_cost = cost + c;
                        if (new_cost < dist[v][wavelength]) {
                            dist[v][wavelength] = new_cost;
                            parent[v][wavelength] = {u, wavelength};
                            pqueue.push({{v, wavelength}, new_cost});
                        }
                    }

                    for (int i = 0; i < k; i++) {
                        if (i == wavelength) continue;

                        float new_cost = cost + change_cost(start, end, u, serv);
                        if (new_cost < dist[u][i]) {
                            dist[u][i] = new_cost;
                            parent[u][i] = {u, wavelength};
                            pqueue.push({{u, i}, new_cost});
                        }
                    }
                    
                }

                if (dist[end][0] == INF) { // there is no path from start to end
                    cout << "No path found" << endl; // TODO !!!
                }

                return get_path_from_parents(parent, start, end); // TODO: update edges' channels & Graph's and update Pi's
            }
    };


    int T;
    vector<pair<Service, vector<Edge>>> replanned_services;

    void replan_services(const Edge &e_failed) {

    }

    void print_replanned_services() {
        cout << replanned_services.size() << endl;
        for (auto &[service, serv_edges] : replanned_services)
        {
            cout << service.id << ' ' << serv_edges.size() << endl;

            for (auto &[edge, e_l, e_r] : serv_edges)
            {
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

        for (int i = 0; i < T; i++)
        {
            int e_failed_idx;
            cin >> e_failed_idx;
            if (e_failed_idx == -1)
                return;

            const Edge &e_failed = input.edges[e_failed_idx - 1];

            replan_services(e_failed);
            print_replanned_services();
        }
    }

}


int main()
{
    ReadInput::read_input();

    Bottleneck::bottleneck();

    Replanning::replanning();

    return 0;
}
