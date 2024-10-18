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
    bool has_failed = false;
    vector<int> channels = vector<int>(40, -1);
    unordered_set<int> services_set = {};

    // const bool has_capacity(int wavelength, int bandwidth) {
    //     for (int i = wavelength; i < wavelength + bandwidth; i++) {
    //         if (channels[i] != -1) return false;
    //     }
    //     return true;
    // }
};

struct Service {
    int id;
    Node source;
    Node dest;
    int Seq_length; // no és molt útlil ja que és igual a edges.size()
    int Left; // channel start 
    int Right; // channel end
    int Value;
    vector<Edge> path; // edge list
    

    /**
     * @return Right - Left + 1
     */
    const int bandwidth() const {
        return Right - Left + 1;
    }

    /**
     * @brief Dimension of the corresponding Hypergraph
     * @return 40 - bandwidth + 1
     */
    const int dim() const {
        return 40 - bandwidth() + 1;
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

            service.path.resize(service.Seq_length);
            for (Edge &edge : service.path) {
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
            // vector<vector<pair<Node, double>>> adj; // TODO: Store the edge to track its capacity, not only to the nodes it connects.
            vector<vector<Edge>> adj;
            static constexpr float CHANGE_COST = 100.0f;

        public:
            Graph() = default;

            Graph (const vector<Edge>& edges) : adj(input.N) {                
                for (const auto& edge : edges) {
                    adj[edge.u].push_back(edge);
                }
            }

            /**
             * @brief Checks if a range of wavelengths is available on a given edge for a service.
             *
             * This function iterates over a specified range of wavelengths on an edge and checks if they are available
             * for a given service. A wavelength is considered available if it is either unoccupied (-1) or already occupied
             * by the same service.
             *
             * @param edge The edge on which to check the availability of wavelengths.
             * @param wavelength The starting wavelength to check.
             * @param k The number of consecutive wavelengths to check.
             * @param serv The service for which the wavelengths are being checked.
             * @return true if all wavelengths in the specified range are available for the service, false otherwise.
             */
            bool is_wavelength_available(Edge edge, int wavelength, int k, Service serv) {
                for (int i = wavelength; i < wavelength + k; i++) {
                    if (edge.channels[i] != -1 and edge.channels[i] != serv.id) {
                        return false;
                    }
                }
                return true;
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
            float get_change_cost(Node start, Node end, Node u, Service Service) {
                if (u == start || u == end) return 0;
                return CHANGE_COST / (input.P[u] * Service.Value); // TODO: review and modify heuristic
                // TODO: make it so that for the last failed_edges have more probability of using the channel change (lower the cost)
                // we know that there are 60 fails at most, so: cost = K * log((61-i))
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

            /**
             */
            vector<HyperNode> dijkstra(Node start, Node end, Service& serv) { // TODO: get the start and end from the service
                int k = serv.dim();
                vector<vector<bool>> visited (input.N, vector<bool>(k, false));
                vector<vector<float>> distances (input.N, vector<float>(k, INF));
                vector<vector<HyperNode>> parent (input.N, vector<HyperNode>(k, {-1, -1}));
                
                auto cmp = [](const pair<HyperNode, float>& a, const pair<HyperNode, float>& b) {
                    return a.second > b.second;
                };
                priority_queue<pair<HyperNode, float>, vector<pair<HyperNode, float>>, decltype(cmp)> pqueue(cmp);

                distances[start][0] = 0;
                pqueue.push({{start, 0}, 0.0f});

                while (!pqueue.empty()) {
                    auto [hnode, dist] = pqueue.top();
                    auto [u, wavelength] = hnode;
                    pqueue.pop();
                    
                    if (visited[u][wavelength]) continue;
                    visited[u][wavelength] = true;
                    
                    if (hnode == HyperNode{end, 0}) break;
                    
                    for (const auto& edge : adj[u]) {
                        if (!is_wavelength_available(edge, wavelength, k, serv)) continue;

                        float new_dist = dist + 1; // WARNING +1 if we don't use heuritics, if not TODO: 
                        if (new_dist < distances[edge.v][wavelength]) {
                            distances[edge.v][wavelength] = new_dist;
                            parent[edge.v][wavelength] = {u, wavelength};
                            pqueue.push({{edge.v, wavelength}, new_dist});
                        }
                    }

                    for (int i = 0; i < k; i++) {
                        if (i == wavelength) continue;

                        float new_dist = dist + get_change_cost(start, end, u, serv);
                        if (new_dist < distances[u][i]) {
                            distances[u][i] = new_dist;
                            parent[u][i] = {u, wavelength};
                            pqueue.push({{u, i}, new_dist});
                        }
                    }
                    
                }

                if (distances[end][0] == INF) { // there is no path from start to end
                    cout << "No path found" << endl; // TODO !!!
                }

                return get_path_from_parents(parent, start, end); // TODO: update edges' channels & Graph's and update Pi's
            }

            /**
             * repeat for all dead services:
                * dijkstra
                * shortest path -> marquem com usats
             *
             * alliberar camins
             */

            void update_graph(Edge& edge_failed) {
                Node u = edge_failed.u;
                auto& edges = adj[u];
                edges.erase(std::remove_if(edges.begin(), edges.end(), [&](const std::pair<Node, int>& edge) {
                    return edge.first == edge_failed.v;
                }), edges.end());
            }
            
    };


    int T;
    vector<pair<Service, vector<Edge>>> replanned_services;
    Graph graph;

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
     * @brief Retrieves a sorted list of services affected by a failed edge.
     *
     * This function identifies the services that are affected by a given failed edge,
     * removes duplicates, and sorts them in @b descending order based on their value.
     *
     * @param failed_edge The edge that has failed, containing channel information.
     * @return A vector of affected services, sorted in descending order of their value.
     */
    vector<Service> get_affected_services_sorted(const Edge& failed_edge) {
        vector<Service> affected_services;
        unordered_set<int> seen_services_id;

        for (const int& id : failed_edge.channels) {
            if (id == -1 || seen_services_id.count(id)) continue;
            affected_services.push_back(input.services[id]);
            seen_services_id.insert(id);
        }
        sort(affected_services.begin(), affected_services.end(), [](const Service& a, const Service& b) {
            return a.Value > b.Value;
        });
        return affected_services;
    }


    void replan_failed_edge(const Edge& failed_edge) {
        // iterate over the affected services in ascending order of serv.Value
        vector<Service> affected_services = get_affected_services_sorted(failed_edge);
        
        for (auto& serv : affected_services) {
            graph.dijkstra(...);
        }
    }


    void replan_scenario() {
        int e_failed_idx;
        cin >> e_failed_idx;
        if (e_failed_idx == -1) {
            return;
        }

        const Edge& failed_edge = input.edges[e_failed_idx - 1];

        replan_failed_edge(failed_edge);
        print_replanned_services();
    }


    /**
     * @brief Handles the replanning process.
     *
     * Reads the number of test cases (T) and for each test case, reads the number of failed services (e_failed),
     * calls the replanned_services function, and prints the result.
     */
    void replanning() {
        cin >> T;
        graph = Graph(input.edges);

        for (int i = 0; i < T; i++)
        {
            replan_scenario();
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
