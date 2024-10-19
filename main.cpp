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
};

struct WavelengthChangeEdge { // TODO:  revise
    int wl_in;
    int wl_out;
};

struct EdgeWithWavelengths {
    Edge edge;
    int Left;
    int Right;
};

// struct WavelengthChange : public Edge {
//     int idx = -1;
//     int wavelength_in;
//     int wavelength_out;
// };

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
            bool is_wavelength_available(const Edge& edge, int wavelength, const Service& serv) {
                int k = serv.dim();
                for (int i = wavelength; i < wavelength + k; i++) {
                    if (edge.channels[i] != -1 and edge.channels[i] != serv.id) {
                        return false;
                    }
                }
                return true;
            }

            /**
             * @brief Calculates the cost of changing the channel for a given node and service.
             *
             * This function computes the cost associated with changing the channel for a node `u` 
             * in the context of a specific service `serv`. If the node `u` is either the source 
             * or the destination of the service, the cost is zero. Otherwise, the cost is 
             * determined by a heuristic formula that takes into account the node's probability 
             * and the service's value.
             *
             * @param u The node for which the channel change cost is being calculated.
             * @param serv The service context which includes the source, destination, and value.
             * @return The cost of changing the channel for the node `u`.
             */
            float get_change_cost(Node u, const Service& serv) {
                if (u == serv.source || u == serv.dest) return 0;
                return CHANGE_COST / (input.P[u] * serv.Value); // TODO: review and modify heuristic
                // TODO: make it so that for the last failed_edges have more probability of using the channel change (lower the cost)
                // we know that there are 60 fails at most, so: cost = K * log((61-i))
            }


            /**
             * @brief Constructs a path from the parent nodes for a given service.
             *
             * This function traces back from the destination node to the source node using the parent nodes
             * and constructs the path taken. It also considers the wavelength used for the service.
             *
             * @param parent A 2D vector containing the parent nodes for each node.
             * @param serv A reference to the Service object containing source and destination nodes.
             * @return A vector of HyperNode representing the path from source to destination.
             */
            vector<HyperNode> get_path_from_parents(vector<vector<HyperNode>> parent, const Service& serv) { // TODO: REMOVE THIS
                vector<HyperNode> path;
                Node v = serv.dest;
                int wavelength = parent[serv.dest][0].wavelength; // Use parent's wavelength if different, if equal this doesn't matter

                while (v != serv.source) {
                    path.push_back({v, wavelength});
                    auto [v, wavelength] = parent[serv.dest][0]; // TODO: THIS IS WRONG!
                }
                path.push_back({serv.dest, wavelength});
                reverse(path.begin(), path.end());
                

                return path; // TODO: update edges' channels & Graph's and update Pi's
            }

            vector<EdgeWithWavelengths> get_path_from_parents(const vector<vector<Edge>>& incoming_edge, const Service& serv) { // TODO: convert this to edges instead of Nodes
                // TODO: review logic for the beggining / end for the path
                vector<EdgeWithWavelengths> path;
                // u → v
                Node v = serv.dest;
                int wavelength = 0;
                Edge edge = incoming_edge[v][wavelength];
                Node u = edge.u;
                
                // TODO: check if edge is a WavelengthChange, if so update u and v (so that the end is at the appropiate wavelength, and not 0)
                
                while (u != serv.source) { // TODO: look for a nice way to store wavelength changes as Edges
                    edge = incoming_edge[edge.v][wavelength];
                    u = edge.u;
                    v = edge.v;
                    path.push_back({edge, wavelength, wavelength + serv.bandwidth()-1}); // TODO: replace with the comment below

                    // if (...)  { //regular edge
                    //     wavelength_path.push_back({edge, wavelength, wavelength + serv.bandwidth()-1});
                    // }
                    // else { // wavelength change
                    //     int wavelength = edge.wl_out
                    // }
                }
                reverse(path.begin(), path.end());

                return path; // TODO: update edges' channels and Graph's, and update Pi's
            }


            /**
             */
            vector<EdgeWithWavelengths> dijkstra(const Service& serv) {
                int k = serv.dim();
                vector<vector<bool>> visited (input.N, vector<bool>(k, false));
                vector<vector<float>> distances (input.N, vector<float>(k, INF));
                vector<vector<Edge>> incoming_edge (input.N, vector<Edge>(k));
                
                auto cmp = [](const pair<HyperNode, float>& a, const pair<HyperNode, float>& b) {
                    return a.second > b.second;
                };
                priority_queue<pair<HyperNode, float>, vector<pair<HyperNode, float>>, decltype(cmp)> pqueue(cmp);

                distances[serv.source][0] = 0.0f;
                pqueue.push({{serv.source, 0}, 0.0f});

                while (!pqueue.empty()) {
                    auto [hnode, dist] = pqueue.top();
                    auto [u, curr_wl] = hnode;
                    pqueue.pop();

                    if (visited[u][curr_wl]) continue;
                    visited[u][curr_wl] = true;
                    
                    if (hnode == HyperNode{serv.dest, 0}) break; // Ensure the algorithm reaches [dest, 0] to correctly construct the path later

                    
                    for (const auto& edge : adj[u]) {
                        if (!is_wavelength_available(edge, curr_wl, serv)) continue;

                        float next_dist = dist + 1; // WARNING +1 if we don't use heuritics, if not TODO: 
                        if (next_dist < distances[edge.v][curr_wl]) {
                            distances[edge.v][curr_wl] = next_dist;
                            // parent[edge.v][wavelength] = {u, wavelength};
                            incoming_edge[edge.v][curr_wl] = edge;
                            pqueue.push({{edge.v, curr_wl}, next_dist});
                        }
                    }

                    for (int next_wl = 0; next_wl < k; next_wl++) { // Wavelength change
                        if (next_wl == curr_wl) continue;

                        float next_dist = dist + get_change_cost(u, serv);
                        if (next_dist < distances[u][next_wl]) {
                            distances[u][next_wl] = next_dist;
                            // parent[u][next_wl] = {u, wavelength};
                            incoming_edge[u][next_wl] = Edge{-1, u, u};
                            // incoming_edge[u][next_wl] = nextObj({u, curr_wl, next_wl}); // TODO: store wavelength in edge change
                            pqueue.push({{u, next_wl}, next_dist});
                        }
                    }
                    
                }

                if (distances[serv.dest][0] == INF) { // there is no path from start to end
                    cout << "No path found" << endl; // TODO !!!
                    return {}; // return empty vector
                }

                return get_path_from_parents(incoming_edge, serv); // TODO: update edges' channels & Graph's and update Pi's
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
