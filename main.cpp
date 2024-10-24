#include <algorithm>
#include <iostream>
#include <limits>
#include <queue>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>
#include <variant>

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

using EdgeIndex = int;

struct Edge {
    EdgeIndex idx;
    Node u;
    Node v;
    bool has_failed = false;
    vector<int> channels = vector<int>(40, -1);
    unordered_set<int> services_set = {};

    /**
     * @brief Returns the node at the other end of the edge.
     * 
     * Given a node `w`, this function returns the node at the opposite end of the edge.
     * If `w` is equal to `u`, it returns `v`; otherwise, it returns `u`.
     * 
     * @param w The node for which the opposite end is to be found.
     * @return Node The node at the other end of the edge.
     */
    Node get_other_end(Node w) const {
        if (w == u) return v;
        return u;
    }

    /**
     * @brief Retrieves the connecting node between the current edge and another edge.
     * 
     * @param other The other edge to compare with.
     * @return The connecting node if found, otherwise returns -1.
     */
    Node get_connecting_node(const Edge& other) const {
        if (u == other.u || u == other.v) return u;
        if (v == other.u || v == other.v) return v;
        return -1; // or throw an exception if no intersection is found
    }
};

struct WavelengthChangeEdge { // TODO: revise
    int wl_in;
    int wl_out;
};

using EdgeVariant = std::variant<Edge, WavelengthChangeEdge>;

struct EdgeWithWavelengths {
    Edge edge;
    int Left;
    int Right;
};


struct Service {
    int id;
    Node source;
    Node dest;
    int Seq_length; // no és molt útlil ja que és igual a edges.size()
    int Left; // channel start
    int Right; // channel end
    int Value;
    vector<EdgeWithWavelengths> path; // edge list
    bool is_dead = false;


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

        for (size_t i = 0; i < input.edges.size(); ++i) {
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
            for (EdgeWithWavelengths &edge_wl : service.path) {
                int ej;
                cin >> ej;
                edge_wl = {input.edges[ej-1], service.Left, service.Right};
            }
        }
    }
}


namespace Bottleneck {

    using Scenario = vector<Edge>;
    vector<Scenario> total_scenarios;

    // double jaccard_similarity(const std::vector<Edge>& set1, const std::vector<Edge>& set2) {
    //     if (set1.empty() && set2.empty()) return 1.0;
    //     if (set1.empty() || set2.empty()) return 0.0;

    //     size_t intersection_size = std::count_if(set1.begin(), set1.end(),
    //         [&set2](const Edge& elem) {
    //             return find(set2.begin(), set2.end(), elem) != set2.end();
    //         });

    //     size_t union_size = set1.size() + set2.size() - intersection_size;

    //     return static_cast<double>(intersection_size) / union_size;
    // }


    void print_bottleneck_results() {
        cout << 0 << endl;
    //     cout << total_scenarios.size() << endl;

    //     for (auto &scenario : total_scenarios)
    //     {
    //         cout << scenario.size() << endl;

    //         for (size_t i = 0; i < scenario.size(); ++i)
    //         {
    //             cout << scenario[i];
    //             if (i < scenario.size() - 1)
    //             {
    //                 cout << " ";
    //             }
    //         }
    //     }
    }

    void random_bottleneck() {

    }

    void bottleneck() {

        print_bottleneck_results();
    }
}


namespace Replanning {

    class Graph {

    private:
        vector<vector<Edge>> adj;
        static constexpr float K_CHANGE_COST = 100.0f;

    public:

        Graph () : adj(input.N) {
            for (const auto &edge : input.edges) {
                adj[edge.u].push_back(edge);
                adj[edge.v].push_back(edge);
            }
        }


        /**
         * @brief Checks if the specified wavelength range is available on the given edge for the service.
         * 
         * @param edge The edge on which to check the wavelength availability.
         * @param wavelength The starting wavelength to check.
         * @param serv The service requesting the wavelength.
         * @return true if the wavelength range is available for the service, false otherwise.
         */
        bool is_wavelength_available(const Edge& edge, int curr_wl, const Service& serv) {
            int W = serv.bandwidth();
            for (int i = curr_wl; i < curr_wl + W; i++) {
                if (input.edges[edge.idx].channels[i] != -1 and input.edges[edge.idx].channels[i] != serv.id) {
                    return false;
                }
            }
            return true;
        }


        /**
         * @brief Checks if an edge is valid for a given wavelength and service.
         * 
         * This function determines whether an edge is valid by checking if it has not failed
         * and if the wavelength is available for the given service.
         * 
         * @param edge The edge to be checked.
         * @param curr_wl The current wavelength to be checked.
         * @param serv The service for which the edge is being validated.
         * @return true if the edge is valid, false otherwise.
         */
        bool is_edge_valid(const Edge& edge, int curr_wl, const Service& serv) {
            return !input.edges[edge.idx].has_failed and is_wavelength_available(edge, curr_wl, serv);
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
        float get_wavelength_change_cost(Node u, const Service& serv) {
            if (u == serv.source || u == serv.dest) return 0;
            return K_CHANGE_COST / (input.P[u] * serv.Value); // TODO: review and modify heuristic
            // TODO: make it so that for the last failed_edges have more probability of using the channel change (lower the cost)
            // we know that there are 60 fails at most, so: cost = K * log((61-i))
        }


        /**
         * @brief Constructs a path of edges with wavelengths from the parent edge information.
         *
         * This function generates a path from the destination node to the source node based on the parent edge information
         * stored in a 2D vector. The path is constructed by tracing back from the destination node to the source node,
         * updating the node and wavelength at each step. The resulting path is a vector of edges with their corresponding
         * wavelengths.
         *
         * @param parent_edge A 2D vector containing the parent edge information for each node and wavelength.
         * @param serv The service object containing the source and destination nodes, as well as the bandwidth.
         * @return A vector of EdgeWithWavelengths representing the path from the source to the destination node.
         */
        vector<EdgeWithWavelengths> get_path_from_parents(const vector<vector<EdgeVariant>> &parent_edge, const Service &serv) {
            // TODO: review logic for the beggining / end for the path
            vector<EdgeWithWavelengths> path;
            const int W = serv.bandwidth();

            Node node = serv.dest;
            int wavelength = 0;
            EdgeVariant parent_edge_variant;

            while (node != serv.source) {
                parent_edge_variant = parent_edge[node][wavelength];

                std::visit([&node, &wavelength, &path, &W](const auto& edge) { // TODO: add Pi and channel update here if code is slow
                    using T = std::decay_t<decltype(edge)>;
                    
                    if constexpr (std::is_same_v<T, Edge>) {
                        path.push_back({edge, wavelength, wavelength + W - 1});
                        node = edge.get_other_end(node);
                    } else if constexpr (std::is_same_v<T, WavelengthChangeEdge>) {
                        wavelength = edge.wl_in;
                    }
                    
                }, parent_edge_variant);
            }
            reverse(path.begin(), path.end());

            return path;
        }

        /**
         * @brief Modified Dijkstra that prevents visiting the same node multiple times across all wavelengths
         */
        vector<EdgeWithWavelengths> find_shortest_path(const Service &serv) {
            int k = serv.dim();
            vector<vector<bool>> visited(input.N, vector<bool>(k, false));
            vector<vector<float>> distances(input.N, vector<float>(k, INF));
            vector<vector<EdgeVariant>> parent_edge(input.N, vector<EdgeVariant>(k));
            
            // Track which nodes have been used in the current best path to each node
            vector<vector<unordered_set<Node>>> used_nodes_to_node(input.N, vector<unordered_set<Node>>(k));

            auto cmp = [](const pair<HyperNode, float>& a, const pair<HyperNode, float>& b) {
                return a.second > b.second;
            };
            priority_queue<pair<HyperNode, float>, vector<pair<HyperNode, float>>, decltype(cmp)> pqueue(cmp);

            // Initialize starting state
            distances[serv.source][0] = 0.0f;
            pqueue.push({{serv.source, 0}, 0.0f});
            used_nodes_to_node[serv.source][0].insert(serv.source);  // Mark source as used

            while (!pqueue.empty()) {
                auto [hnode, dist] = pqueue.top();
                auto [curr_node, curr_wl] = hnode;
                pqueue.pop();

                if (visited[curr_node][curr_wl]) continue;
                visited[curr_node][curr_wl] = true;

                if (hnode == HyperNode{serv.dest, 0}) break;

                // Process edges in current wavelength
                for (const auto& edge : adj[curr_node]) {
                    Node next_node = edge.get_other_end(curr_node);

                    // Skip if this node has been used in any wavelength in the path to current node
                    if (used_nodes_to_node[curr_node][curr_wl].count(next_node)) continue;

                    if (!is_edge_valid(edge, curr_wl, serv)) continue;

                    float next_dist = dist + 1;

                    // Only update if we find a better path AND the node hasn't been used
                    if (next_dist < distances[next_node][curr_wl]) {
                        // Create new set of used nodes for this path
                        auto new_used_nodes = used_nodes_to_node[curr_node][curr_wl];
                        new_used_nodes.insert(next_node);

                        // Update the path
                        distances[next_node][curr_wl] = next_dist;
                        parent_edge[next_node][curr_wl] = edge;
                        used_nodes_to_node[next_node][curr_wl] = new_used_nodes;
                        pqueue.push({{next_node, curr_wl}, next_dist});
                    }
                }

                // Process wavelength changes
                if (input.P[curr_node] == 0) continue;

                for (int next_wl = 0; next_wl < k; next_wl++) {
                    if (next_wl == curr_wl) continue;

                    float next_dist = dist + get_wavelength_change_cost(curr_node, serv);

                    if (next_dist < distances[curr_node][next_wl]) {
                        // When changing wavelength, we keep the same used nodes
                        distances[curr_node][next_wl] = next_dist;
                        parent_edge[curr_node][next_wl] = WavelengthChangeEdge{curr_wl, next_wl};
                        used_nodes_to_node[curr_node][next_wl] = used_nodes_to_node[curr_node][curr_wl];
                        pqueue.push({{curr_node, next_wl}, next_dist});
                    }
                }
            }

            if (distances[serv.dest][0] == INF) {
                return {};
            }

            return get_path_from_parents(parent_edge, serv);
        }
    };



    int T;
    // vector<pair<Service, vector<Edge>>> replanned_services;
    vector<Edge> failed_edges; // keep track of failed edges to revert after scenario change
    Graph graph;


void print_replanned_services(const vector<vector<EdgeWithWavelengths>>& paths, const vector<Service>& affected_services) {
        
        int num_Replanned_services = 0;
        for (const auto& path : paths) {
            if (!path.empty()) num_Replanned_services++;
        }
        cout << num_Replanned_services << '\n';

        // cout << "replanned successfully: " << num_Replanned_services << " / " << paths.size() << endl;

        for (size_t i = 0; i < paths.size(); i++) {
            if (paths[i].empty()) continue;
            vector<EdgeWithWavelengths>path = paths[i];
            Service serv = affected_services[i];
            cout << serv.id+1 << ' ' << path.size() << '\n';

            for (const auto &[edge, Left, Right] : path) {
                cout << edge.idx+1 << ' ' << Left+1 << ' ' << Right+1 << ' ';
            }
            // // the same but without the trailing whitespace
            // for (size_t j = 0; j < path.size(); ++j) {
            //     const auto &[edge, Left, Right] = path[j];
            //     cout << edge.idx+1 << ' ' << Left+1 << ' ' << Right+1;
            //     if (j < path.size() - 1) {
            //         cout << ' ';
            //     }
            // }
            cout << endl;
            fflush(stdout);
        }

    }


    /**
     * @brief Updates the channels of the edges in the given path and adjusts the Pi values of the nodes.
     *
     * This function iterates through the provided path of edges and updates the channels for each edge
     * based on the given service. Additionally, it adjusts the Pi values of the nodes when necessary.
     *
     * @param path A vector of EdgeWithWavelengths representing the path of edges to be updated.
     * @param serv A constant reference to a Service object containing the service information.
     */
    void update_Pis_and_channels(const vector<EdgeWithWavelengths>& path, const Service& serv) {
        if (path.empty()) return;

        int previous_Right = path[0].Right;
        Edge previous_edge;

        for (const auto &[edge, Left, Right] : path) {
            for (int i = Left; i <= Right; i++) { // update channels
                input.edges[edge.idx].channels[i] = serv.id;
            }
            if (Right != previous_Right) { // update Pi
                Node node = edge.get_connecting_node(previous_edge);
                input.P[node] -= 1;
            }
            previous_edge = edge;
            previous_Right = Right;
        }

    }


    /**
     * @brief Frees old used channels by comparing old service paths with new paths.
     *
     * This function iterates through the old service paths and compares them with the new paths.
     * It performs the following operations:
     * - If an edge in the old path is not present in the new path, it marks the channels of that edge as free (-1).
     * - If an edge is present in both old and new paths, it updates the channels based on the new path's wavelength.
     *
     * @param new_paths A vector of vectors containing the new paths with edges and their wavelengths.
     * @param old_services A vector of services, each containing a path with edges and their wavelengths.
     */
    void free_old_used_channels(const vector<vector<EdgeWithWavelengths>>& new_paths, const vector<Service>& old_services) {
        //old_path=serv.path and new_path
        //els que estan en el old i no en el new, eliminem
        //els que estan als dos quedar-me amb la wave length del nou i posar -1 al antic

        //TODO: change everything to reference in order to change old_edges
        //TODO: change service.path to type EdgeWithWavelengths and then optimize delete channel and change wavelenght

        for (size_t i = 0; i < old_services.size(); i++) {

            for (const auto& [old_edge, old_Left, old_Right] : old_services[i].path) {
                bool is_edge_in_new_path = false; // Flag to track if a new edge is found

                for (const auto &[new_edge, new_Left, new_Right] : new_paths[i]) {
                    if (old_edge.idx == new_edge.idx) {
                        is_edge_in_new_path = true; // Set the flag to true if a matching new edge is found
                        // Change wavelength if needed
                        for (int j = old_Left; j <= min(new_Left + 1, old_Right); j++) {
                            input.edges[old_edge.idx].channels[j] = -1;
                        }
                        for (int j = max(new_Right + 1, old_Left); j <= old_Right; j++) {
                            input.edges[old_edge.idx].channels[j] = -1;
                        }
                        break; // Exit the loop early since a matching new edge is found
                    }
                }

                if (!is_edge_in_new_path && !new_paths[i].empty()) {
                    // Delete channel (i.e., change to -1 our service)
                    for (int j = old_Left; j <= old_Right; j++) {
                        input.edges[old_edge.idx].channels[j] = -1;
                    }
                }

            }
        }
    }


    void update_services_path(const vector<vector<EdgeWithWavelengths>>& new_paths, const vector<Service>& old_services) {
        for (size_t i = 0; i < new_paths.size(); i++) {
            input.services[old_services[i].id].path = new_paths[i];
        }
    }


    /**
     * @brief Retrieves a sorted list of services affected by a failed edge.
     *
     * This function identifies the services affected by a given failed edge,
     * filters out any services that are already marked as dead, and sorts
     * the remaining services in descending order based on their value.
     *
     * @param failed_edge The edge that has failed, containing channel information.
     * @return A vector of affected services, sorted by their value in descending order.
     */
    vector<Service> get_affected_services_sorted(const Edge& failed_edge) {
        vector<Service> affected_services;
        unordered_set<int> seen_services_id;

        for (const int& serv_id : failed_edge.channels) {
            if (serv_id == -1 || seen_services_id.count(serv_id) || input.services[serv_id].is_dead) continue;

            affected_services.push_back(input.services[serv_id]);
            seen_services_id.insert(serv_id);
        }
        sort(affected_services.begin(), affected_services.end(), [](const Service& a, const Service& b) {
            return a.Value > b.Value;
        });
        // for (const auto& service : affected_services) {
        //     cout << "Service ID: " << service.id + 1 << "\n";
        //     cout << "Value: " << service.Value << "\n";
        //     cout << "\n";
        // }
        return affected_services;
    }


    /**
     * @brief Replans the network paths for services affected by a failed edge.
     *
     * This function handles the re-routing of services that are impacted by a failure in a specific edge.
     * It first retrieves the affected services sorted in ascending order of their service value.
     * For each affected service, it finds the shortest path in the graph and updates the Pis and channels if a valid path is found.
     * Finally, it updates the dead channels based on the newly computed shortest paths.
     *
     * @param failed_edge The edge that has failed and needs re-routing.
     */
    void replan_failed_edge(const Edge& failed_edge) {
        const vector<Service>& affected_services = get_affected_services_sorted(failed_edge);
        vector<vector<EdgeWithWavelengths>> shortest_paths;

        for (const Service &serv : affected_services) {
            vector<EdgeWithWavelengths> path = graph.find_shortest_path(serv);
            update_Pis_and_channels(path, serv);
            // cout << "replanning Serv.id " << serv.id+1 << " was succes: " << !path.empty() << endl;
            shortest_paths.push_back(path);
            // TODO: remove the empty paths and their corresponding affected_service
        }

        free_old_used_channels(shortest_paths, affected_services);

        update_services_path(shortest_paths, affected_services);

        print_replanned_services(shortest_paths, affected_services);
    }


    void replan_scenario() {
        int e_failed_idx;

        // cout << "\nFailed edges" << endl;
        // for (const auto& edge : input.edges) {
        //     cout << "Edge ID: " << edge.idx + 1 << "\n";
        //     cout << "Nodes: " << edge.u + 1 << " - " << edge.v + 1 << "\n";
        //     cout << "Used Channels: ";
        //     cout << "has_failed? " << edge.has_failed << '\n';
        //     for (size_t i = 0; i < edge.channels.size(); ++i) {
        //         cout << edge.channels[i] << " ";
        //     }
        //     cout << "\n\n";
        // }

        while (cin >> e_failed_idx) {
            if (e_failed_idx == -1) return;

            Edge& failed_edge = input.edges[e_failed_idx - 1];
            failed_edge.has_failed = true;
            failed_edges.push_back(failed_edge);

            replan_failed_edge(failed_edge);
            
            // cout << "\nFailed edges" << endl;
            // for (const auto& edge : input.edges) {
            //     cout << "Edge ID: " << edge.idx + 1 << "\n";
            //     cout << "Nodes: " << edge.u + 1 << " - " << edge.v + 1 << "\n";
            //     cout << "Used Channels: ";
            //     cout << "has_failed? " << edge.has_failed << '\n';
            //     for (size_t i = 0; i < edge.channels.size(); ++i) {
            //         cout << edge.channels[i] << " ";
            //     }
            //     cout << "\n\n";
            // }


        }
    }

    /**
     * @brief Resets the scenario to its initial state.
     * 
     * This function resets the input state to the provided initial state and clears
     * the list of failed edges in preparation for the next scenario.
     * 
     * @param initial_input_state The initial state to reset the input to.
     */
    void reset_scenario(const Input initial_input_state) {
        input = initial_input_state;
        failed_edges = {}; // reset failed edges for the next scenario
    }

    void setup_edges_channels() {
        for (Service &serv : input.services) {
            for (auto &[edge, L, R] : serv.path) {

                    for (int i = L ; i <= R; i++) {
                    input.edges[edge.idx].channels[i] = serv.id;
                }

            }
        }
    }

    /**
     * @brief Handles the replanning process.
     *
     * Reads the number of test cases (T) and for each test case, reads the number of failed services (e_failed),
     * calls the replanned_services function, and prints the result.
     */
    void replanning() {
        cin >> T;
        graph = Graph();

        // save initial services paths
        setup_edges_channels();
        failed_edges = {};
        Input initial_input_state = input;

        for (int i = 0; i < T; i++) {
            // for (const auto& service : input.services) {
            //     cout << "Service ID: " << service.id + 1 << "\n";
            //     cout << "Source: " << service.source + 1 << "\n";
            //     cout << "Destination: " << service.dest + 1 << "\n";
            //     cout << "Sequence Length: " << service.Seq_length << "\n";
            //     cout << "Left: " << service.Left + 1 << "\n";
            //     cout << "Right: " << service.Right + 1 << "\n";
            //     cout << "Value: " << service.Value << "\n";
            //     cout << "Path: ";
            //     for (const auto& edge_wl : service.path) {
            //     cout << "(" << edge_wl.edge.idx + 1 << ", " << edge_wl.Left + 1 << "-" << edge_wl.Right + 1 << ") ";
            //     }
            //     cout << "\n\n";
            // }
            // for (size_t i = 0; i < input.P.size(); i++) {
            //     cout << "p-" << i << ' ' << input.P[i] << endl;
                
            // }
            // cout << "failed-edges:" << endl;
            // for (const Edge &edge: failed_edges) {
            //     cout << edge.idx << endl;
            // }

            replan_scenario();
            reset_scenario(initial_input_state);
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
