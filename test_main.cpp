// test_main.cpp
#include <iostream>
#include <vector>
#include <unordered_set>
#include "main.cpp" // Include the main file to access the functions and structures

using namespace std;

void print_channels(const vector<Edge>& edges) {
    for (const auto& edge : edges) {
        cout << "Edge " << edge.idx << ": ";
        for (const auto& channel : edge.channels) {
            cout << channel << " ";
        }
        cout << endl;
    }
}

void test_free_old_used_channels() {
    // Mock data
    input.N = 3;
    input.M = 2;
    input.edges = {
        {0, 0, 1, false, vector<int>(40, -1), {}},
        {1, 1, 2, false, vector<int>(40, -1), {}}
    };
    input.services = {
        {0, 0, 2, 2, 0, 1, 10, {{input.edges[0], 0, 1}, {input.edges[1], 0, 1}}, false},
        {1, 0, 2, 2, 2, 3, 20, {{input.edges[0], 2, 3}, {input.edges[1], 2, 3}}, false}
    };

    // Set initial channels
    for (int i = 0; i <= 1; ++i) {
        input.edges[0].channels[i] = 0;
        input.edges[1].channels[i] = 0;
    }
    for (int i = 2; i <= 3; ++i) {
        input.edges[0].channels[i] = 1;
        input.edges[1].channels[i] = 1;
    }

    // New paths with changed frequencies
    vector<vector<EdgeWithWavelengths>> new_paths = {
        {{input.edges[0], 6, 7}, {input.edges[1], 6, 7}},
        {{input.edges[0], 8, 9}, {input.edges[1], 8, 9}}
    };

    // Call the function
    Replanning::free_old_used_channels(new_paths, input.services);

    // Print results
    print_channels(input.edges);
}

int main() {
    test_free_old_used_channels();
    return 0;
}