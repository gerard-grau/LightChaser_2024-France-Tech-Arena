#include "Utils.cpp"
#include <iostream>

using namespace std;

extern Input input;
vector<vector<int>> failures;

void Bottleneck() {
    
    print_bottleneck_results();
}

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

}