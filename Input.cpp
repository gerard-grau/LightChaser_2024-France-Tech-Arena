#include "Utils.cpp"
#include <iostream>

extern Input input;

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
