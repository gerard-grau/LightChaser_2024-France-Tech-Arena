#include "Input.cpp"
#include "Bottleneck.cpp"
#include "Replanning.cpp"

Input input;

int main() {
    read_input();
    
    Bottleneck();
    
    Replanning();

    return 0;
}
