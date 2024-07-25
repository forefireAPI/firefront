#include "../src/ANN.h"
#include <iostream>
#include <vector>
#include <random>
#include <chrono>

int main() {
    Network network;
    network.loadFromFile("mbase.ffann");  // Ensure this path is correct

    std::cout << "Network Configuration:\n" << network.toString() << std::endl;

    // Constants for testing
    const size_t numInputs = 10000000;
    const size_t inputSize = network.layers.front().weights[0].size();

    // Initialize random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(0.0, 1.0);

    // Prepare inputs
    std::vector<std::vector<float>> inputs(numInputs, std::vector<float>(inputSize));
    for (auto& inp : inputs) {
        std::generate(inp.begin(), inp.end(), [&]() { return dis(gen); });
    }

    // Prepare to store results
    std::vector<float> results;
    results.reserve(numInputs);

    // Timing the feedforward process
    auto start = std::chrono::high_resolution_clock::now();

    // Process each input through the network
    for (const auto& input : inputs) {
        results.push_back(network.processInput(input)[0]);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    // Output timing results
    std::cout << "Time taken for processing " << numInputs << " inputs: " << elapsed.count() << " seconds\n";

    // Specific inputs for verification
    std::vector<std::vector<float>> testInputs = {
        {1, 1, 1, 0, 0},
        {1, 0, 0, 0, 0},
        {0, 1, 0, 0, 0},
        {0, 0, 1, 0, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 0, 0, 1},
        {1, 1, 0, 1, 1}
    };

    // Process each test input through the network and print the results
    for (const auto& input : testInputs) {
        auto outputs = network.processInput(input);
        std::cout << "[";
        for (auto i : input) std::cout << i << " ";
        std::cout << "] result: " << outputs[0] << std::endl;
    }

    return 0;
}
