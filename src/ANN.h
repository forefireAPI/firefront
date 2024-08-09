/*

Copyright (C) 2024 ForeFire Team, SPE, Universit� de Corse.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 US

 */

#ifndef ANN_H
#define ANN_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <cstring>
#include <sstream>
#include <stdexcept>

// Activation functions
inline float sigmoid(float x) {
    
    return 1.0f / (1.0f + std::exp(-x));
}

inline float relu(float x) {
    return std::max(0.0f, x);
}

inline float linear(float x) {
    return x;
}

// Activation function lookup
inline float (*getActivationFunction(const std::string& name))(float) {
    if (name == "RELU") return relu;
    if (name == "SIGM") return sigmoid;
    if (name == "LINE") return linear;
    return nullptr;
}

// Layer structure
struct Layer {
    std::vector<float> neurons;
    std::vector<std::vector<float>> weights;
    std::vector<float> biases;
    float (*activation)(float);

    // Constructor
    Layer(int inputSize, int outputSize, float (*actFunc)(float))
        : neurons(outputSize, 0.0f), weights(outputSize, std::vector<float>(inputSize)), biases(outputSize), activation(actFunc) {}



    void loadWeightsAndBiases(const std::vector<float>& weightData, const std::vector<float>& biasData) {
        for (size_t i = 0; i < weights.size(); ++i) {
            for (size_t j = 0; j < weights[0].size(); ++j) {
                weights[i][j] = weightData[j * weights.size() + i];
            }
        }
        std::copy(biasData.begin(), biasData.end(), biases.begin());
    }

    std::vector<float> feedforward(const std::vector<float>& inputs) {
        std::vector<float> output(neurons.size(), 0.0f);
        for (size_t i = 0; i < neurons.size(); ++i) {
            float neuronOutput = biases[i];
            for (size_t j = 0; j < inputs.size(); ++j) {
                neuronOutput += weights[i][j] * inputs[j];
            }
            output[i] = activation(neuronOutput);
        }
        return output;
    }
};

// Network structure
struct Network {
    std::vector<Layer> layers;
    std::vector<std::string> inputNames;
    std::vector<std::string> outputNames;

    std::vector<std::string> splitNames(const std::string& names) {
        std::vector<std::string> result;
        std::istringstream iss(names);
        std::string token;
        while (getline(iss, token, ',')) {
            result.push_back(token);
        }
        return result;
    }

    void loadFromFile(const char* filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open network structure file.");
        }

        char header[9] = {0};
        int numLayers;
        file.read(header, 8);
        file.read(reinterpret_cast<char*>(&numLayers), sizeof(int));

        if (std::string(header) != "FFANN001") {
            throw std::runtime_error("Invalid file format.");
        }

        int inputSize = 0;  // This will be set based on the first layer's weight matrix dimensions

        for (int i = 0; i < numLayers; ++i) {
            char activation[5] = {0};
            int width, height;
            file.read(activation, 4);
            file.read(reinterpret_cast<char*>(&width), sizeof(int));
            file.read(reinterpret_cast<char*>(&height), sizeof(int));

            if (i == 0) {
                inputSize = width;
            }

            float (*actFunc)(float) = getActivationFunction(std::string(activation));
            if (!actFunc) {
                throw std::runtime_error("Unsupported activation function.");
            }

            Layer layer(inputSize, height, actFunc);
            std::vector<float> weightData(height * width);
            std::vector<float> biasData(height);
            file.read(reinterpret_cast<char*>(weightData.data()), height * width * sizeof(float));
            file.read(reinterpret_cast<char*>(biasData.data()), height * sizeof(float));

            layer.loadWeightsAndBiases(weightData, biasData);
            layers.push_back(layer);
            inputSize = height;
        }

    int input_names_length, output_names_length;
    std::string nameBuffer;
    file.read(reinterpret_cast<char*>(&input_names_length), sizeof(int));
    nameBuffer.resize(input_names_length);
    file.read(&nameBuffer[0], input_names_length);
    inputNames = splitNames(nameBuffer);

    file.read(reinterpret_cast<char*>(&output_names_length), sizeof(int));
    nameBuffer.resize(output_names_length);
    file.read(&nameBuffer[0], output_names_length);
    outputNames = splitNames(nameBuffer);

    std::cout << "name sizes : " << input_names_length <<  "  " << output_names_length << std::endl;
        printLayerInfo();    
    }
    void printLayerInfo() {
        for (const auto& layer : layers) {
            std::cout << "Layer: Input Size = " << layer.weights[0].size() 
                    << ", Output Size = " << layer.neurons.size()
                    << ", Activation Function = ";
            if (layer.activation == relu) {
                std::cout << "ReLU";
            } else if (layer.activation == sigmoid) {
                std::cout << "Sigmoid";
            } else {
                std::cout << "Linear";
            }
            std::cout << std::endl;  // End the line after printing each layer's details
        }
    }
    std::vector<float> processInput(const std::vector<float>& input) {
        std::vector<float> result = input;
        for (auto& layer : layers) {
            result = layer.feedforward(result);
        }
        return result;
    }
    // Assuming this method processes only one input set and expects exactly one output
    bool processDirectInput(const std::vector<double>& inputs, double& output) {
        std::vector<float> currentOutput(inputs.begin(), inputs.end());

        for (auto& layer : layers) {
            std::vector<float> nextOutput(layer.neurons.size(), 0.0f);

            for (size_t i = 0; i < layer.neurons.size(); ++i) {
                float neuronOutput = layer.biases[i];
                for (size_t j = 0; j < currentOutput.size(); ++j) {
                    neuronOutput += layer.weights[i][j] * currentOutput[j];
                }
                nextOutput[i] = layer.activation(neuronOutput);
            }

            currentOutput = std::move(nextOutput);
        }

        // Assumes the final layer has exactly one output
        if (currentOutput.size() == 1) {
            output = static_cast<double>(currentOutput[0]);
            return true;
        }
        return false;
    }
    std::string toString() const {
        std::ostringstream ss;
        for (const auto& layer : layers) {
            ss << "Layer: Input Size = " << layer.weights[0].size() << ", Output Size = " << layer.neurons.size()
               << ", Activation Function = " << (layer.activation == relu ? "ReLU" :
                                                layer.activation == sigmoid ? "Sigmoid" : "Linear") << "\n";
        }
        return ss.str();
    }
};

#endif // ANN_H

