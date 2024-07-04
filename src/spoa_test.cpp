#include <iostream>
#include <chrono>
#include <fstream> 
#include <iostream> 
#include <string> 
#include "spoa/spoa.hpp"

int main(int argc, char** argv) {
    // Check if the filename is provided as a command-line argument
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1; // Return an error code
    }

    std::string filename = argv[1]; // Get the filename from the command-line argument
    std::ifstream file(filename); // Open the file for reading

    // Check if the file was successfully opened
    if (!file.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return 1; // Return an error code
    }
    std::vector<std::string> sequences;
    std::string line; // Temporary string to hold each line
    int line_index = 0;
    // Read the file line by line
    while (std::getline(file, line)) {
        if line_index % 2 == 1 {
            sequences.push_back(line); // Add each line to the vector
        }
        line_index++;
    }
    file.close(); // Close the file

    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 2, -2, -2);  // linear gaps
    spoa::Graph graph{};

    int run_index = 0;
    for (const auto& it : sequences) {
        auto alignment = alignment_engine->Align(it, graph);
        graph.AddAlignment(alignment, it);
        if run_index == 1 { //start time measurement
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        }
        run_index++;
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    auto consensus = graph.GenerateConsensus();
    std::cerr << ">Consensus LN:i:" << consensus.size() << std::endl << consensus << std::endl;
    return 0;
}