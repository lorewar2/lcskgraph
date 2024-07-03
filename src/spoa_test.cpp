#include <iostream>
#include <chrono>
#include <fstream> 
#include <iostream> 
#include <string> 

#include "spoa/spoa.hpp"

int main(int argc, char** argv) {
    std::vector<std::string> sequences;
    sequences.push_back("CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT");
    sequences.push_back("ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT");
    sequences.push_back("ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG");
    sequences.push_back("CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA");
    sequences.push_back("GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT");
    sequences.push_back("CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT");
    ifstream file("Synthetic10-0.fa");
    auto alignment_engine = spoa::AlignmentEngine::Create(
        spoa::AlignmentType::kNW, 2, -2, -2);  // linear gaps

    spoa::Graph graph{};
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (const auto& it : sequences) {
    auto alignment = alignment_engine->Align(it, graph);
        graph.AddAlignment(alignment, it);
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    auto consensus = graph.GenerateConsensus();

    std::cerr << ">Consensus LN:i:" << consensus.size() << std::endl << consensus << std::endl;

    auto msa = graph.GenerateMultipleSequenceAlignment();

    for (const auto& it : msa) {
    std::cerr << it << std::endl;
    }

    return 0;
}