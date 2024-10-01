/******************************************************************************
 * RunLengthEncodedVector.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Adil Chhabra <hg269@uni-heidelberg.de>
 *****************************************************************************/

#ifndef RUN_LENGTH_ENCODED_VECTOR_H
#define RUN_LENGTH_ENCODED_VECTOR_H

#include "definitions.h"
#include <vector>
#include <iostream>
#include <fstream>

class RunLengthEncodedVector {
public:
    RunLengthEncodedVector() = default;

    void runLengthEncoding(const std::vector<NodeID>& input) {
        if (input.empty()) {
            return;
        }

        NodeID currentValue = input[0];
        NodeID run_count = 1;
        NodeID cumulative_count = 1;

        for (size_t i = 1; i < input.size(); ++i) {
            if (input[i] == currentValue) {
                run_count++;
            } else {
                encodedVector.push_back({currentValue, run_count, cumulative_count});
                currentValue = input[i];
                run_count = 1;
            }
            cumulative_count++;
        }
        encodedVector.push_back({currentValue, run_count, cumulative_count});
    }

    void reserveEncodedVectorSize(NodeID reservation_size) {
        encodedVector.reserve(reservation_size);
    }

    void appendNodeToEncodedVector(NodeID newNodeValue) {
        if (!encodedVector.empty()) {
            Run& lastRun = encodedVector.back();
            if (lastRun.value == newNodeValue) {
                lastRun.count++;
                lastRun.cumulativeCount++;
            } else {
                this->encodedVector.push_back({newNodeValue, 1, lastRun.cumulativeCount+1});
            }
        } else {
            this->encodedVector.push_back({newNodeValue, 1,1});
        }
    }

    std::vector<NodeID> runLengthDecoding() {
        std::vector<NodeID> result;

        for (const Run& run : encodedVector) {
            result.insert(result.end(), run.count, run.value);
        }

        return result;
    }

    void displayEncodedVector() const {
        for (const Run& run : encodedVector) {
            std::cout << "(" << run.value << ", " << run.count << ", " << run.cumulativeCount << ") ";
        }
        std::cout << std::endl;
    }

    void writeEncodedVectorToFile(const std::string& fileName) const {
        std::ofstream outputFile(fileName);
        if (!outputFile.is_open()) {
            std::cerr << "Error: Unable to open file " << fileName << " for writing." << std::endl;
            return;
        }

        for (const Run& run : encodedVector) {
            outputFile << "(" << run.value << ", " << run.count << ", " << run.cumulativeCount << ") ";
        }
        outputFile << std::endl;

        outputFile.close();
    }

    NodeID findValueAtIndex(NodeID index) const {
        if (encodedVector.empty()) {
            return -1;
        }

        if(index == 0) {
            return 1;
        }

        size_t left = 0;
        size_t right = encodedVector.size() - 1;

        while (left <= right) {
            size_t mid = left + (right - left) / 2;
            const Run& currentRun = encodedVector[mid];
            //std::cout << "At mid: " << mid << " : " << currentRun.value << " and " << currentRun.cumulativeCount << " and " << currentRun.count << std::endl;
            if (index <= currentRun.cumulativeCount - currentRun.count) {
                // The index is in the current run
                right = mid - 1;
                //return currentRun.value;
            } else if (index > currentRun.cumulativeCount) {
                // Search in the right half
                left = mid + 1;
            } else {
                //The index is within the current count range
                return currentRun.value;
            }
        }
        return -1; //if index out of range
        //throw std::out_of_range("Index out of range");
    }

    NodeID memoryConsumption() const {
        return (sizeof(Run) * encodedVector.size()) / 1000;
    }

    NodeID getEncodedVectorSize() {
        return encodedVector.size();
    }

private:
    struct Run {
        NodeID value;
        NodeID count;
        NodeID cumulativeCount;

        Run(NodeID v, NodeID c, NodeID cc) : value(v), count(c), cumulativeCount(cc) {}
    };

    std::vector<Run> encodedVector;
};


#endif /* end of include guard: RUN_LENGTH_ENCODED_VECTOR_H */
