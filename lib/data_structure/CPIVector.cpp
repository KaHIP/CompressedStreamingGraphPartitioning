#include "CPIVector.h"
#include "cpi/run_length_compression.hpp"

void CPIVector::appendNodeToEncodedVector(NodeID newNodeValue) {
    encodedVector.push_back(newNodeValue);
}

NodeID CPIVector::findValueAtIndex(NodeID index) const {
    return encodedVector[index];
}