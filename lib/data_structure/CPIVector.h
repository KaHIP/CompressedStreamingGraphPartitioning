/******************************************************************************
 * CPIVector.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Adil Chhabra <hg269@uni-heidelberg.de>
 *****************************************************************************/
#pragma once

#ifndef CPI_VECTOR_H
#define CPI_VECTOR_H

#include "definitions.h"
#include <vector>
#include <iostream>
#include <fstream>

namespace {

#include "cpi/run_length_compression.hpp"

}
class CPIVector {
public:
    CPIVector() = default;

    void appendNodeToEncodedVector(NodeID newNodeValue);
    NodeID findValueAtIndex(NodeID index) const;

private:
    cpi::RunLengthCompression<PartitionID> encodedVector;
};


#endif /* end of include guard: CPI_VECTOR_H */
