/******************************************************************************
 * External_PQ.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Adil Chhabra <hg269@uni-heidelberg.de>
 *****************************************************************************/

#pragma once

#ifndef EXTERNAL_PQ_H
#define EXTERNAL_PQ_H

#include "definitions.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <stxxl/priority_queue>
#include <limits>

struct BlockNeighbor {
    PartitionID block;
    LongNodeID neighbor;

    // Overload the << operator to enable printing BlockNeighbor objects
    friend std::ostream &operator<<(std::ostream &os, const BlockNeighbor &bn) {
        os << "(" << bn.block << ", " << bn.neighbor << ")";
        return os;
    }
};

struct ComparatorGreater {
    bool operator()(const BlockNeighbor &a, const BlockNeighbor &b) const {
        return (a.neighbor > b.neighbor);
    }

    BlockNeighbor min_value() const {
        return {std::numeric_limits<PartitionID>::max(), std::numeric_limits<LongNodeID>::max()};
    }
};

class ExternalPQ {
public:
    ExternalPQ() = default;

    void appendToExtPQ(PartitionID block, LongNodeID target) {
        BlockNeighbor to_insert;
        to_insert.block = block;
        to_insert.neighbor = target;
        part_id_pq.push(to_insert);
    }

    PartitionID getBlock(LongNodeID source) const {
        // use time forward processing to fetch next part ID
        if (part_id_pq.empty()) {
            std::cout << "Attempted to pop from empty PQ." << std::endl;
            exit(1);
        }
        BlockNeighbor next = part_id_pq.top();
        if (source != next.neighbor) {
            std::cout << "Node mismatch in external PQ." << std::endl;
            std::cout << "Source = " << source << std::endl;
            std::cout << "PQ Entry: (" << next.block << ", " << next.neighbor << ")" << std::endl;
            exit(1);
        }
        part_id_pq.pop();
        return next.block;
    }

private:
    typedef stxxl::PRIORITY_QUEUE_GENERATOR<BlockNeighbor, ComparatorGreater,
            3 * 1048576, 1024 * 1024>::result pqueue_type;
    typedef pqueue_type::block_type block_type;

    //use pairs - block, neighbor for PQ with time forward processing
    const unsigned int mem_for_pools = 1024 * 1024; // restricts memory consumption of the pools
    stxxl::read_write_pool <block_type> pool{
            ((mem_for_pools / 2) / block_type::raw_size, (mem_for_pools / 2) / block_type::raw_size)};
    pqueue_type part_id_pq{(pool)}; // creates priority queue object with read-write-pool
};


#endif /* end of include guard: EXTERNAL_PQ_H */
