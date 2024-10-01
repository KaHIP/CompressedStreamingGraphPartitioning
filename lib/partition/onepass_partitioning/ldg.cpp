/******************************************************************************
 * ldg.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj
 *****************************************************************************/

#include "ldg.h"
#include "partition/onepass_partitioning/floating_block.h"
#include "partition/onepass_partitioning/vertex_partitioning.h"

onepass_ldg::onepass_ldg(PartitionID k0, PartitionID kf, PartitionID max_blocks,
                         NodeID n_threads, bool hashing /*=false*/)
    : vertex_partitioning(k0, kf, max_blocks, n_threads, hashing) {}

onepass_ldg::~onepass_ldg() {}

float onepass_ldg::compute_score(floating_block &block, int my_thread) {
  return block.get_ldg_obj(my_thread);
}

void onepass_ldg::instantiate_blocks(LongNodeID n, LongEdgeID m, PartitionID k,
                                     ImbalanceType epsilon) {
  vertex_partitioning::instantiate_blocks(n, m, k, epsilon);
}
