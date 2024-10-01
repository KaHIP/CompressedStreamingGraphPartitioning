/******************************************************************************
 * vertex_partitioning.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj
 *****************************************************************************/

#include "vertex_partitioning.h"
#include "partition/onepass_partitioning/floating_block.h"

vertex_partitioning::vertex_partitioning(PartitionID k0, PartitionID kf,
                                         PartitionID max_blocks,
                                         NodeID n_threads,
                                         bool hashing /*=false*/) {
  this->k0 = k0;
  this->kf = kf;
  this->max_blocks = max_blocks;
  this->n_threads = n_threads;
  this->neighbor_blocks.resize(n_threads);
  original_problem = NULL;
  subproblem_tree_root = NULL;
  parent_block = NULL;
  this->hashing = hashing;
  this->use_self_sorting_array = false;
  int largest_dim = max_blocks;

  blocks.reserve(largest_dim); // This is very important to ensure that the
                               // memory address of each block
  // will never change in the subproblems. When the memory addresses change,
  // random errors occurr.
}

void vertex_partitioning::instantiate_blocks(LongNodeID n, LongEdgeID m,
                                             PartitionID k,
                                             ImbalanceType epsilon) {
  if (blocks.size() > 0)
    return;
  base_size_constraint = ceil(((100 + epsilon) / 100.) * (n / (double)k));
  for (PartitionID id = 0; id < k; id++) {
    blocks.push_back(floating_block(this, id, n_threads));
    blocks[id].set_capacity(base_size_constraint);
  }
  if (this->use_self_sorting_array == true) {
    this->sorted_blocks.initialize(k, (NodeWeight)0);
  }
}

float vertex_partitioning::compute_score(floating_block &block, int my_thread) {
  return 0;
}
