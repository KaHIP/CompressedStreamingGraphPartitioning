/******************************************************************************
 * vertex_partitioning.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj
 *****************************************************************************/

#ifndef VERTEX_PARTITIONING_7I4IR31Y
#define VERTEX_PARTITIONING_7I4IR31Y

#include "definitions.h"
#include "random_functions.h"
#include "timer.h"
#include <algorithm>
#include <omp.h>

#include "data_structure/priority_queues/self_sorting_monotonic_vector.h"
#include "partition/onepass_partitioning/floating_block.h"

#define MIN(A, B) (((A) < (B)) ? (A) : (B))
#define MAX(A, B) (((A) > (B)) ? (A) : (B))

class vertex_partitioning {
public:
  vertex_partitioning(PartitionID k0, PartitionID kf, PartitionID max_blocks,
                      NodeID n_threads, bool hashing = false);
  ~vertex_partitioning() {}
  virtual void instantiate_blocks(LongNodeID n, LongEdgeID m, PartitionID k,
                                  ImbalanceType epsilon);
  void load_edge(PartitionID block, EdgeWeight e_weight, int my_thread);
  void clear_edgeweight(PartitionID block, int my_thread);
  void propagate_load_edge(EdgeWeight e_weight, PartitionID block,
                           int my_thread);
  void propagate_clear_edgeweight(int my_thread);
  void remove_nodeweight(PartitionID block, NodeWeight n_weight);
  void remove_parent_nodeweight(NodeWeight n_weight);
  PartitionID solve_node(LongNodeID curr_node_id, NodeWeight curr_node_weight,
                         PartitionID previous_assignment, double kappa,
                         int my_thread);
  PartitionID solve_hashing(LongNodeID curr_node_id,
                            NodeWeight curr_node_weight, int my_thread);
  PartitionID set_decision(PartitionID block, LongNodeID curr_node_id,
                           NodeWeight curr_node_weight, int my_thread);
  void set_original_problem(vertex_partitioning *original_problem);
  void set_parent_block(floating_block *parent_block);
  void clear_edgeweight_blocks(
      std::vector<std::pair<PartitionID, EdgeWeight>> &neighbor_blocks_thread,
      PartitionID n_elements, int my_thread);

  // methods for dealing with blocks
  void load_edge_block(floating_block &real_block, EdgeWeight e_weight,
                       int my_thread);
  void clear_edgeweight_block(floating_block &real_block, int my_thread);
  void remove_nodeweight_block(floating_block &real_block, NodeWeight n_weight);
  PartitionID assign_node_block(floating_block &real_block,
                                LongNodeID curr_node_id,
                                NodeWeight curr_node_weight, int my_thread);
  PartitionID force_decision_block(floating_block &real_block,
                                   LongNodeID curr_node_id,
                                   NodeWeight curr_node_weight, int my_thread);
  floating_block *get_block_address(PartitionID block_id) const;
  void enable_self_sorting_array();

  std::vector<floating_block> blocks;

protected:
  PartitionID solve(LongNodeID curr_node_id, NodeWeight curr_node_weight,
                    PartitionID previous_assignment, double kappa,
                    int my_thread);
  PartitionID solve_linear_complexity(LongNodeID curr_node_id,
                                      NodeWeight curr_node_weight,
                                      int my_thread);
  virtual float compute_score(floating_block &block, int my_thread);

  void gothrough_neighborhood(PartitionID &decision, float &best, int my_thread);
  void check_best_nonneighbor(PartitionID &decision, float &best, int my_thread);

  NodeID n_threads;
  PartitionID k0;
  PartitionID kf;
  PartitionID max_blocks;
  ImbalanceType base_size_constraint;
  bool hashing;
  std::vector<std::vector<PartitionID>> neighbor_blocks;
  random_functions::fastRandBool<uint64_t> random_obj;
  self_sorting_monotonic_vector<PartitionID, NodeWeight> sorted_blocks;
  bool use_self_sorting_array;

  vertex_partitioning *original_problem;
  vertex_partitioning *subproblem_tree_root;
  floating_block *parent_block;
};

inline void vertex_partitioning::clear_edgeweight_blocks(
    std::vector<std::pair<PartitionID, EdgeWeight>> &neighbor_blocks_thread,
    PartitionID n_elements, int my_thread) {
  for (PartitionID key = 0; key < n_elements; key++) {
    auto &element = neighbor_blocks_thread[key];
    clear_edgeweight(element.first, my_thread);
  }
}

inline void vertex_partitioning::load_edge(PartitionID block,
                                           EdgeWeight e_weight, int my_thread) {
  load_edge_block(blocks[block], e_weight, my_thread);
}

inline void vertex_partitioning::load_edge_block(floating_block &real_block,
                                                 EdgeWeight e_weight,
                                                 int my_thread) {
  real_block.increment_edgeweight(my_thread, e_weight);
  real_block.parent_problem->propagate_load_edge(
      e_weight, real_block.get_block_id(), my_thread);
}

inline void vertex_partitioning::propagate_load_edge(EdgeWeight e_weight,
                                                     PartitionID block,
                                                     int my_thread) {
  neighbor_blocks[my_thread].push_back(block);
}

inline void vertex_partitioning::clear_edgeweight(PartitionID block,
                                                  int my_thread) {
  clear_edgeweight_block(blocks[block], my_thread);
}

inline void
vertex_partitioning::clear_edgeweight_block(floating_block &real_block,
                                            int my_thread) {
  real_block.set_edgeweight(my_thread, 0);
  real_block.parent_problem->propagate_clear_edgeweight(my_thread);
}

inline void vertex_partitioning::propagate_clear_edgeweight(int my_thread) {
  neighbor_blocks[my_thread].clear();
}

inline void vertex_partitioning::remove_nodeweight(PartitionID block,
                                                   NodeWeight n_weight) {
  remove_nodeweight_block(blocks[block], n_weight);
}

inline void
vertex_partitioning::remove_nodeweight_block(floating_block &real_block,
                                             NodeWeight n_weight) {
  real_block.increment_curr_weight(-n_weight);
}

inline void vertex_partitioning::remove_parent_nodeweight(NodeWeight n_weight) {
  if (parent_block != NULL) {
    remove_nodeweight_block(*parent_block, n_weight);
  }
}

inline PartitionID
vertex_partitioning::set_decision(PartitionID block, LongNodeID curr_node_id,
                                  NodeWeight curr_node_weight, int my_thread) {
  return assign_node_block(blocks[block], curr_node_id, curr_node_weight,
                           my_thread);
}

inline PartitionID vertex_partitioning::assign_node_block(
    floating_block &real_block, LongNodeID curr_node_id,
    NodeWeight curr_node_weight, int my_thread) {
  PartitionID decision = real_block.get_block_id();
  real_block.increment_curr_weight(curr_node_weight);
  return decision;
}

inline PartitionID vertex_partitioning::force_decision_block(
    floating_block &real_block, LongNodeID curr_node_id,
    NodeWeight curr_node_weight, int my_thread) {
  return real_block.parent_problem->set_decision(
      real_block.get_block_id(), curr_node_id, curr_node_weight, my_thread);
}

inline PartitionID vertex_partitioning::solve_node(
    LongNodeID curr_node_id, NodeWeight curr_node_weight,
    PartitionID previous_assignment, double kappa, int my_thread) {
  PartitionID decision;
  if (hashing) {
    decision = solve_hashing(curr_node_id, curr_node_weight, my_thread);
  } else {
    decision = solve(curr_node_id, curr_node_weight, previous_assignment, kappa,
                     my_thread);
  }

  return decision;
}

inline PartitionID
vertex_partitioning::solve_hashing(LongNodeID curr_node_id,
                                   NodeWeight curr_node_weight, int my_thread) {
  PartitionID k = blocks.size();
  PartitionID random_add = random_functions::nextIntHashing(k);
  /* PartitionID random_add = random_functions::nextInt(0, k-1); */
  PartitionID decision = crc32(curr_node_id + random_add) % k;
  for (PartitionID i = 0; blocks[decision].fully_loaded() && i < 100; i++) {
    decision = crc32(curr_node_id + random_add + i * k) % (blocks.size());
  }
  return set_decision(decision, curr_node_id, curr_node_weight, my_thread);
}

inline PartitionID vertex_partitioning::solve(LongNodeID curr_node_id,
                                              NodeWeight curr_node_weight,
                                              PartitionID previous_assignment,
                                              double kappa, int my_thread) {
  float best = std::numeric_limits<float>::lowest();
  float score;
  PartitionID k = blocks.size();
  PartitionID decision = random_functions::nextIntHashing(k);
  /* PartitionID decision = random_functions::nextInt(0, k-1); */
  for (auto &block : blocks) {
    if (block.fully_loaded()) {
      continue;
    }
    score = compute_score(block, my_thread);
    // float old_score = score;
    if (block.get_block_id() == previous_assignment) {
      if (score > 0) {
        score = score * kappa; // mult not good enough
      } else {
        score = score / kappa;
      }
    }
    if (score > best || (random_obj.nextBool() && score == best)) {
      decision = block.get_block_id();
      best = score;
    }
  }
  return set_decision(decision, curr_node_id, curr_node_weight, my_thread);
}

inline PartitionID vertex_partitioning::solve_linear_complexity(
    LongNodeID curr_node_id, NodeWeight curr_node_weight, int my_thread) {
  float best = std::numeric_limits<float>::lowest();
  PartitionID decision = random_functions::nextIntHashing(blocks.size());
  /* NodeID decision = crc32(curr_node_id)% (blocks.size()); */
  gothrough_neighborhood(decision, best, my_thread);
  check_best_nonneighbor(decision, best, my_thread);
  /* amortized_sampling_for_feasibility(decision); */
  return set_decision(decision, curr_node_id, curr_node_weight, my_thread);
}

inline void vertex_partitioning::gothrough_neighborhood(PartitionID &decision,
                                                        float &best,
                                                        int my_thread) {
  float score;
  for (auto id : neighbor_blocks[my_thread]) {
    auto &block = blocks[id];
    if (block.fully_loaded()) {
      continue;
    }
    score = compute_score(block, my_thread);
    if (score > best) {
      decision = block.get_block_id();
      best = score;
    }
  }
}

inline void vertex_partitioning::check_best_nonneighbor(PartitionID &decision,
                                                        float &best,
                                                        int my_thread) {
  float score;
  PartitionID id = this->sorted_blocks[0]; // it does not matter whether or not
                                           // it is a nonneighbor
  auto &block = blocks[id];
  score = compute_score(block, my_thread);
  if (score > best) {
    decision = block.get_block_id();
    best = score;
  }
}

inline void
vertex_partitioning::set_parent_block(floating_block *parent_block) {
  this->parent_block = parent_block;
}

inline void vertex_partitioning::set_original_problem(
    vertex_partitioning *original_problem) {
  this->original_problem = original_problem;
}

inline floating_block *
vertex_partitioning::get_block_address(PartitionID block_id) const {
  floating_block *block_pointer = &(blocks[block_id]);
  return block_pointer;
}

inline void vertex_partitioning::enable_self_sorting_array() {
  this->use_self_sorting_array = true;
}

#endif /* end of include guard: VERTEX_PARTITIONING_7I4IR31Y */
