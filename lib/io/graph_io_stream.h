/******************************************************************************
 * graph_io_stream.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj
 *****************************************************************************/

#ifndef GRAPHIOSTREAM_H_
#define GRAPHIOSTREAM_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <ostream>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <unordered_map>
#include <vector>

#include "cpi/run_length_compression.hpp"
#include "data_structure/ExternalPQ.h"
#include "data_structure/buffered_map.h"
#include "data_structure/graph_access.h"
#include "data_structure/single_adj_list.h"
#include "definitions.h"
#include "partition/onepass_partitioning/fennel.h"
#include "partition/onepass_partitioning/vertex_partitioning.h"
#include "partition/partition_config.h"
#include "random_functions.h"
#include "timer.h"
#include "data_structure/compression_vectors/CompressionDataStructure.h"
#include "data_structure/compression_vectors/RunLengthCompressionVector.h"
#include "data_structure/compression_vectors/BatchRunLengthCompression.h"
//#include <kagen.h>

#define MIN(A, B) (((A) < (B)) ? (A) : (B))
#define MAX(A, B) (((A) > (B)) ? (A) : (B))

typedef std::vector<std::string> *LINE_BUFFER;

class graph_io_stream {
public:
  graph_io_stream();

  virtual ~graph_io_stream();

  static void readFirstLineStream(PartitionConfig &partition_config,
                                  std::string graph_filename,
                                  EdgeWeight &total_edge_cut, EdgeWeight &qap);

  static void
  loadRemainingLinesToBinary(PartitionConfig &partition_config,
                             std::vector<std::vector<LongNodeID>> *&input);

  static void
  loadBufferLinesToBinary(PartitionConfig &partition_config,
                          std::vector<std::vector<LongNodeID>> *&input,
                          LongNodeID num_lines);

  static std::vector<std::vector<LongNodeID>> *
  loadLinesFromStreamToBinary(PartitionConfig &partition_config,
                              LongNodeID num_lines);

  static void
  generateNeighborhood(PartitionConfig &partition_config,
                       std::vector<NodeID> &generated_neighbors,
                       std::vector<std::vector<LongNodeID>> *&input);

  static void readNodeOnePass(PartitionConfig &config, LongNodeID curr_node,
                              int my_thread,
                              std::vector<std::vector<LongNodeID>> *&input,
                              const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments,
                              vertex_partitioning *onepass_partitioner);

  static void streamEvaluatePartition(PartitionConfig &config,
                                      const std::string &filename,
                                      EdgeWeight &edgeCut, EdgeWeight &qap,
                                      const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments);

  static void writePartitionStream(PartitionConfig &config,
                                   const std::string &filename,
                                   const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments);

  static void readPartition(PartitionConfig &config,
                            const std::string &filename);

  static void configureGeneratedStream(PartitionConfig &partition_config,
                                                   std::string graph_filename,
                                                   uint64_t num_edges, EdgeWeight &total_edge_cut,EdgeWeight &qap);

//    static void streamEvaluatePartitionGenerated(PartitionConfig &config,
//                                        const std::string &filename,
//                                        EdgeWeight &edgeCut, EdgeWeight &qap,
//                                        kagen::StreamingGenerator &streamGenerator,
//                                        const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments);
};

inline void
graph_io_stream::readNodeOnePass(PartitionConfig &config, LongNodeID curr_node,
                                 int my_thread,
                                 std::vector<std::vector<LongNodeID>> *&input,
                                 const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments,
                                 vertex_partitioning *onepass_partitioner) {
  /* NodeWeight total_nodeweight = 0; */
  auto &read_ew = config.read_ew;
  auto &read_nw = config.read_nw;
  LongNodeID target;
  NodeWeight weight;
  /* LongNodeID nmbNodes = 1; */

  LongNodeID cursor = (config.ram_stream) ? curr_node : 0;

  if ((config.one_pass_algorithm == ONEPASS_HASHING) ||
      (config.one_pass_algorithm == ONEPASS_HASHING_CRC32)) {
    return;
  }

  auto &all_blocks_to_keys = config.all_blocks_to_keys[my_thread];
  auto &next_key = config.next_key[my_thread];
  auto &neighbor_blocks = config.neighbor_blocks[my_thread];

  onepass_partitioner->clear_edgeweight_blocks(neighbor_blocks, next_key,
                                               my_thread);
  next_key = 0;

  std::vector<LongNodeID> &line_numbers = (*input)[cursor];
  LongNodeID col_counter = 0;
  weight = (read_nw) ? line_numbers[col_counter++] : 1;
  /* total_nodeweight += weight; */

#ifdef MODE_STREAMMULTISECTION
  if (config.restream_number) {
    NodeID node = 0;
    LongNodeID lower_global_node =
        curr_node + 1; // Bounds below start from 1 instead of 0
    LongNodeID global_node = lower_global_node + (LongNodeID)node - 1;
    PartitionID nodeGlobalPar = (*config.stream_nodes_assign)[global_node];
    (*config.stream_blocks_weight)[nodeGlobalPar] -= weight;
    onepass_partitioner->remove_nodeweight(nodeGlobalPar, weight);
  }
#endif

  PartitionID selecting_factor = (1 + (PartitionID)read_ew);
  config.edges = (line_numbers.size() - col_counter) / selecting_factor;
  float scaling_factor = 1;
  LongNodeID deg = 0;
  while (col_counter < line_numbers.size()) {
    target = line_numbers[col_counter++];
    deg++;
    EdgeWeight edge_weight = (read_ew) ? line_numbers[col_counter++] : 1;

    PartitionID targetGlobalPar, targetGlobalParTest;
    if (config.rle_length == -1) {
      targetGlobalPar = (*config.stream_nodes_assign)[target - 1];
    } else if (config.rle_length == 0) {
      if ((target - 1) < curr_node) {
          targetGlobalPar = block_assignments->GetValueByIndex(target-1);
      } else {
        targetGlobalPar = INVALID_PARTITION;
      }
    } else if (config.rle_length == -2) {
      if ((target - 1) < curr_node) {
        targetGlobalPar =
            (*config.external_pq_partition_assign).getBlock(curr_node);
      } else {
        targetGlobalPar = INVALID_PARTITION;
      }
    } else {
      if ((target - 1) < curr_node) {
        targetGlobalPar = block_assignments->GetValueByBatchIndex((target - 1) / config.rle_length, (target - 1) % config.rle_length);
      } else {
        targetGlobalPar = INVALID_PARTITION;
      }
    }
    if (targetGlobalPar != INVALID_PARTITION) {
      PartitionID key = all_blocks_to_keys[targetGlobalPar];
      if (key >= next_key || neighbor_blocks[key].first != targetGlobalPar) {
        all_blocks_to_keys[targetGlobalPar] = next_key;
        auto &new_element = neighbor_blocks[next_key];
        new_element.first = targetGlobalPar;
        new_element.second = edge_weight;
        next_key++;
      } else {
        neighbor_blocks[key].second += edge_weight;
      }
    }
  }

  if (deg > config.max_degree)
    config.max_degree = deg;

  for (PartitionID key = 0; key < next_key; key++) {
    auto &element = neighbor_blocks[key];
    onepass_partitioner->load_edge(element.first,
                                   element.second * scaling_factor, my_thread);
  }

  config.remaining_stream_nodes--;
}

inline void graph_io_stream::loadRemainingLinesToBinary(
    PartitionConfig &partition_config,
    std::vector<std::vector<LongNodeID>> *&input) {
  if (partition_config.ram_stream) {
    input = graph_io_stream::loadLinesFromStreamToBinary(
        partition_config, partition_config.remaining_stream_nodes);
  }
}

inline void graph_io_stream::loadBufferLinesToBinary(
    PartitionConfig &partition_config,
    std::vector<std::vector<LongNodeID>> *&input, LongNodeID num_lines) {
  if (!partition_config.ram_stream) {
    input = graph_io_stream::loadLinesFromStreamToBinary(partition_config,
                                                         num_lines);
  }
}

inline std::vector<std::vector<LongNodeID>> *
graph_io_stream::loadLinesFromStreamToBinary(PartitionConfig &partition_config,
                                             LongNodeID num_lines) {
  std::vector<std::vector<LongNodeID>> *input;
  input = new std::vector<std::vector<LongNodeID>>(num_lines);
  std::vector<std::string> *lines;
  lines = new std::vector<std::string>(1);
  LongNodeID node_counter = 0;
  buffered_input *ss2 = NULL;
  while (node_counter < num_lines) {
    std::getline(*(partition_config.stream_in), (*lines)[0]);
    if ((*lines)[0][0] == '%') { // a comment in the file
      continue;
    }
    ss2 = new buffered_input(lines);
    ss2->simple_scan_line((*input)[node_counter++]);
    (*lines)[0].clear();
    delete ss2;
  }
  delete lines;
  return input;
}

inline void
graph_io_stream::generateNeighborhood(PartitionConfig &partition_config,
                                      std::vector<NodeID> &generated_neighbors,
                                      std::vector<std::vector<LongNodeID>> *&input) {
    input = new std::vector<std::vector<LongNodeID>>(1);
    for(auto & neighbor: generated_neighbors) {
        //std::cout << neighbor << ", ";
        (*input)[0].push_back(neighbor);
    }
}

#endif /*GRAPHIOSTREAM_H_*/
