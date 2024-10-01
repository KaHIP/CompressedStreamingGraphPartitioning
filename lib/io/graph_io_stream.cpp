/******************************************************************************
 * graph_io_stream.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj
 *****************************************************************************/

#include "graph_io_stream.h"
#include "cpi/run_length_compression.hpp"
#include "timer.h"
#include <math.h>
#include <sstream>
//#include <kagen.h>
//#include <mpi.h>

#define MIN(A, B) (((A) < (B)) ? (A) : (B))
#define MAX(A, B) (((A) > (B)) ? (A) : (B))

graph_io_stream::graph_io_stream() {}

graph_io_stream::~graph_io_stream() {}

void graph_io_stream::readFirstLineStream(PartitionConfig &partition_config,
                                          std::string graph_filename,
                                          EdgeWeight &total_edge_cut,
                                          EdgeWeight &qap) {
  if (partition_config.stream_in != NULL) {
    delete partition_config.stream_in;
  }
  partition_config.stream_in = new std::ifstream(graph_filename.c_str());
  if (!(*(partition_config.stream_in))) {
    std::cerr << "Error opening " << graph_filename << std::endl;
    exit(1);
  }
  std::vector<std::string> *lines;

  lines = new std::vector<std::string>(1);
  std::getline(*(partition_config.stream_in), (*lines)[0]);

  // skip comments
  while ((*lines)[0][0] == '%') {
    std::getline(*(partition_config.stream_in), (*lines)[0]);
  }

  std::stringstream ss((*lines)[0]);
  ss >> partition_config.remaining_stream_nodes;
  ss >> partition_config.remaining_stream_edges;
  ss >> partition_config.remaining_stream_ew;

  switch (partition_config.remaining_stream_ew) {
  case 1:
    partition_config.read_ew = true;
    break;
  case 10:
    partition_config.read_nw = true;
    break;
  case 11:
    partition_config.read_ew = true;
    partition_config.read_nw = true;
    break;
  }

  partition_config.total_edges = partition_config.remaining_stream_edges;
  partition_config.total_nodes = partition_config.remaining_stream_nodes;

    if (((partition_config.one_pass_algorithm == ONEPASS_HASHING) ||
         (partition_config.one_pass_algorithm == ONEPASS_HASHING_CRC32)) && !partition_config.evaluate) {
    } else {
        if (partition_config.stream_nodes_assign == NULL &&
            partition_config.rle_length == -1) {
            std::cout << "Using std::vector of size " << partition_config.total_nodes
                      << std::endl;
            partition_config.stream_nodes_assign = new std::vector<PartitionID>(
                    partition_config.remaining_stream_nodes, INVALID_PARTITION);
        } else if (partition_config.rle_length == -2) {
            std::cout << "Using external memory PQ with time-forward processing"
                      << std::endl;
            partition_config.external_pq_partition_assign = new ExternalPQ();
            if (partition_config.evaluate) {
                std::cout << "Running eval.. (not accurate for memory consumption)"
                          << std::endl;
                partition_config.stream_nodes_assign = new std::vector<PartitionID>(
                        partition_config.remaining_stream_nodes, INVALID_PARTITION);
            }
        }
    }

  if (partition_config.stream_blocks_weight == NULL) {
    partition_config.stream_blocks_weight =
        new std::vector<NodeWeight>(partition_config.k, 0);
  }
  partition_config.total_stream_nodeweight = 0;
  partition_config.total_stream_nodecounter = 0;
  partition_config.stream_n_nodes = partition_config.remaining_stream_nodes;

  if (partition_config.num_streams_passes >
      1 + partition_config.restream_number) {
    partition_config.stream_total_upperbound = ceil(
        ((100 + 1.5 * partition_config.imbalance) / 100.) *
        (partition_config.remaining_stream_nodes / (double)partition_config.k));
  } else {
    partition_config.stream_total_upperbound = ceil(
        ((100 + partition_config.imbalance) / 100.) *
        (partition_config.remaining_stream_nodes / (double)partition_config.k));
  }

  partition_config.fennel_alpha =
      partition_config.remaining_stream_edges *
      std::pow(partition_config.k, partition_config.fennel_gamma - 1) /
      (std::pow(partition_config.remaining_stream_nodes,
                partition_config.fennel_gamma));

  partition_config.fennel_alpha_gamma =
      partition_config.fennel_alpha * partition_config.fennel_gamma;

  if (partition_config.full_stream_mode && !partition_config.restream_number) {
    partition_config.quotient_nodes = 0;
  } else {
    partition_config.quotient_nodes = partition_config.k;
  }

  total_edge_cut = 0;
  qap = 0;
  if (partition_config.stream_buffer_len ==
      0) { // signal of partial restream standard buffer size
    partition_config.stream_buffer_len = (LongNodeID)ceil(
        partition_config.remaining_stream_nodes / (double)partition_config.k);
  }
  partition_config.nmbNodes = MIN(partition_config.stream_buffer_len,
                                  partition_config.remaining_stream_nodes);
  partition_config.n_batches = ceil(partition_config.remaining_stream_nodes /
                                    (double)partition_config.nmbNodes);
  partition_config.curr_batch = 0;

  delete lines;
}

void graph_io_stream::configureGeneratedStream(PartitionConfig &partition_config,
                                          std::string graph_filename,
                                          uint64_t num_edges,
                                          EdgeWeight &total_edge_cut,
                                          EdgeWeight &qap) {
    if (partition_config.stream_in != NULL) {
        delete partition_config.stream_in;
    }
    partition_config.remaining_stream_nodes = partition_config.nodes_to_generate;
    partition_config.remaining_stream_edges = num_edges;

    partition_config.total_edges = partition_config.remaining_stream_edges;
    partition_config.total_nodes = partition_config.remaining_stream_nodes;

    std::cout << "No. of generated edges (estimated): " << num_edges << std::endl;

    if (((partition_config.one_pass_algorithm == ONEPASS_HASHING) ||
         (partition_config.one_pass_algorithm == ONEPASS_HASHING_CRC32)) && !partition_config.evaluate) {
    } else {
        if (partition_config.stream_nodes_assign == NULL &&
            partition_config.rle_length == -1) {
            std::cout << "Using std::vector of size " << partition_config.total_nodes
                      << std::endl;
            partition_config.stream_nodes_assign = new std::vector<PartitionID>(
                    partition_config.remaining_stream_nodes, INVALID_PARTITION);
        } else if (partition_config.rle_length == -2) {
            std::cout << "Using external memory PQ with time-forward processing"
                      << std::endl;
            partition_config.external_pq_partition_assign = new ExternalPQ();
            if (partition_config.evaluate) {
                std::cout << "Running eval.. (not accurate for memory consumption)"
                          << std::endl;
                partition_config.stream_nodes_assign = new std::vector<PartitionID>(
                        partition_config.remaining_stream_nodes, INVALID_PARTITION);
            }
        }
    }

    if (partition_config.stream_blocks_weight == NULL) {
        partition_config.stream_blocks_weight =
                new std::vector<NodeWeight>(partition_config.k, 0);
    }
    partition_config.total_stream_nodeweight = 0;
    partition_config.total_stream_nodecounter = 0;
    partition_config.stream_n_nodes = partition_config.remaining_stream_nodes;

    if (partition_config.num_streams_passes >
        1 + partition_config.restream_number) {
        partition_config.stream_total_upperbound = ceil(
                ((100 + 1.5 * partition_config.imbalance) / 100.) *
                (partition_config.remaining_stream_nodes / (double)partition_config.k));
    } else {
        partition_config.stream_total_upperbound = ceil(
                ((100 + partition_config.imbalance) / 100.) *
                (partition_config.remaining_stream_nodes / (double)partition_config.k));
    }

    partition_config.fennel_alpha =
            partition_config.remaining_stream_edges *
            std::pow(partition_config.k, partition_config.fennel_gamma - 1) /
            (std::pow(partition_config.remaining_stream_nodes,
                      partition_config.fennel_gamma));

    partition_config.fennel_alpha_gamma =
            partition_config.fennel_alpha * partition_config.fennel_gamma;

    if (partition_config.full_stream_mode && !partition_config.restream_number) {
        partition_config.quotient_nodes = 0;
    } else {
        partition_config.quotient_nodes = partition_config.k;
    }

    total_edge_cut = 0;
    qap = 0;
    if (partition_config.stream_buffer_len ==
        0) { // signal of partial restream standard buffer size
        partition_config.stream_buffer_len = (LongNodeID)ceil(
                partition_config.remaining_stream_nodes / (double)partition_config.k);
    }
    partition_config.nmbNodes = MIN(partition_config.stream_buffer_len,
                                    partition_config.remaining_stream_nodes);
    partition_config.n_batches = ceil(partition_config.remaining_stream_nodes /
                                      (double)partition_config.nmbNodes);
    partition_config.curr_batch = 0;
}

void graph_io_stream::streamEvaluatePartition(PartitionConfig &config,
                                              const std::string &filename,
                                              EdgeWeight &edgeCut,
                                              EdgeWeight &qap,
                                              const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments) {
  std::vector<std::vector<LongNodeID>> *input;
  std::vector<std::string> *lines;
  lines = new std::vector<std::string>(1);
  LongNodeID node_counter = 0;
  buffered_input *ss2 = NULL;
  std::string line;
  std::ifstream in(filename.c_str());
  if (!in) {
    std::cerr << "Error opening " << filename << std::endl;
    return 1;
  }
  long nmbNodes;
  long nmbEdges;
  int ew = 0;
  std::getline(in, (*lines)[0]);
  while ((*lines)[0][0] == '%')
    std::getline(in, (*lines)[0]); // a comment in the file

  std::stringstream ss((*lines)[0]);
  ss >> nmbNodes;
  ss >> nmbEdges;
  ss >> ew;
  bool read_ew = false;
  bool read_nw = false;
  if (ew == 1) {
    read_ew = true;
  } else if (ew == 11) {
    read_ew = true;
    read_nw = true;
  } else if (ew == 10) {
    read_nw = true;
  }
  NodeID target;
  NodeWeight total_nodeweight = 0;
  EdgeWeight total_edgeweight = 0;
  edgeCut = 0;
  qap = 0;
  std::vector <NodeID> block_weights(config.k, 0);

  while (std::getline(in, (*lines)[0])) {
    if ((*lines)[0][0] == '%')
      continue; // a comment in the file
    NodeID node = node_counter++;

    PartitionID partitionIDSource;

    if (config.rle_length == -1 ||
        (config.evaluate && config.rle_length == -2)) {
      partitionIDSource = (*config.stream_nodes_assign)[node];
    } else if (config.rle_length == 0) {
      partitionIDSource = block_assignments->GetValueByIndex(node);
    } else {
        partitionIDSource = block_assignments->GetValueByBatchIndex(node / config.rle_length,
                                                                    (node % config.rle_length));
    }
    input = new std::vector<std::vector<LongNodeID>>(1);
    ss2 = new buffered_input(lines);
    ss2->simple_scan_line((*input)[0]);
    std::vector<LongNodeID> &line_numbers = (*input)[0];
    LongNodeID col_counter = 0;
    block_weights[partitionIDSource]++;

    NodeWeight weight = 1;
    if (read_nw) {
      weight = line_numbers[col_counter++];
      total_nodeweight += weight;
    }
    while (col_counter < line_numbers.size()) {
      target = line_numbers[col_counter++];
      target = target - 1;
      EdgeWeight edge_weight = 1;
      if (read_ew) {
        edge_weight = line_numbers[col_counter++];
      }
      total_edgeweight += edge_weight;
      PartitionID partitionIDTarget;
      if (config.rle_length == -1 ||
          (config.evaluate && config.rle_length == -2)) {
        partitionIDTarget = (*config.stream_nodes_assign)[target];
      } else if (config.rle_length == 0) {
        partitionIDTarget = block_assignments->GetValueByIndex(target);
      } else {
        partitionIDTarget = block_assignments->GetValueByBatchIndex(target / config.rle_length, (target % config.rle_length));
      }

      if (partitionIDSource != partitionIDTarget) {
        edgeCut += edge_weight;
      }
    }
    (*lines)[0].clear();
    delete ss2;
    delete input;
    if (in.eof()) {
      break;
    }
  }
  edgeCut = edgeCut / 2; // Since every edge is counted twice
  qap = edgeCut;
  delete lines;

    double max = -1;
    double total_weight = 0;
    double balance;

    for (int i = 0; i < config.k; i++) {
        if (block_weights[i] > max) {
            max = block_weights[i];
        }
        if (block_weights[i] > config.stream_total_upperbound) {
            std::cout << "Partition is imbalanced" << std::endl;
            exit(1);
        }
        total_weight += block_weights[i];
    }
    double balance_part_weight = ceil(total_weight / (double) config.k);
    balance = max / balance_part_weight;
    std::cout << "Computed balance: " << balance << std::endl;
}

void graph_io_stream::writePartitionStream(PartitionConfig &config,
                                           const std::string &filename,
                                           const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments) {
  std::ofstream f(filename.c_str());
  std::cout << "writing partition to " << filename << " ... " << std::endl;

  for (int node = 0; node < config.total_nodes; node++) {
    if (config.rle_length == -1 ||
        (config.evaluate && config.rle_length == -2)) {
      f << (*config.stream_nodes_assign)[node] << "\n";
    } else if (config.rle_length == 0) {
      f << block_assignments->GetValueByIndex(node) << "\n";
    } else {
      f << block_assignments->GetValueByBatchIndex(node / config.rle_length, (node % config.rle_length))
        << "\n";
    }
  }

  f.close();
}

void graph_io_stream::readPartition(PartitionConfig &config,
                                    const std::string &filename) {
  std::string line;

  // open file for reading
  std::ifstream in(filename.c_str());
  if (!in) {
    std::cerr << "Error opening file" << filename << std::endl;
    return 1;
  }

  PartitionID max = 0;
  for (auto &node : (*config.stream_nodes_assign)) {
    // fetch current line
    std::getline(in, line);
    while (line[0] == '%') { // Comments
      std::getline(in, line);
    }
    node = (PartitionID)atol(line.c_str());
    (*config.stream_blocks_weight)[node] += 1;

    if (node > max)
      max = node;
  }

  config.k = max + 1;
  in.close();
}
