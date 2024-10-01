/******************************************************************************
 * streammultisection.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj <marcelofaraj@gmail.com>
 *****************************************************************************/

#include <argtable3.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <vector>
#include <cmath>

#include "data_structure/ExternalPQ.h"
#include "data_structure/graph_access.h"
#include "graph_io_stream.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "partition/partition_config.h"
#include "quality_metrics.h"
#include "timer.h"
#include "tools/random_functions.h"

#include "partition/onepass_partitioning/fennel.h"
#include "partition/onepass_partitioning/fennel_approx_sqrt.h"
#include "partition/onepass_partitioning/ldg.h"
#include "partition/onepass_partitioning/vertex_partitioning.h"

#include "FlatBufferWriter.h"
#include "Stream_CPI_Info_generated.h"
#include "cpi/run_length_compression.hpp"

#include "data_structure/compression_vectors/CompressionDataStructure.h"
#include "data_structure/compression_vectors/RunLengthCompressionVector.h"
#include "data_structure/compression_vectors/BatchRunLengthCompression.h"

#include <kagen.h>
#include <mpi.h>
#include <unistd.h>

#define MIN(A, B) (((A) < (B)) ? (A) : (B))
#define MAX(A, B) (((A) > (B)) ? (A) : (B))

// Struct to store captured values
struct CapturedValues {
    std::size_t space_in_bytes;
    std::size_t uncompressed_space_in_bytes;
    double space_in_mib;
    double relative;
};

void initialize_onepass_partitioner(PartitionConfig &config,
                                    vertex_partitioning *&onepass_partitioner);

long getMaxRSS();

void streamEvaluatePartitionGenerated(PartitionConfig &config,
                                      const std::string &filename,
                                      EdgeWeight &edgeCut,
                                      EdgeWeight &qap,
                                      kagen::StreamingGenerator &streamGenerator,
                                      const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments);

std::string extractBaseFilename(const std::string &fullPath);

std::ostream &cout_redirect();

CapturedValues parseCapturedValues(const std::string &output_str);

std::ostringstream redirected_cout;

int main(int argn, char **argv) {
    PartitionConfig config;
    std::string graph_filename;
    /* LINE_BUFFER lines = NULL; */
    std::vector<std::vector<LongNodeID>> *input = NULL;
    timer t, processing_t, io_t;
    EdgeWeight total_edge_cut = 0;
    int counter = 0;
    double global_mapping_time = 0;
    double buffer_mapping_time = 0;
    double buffer_io_time = 0;
    double total_time = 0;
    quality_metrics qm;
    EdgeWeight qap = 0;
    double balance = 0;
    int full_stream_count = 0;
    double total_nodes = 0;

    bool is_graph_weighted = false;
    bool suppress_output = false;
    bool recursive = false;

    int ret_code =
            parse_parameters(argn, argv, config, graph_filename, is_graph_weighted,
                             suppress_output, recursive);

    if (ret_code) {
        return 0;
    }

    std::streambuf *backup = std::cout.rdbuf();
    std::ofstream ofs;
    ofs.open("/dev/null");
    if (suppress_output) {
        std::cout.rdbuf(ofs.rdbuf());
    }

    srand(config.seed);
    random_functions::setSeed(config.seed);

    config.LogDump(stdout);
    config.stream_input = true;

    bool already_fully_partitioned;
    pid_t processID = getpid();
    config.nodes_to_generate = (pow(2, config.nodes_to_generate));
    config.kagen_chunk_count = (pow(2, config.kagen_chunk_count));
    //config.nodes_to_generate = 10000000000;
    if (config.nodes_to_generate < 0 || config.nodes_to_generate > std::numeric_limits<NodeID>::max()) {
	throw std::out_of_range("Value is out of range.");
    }


    vertex_partitioning *onepass_partitioner = NULL;
    initialize_onepass_partitioner(config, onepass_partitioner);

    // container for storing block assignments used by Fennel
    std::shared_ptr<CompressionDataStructure<PartitionID>> block_assignments;
    config.streaming_graph_generation = true;

    MPI_Init(NULL, NULL);
    std::cout << "Nodes to generate: " << config.nodes_to_generate << std::endl;
    kagen::StreamingGenerator streamGenerator(MPI_COMM_WORLD, config.kagen_chunk_count);
    if(config.rgg2d) {
        std::cout << "KaGen radius: " << config.kagen_r << std::endl;
        streamGenerator.setupConfig_RGG2D(config.nodes_to_generate, 0, config.kagen_r, false);
    } else if (config.ba) {
		std::cout << "KaGen average degree: " << config.kagen_d_ba << std::endl; 
        streamGenerator.setupConfig_BA(config.nodes_to_generate, 0, config.kagen_d_ba, false, false);
    } else if (config.rgg3d) {
	    std::cout << "KaGen radius: " << config.kagen_r << std::endl;
        streamGenerator.setupConfig_RGG3D(config.nodes_to_generate, 0, config.kagen_r, false);
	} else { // (config.rhg) {
	    std::cout << "KaGen average degree : " << config.kagen_r << ", KaGen gamma: " << config.kagen_gamma << std::endl;
        streamGenerator.setupConfig_RHG(config.nodes_to_generate, 0, config.kagen_d_rhg, config.kagen_gamma, false);
	}
    streamGenerator.setRandomSeed(config.seed);
    streamGenerator.setupChunkGeneration(MPI_COMM_WORLD);
    uint64_t estimated_edge_count = streamGenerator.estimate_edges();


    int &passes = config.num_streams_passes;
    for (config.restream_number = 0; config.restream_number < passes;
         config.restream_number++) {

        io_t.restart();
        graph_io_stream::configureGeneratedStream(config, graph_filename, estimated_edge_count, total_edge_cut, qap);
        buffer_io_time += io_t.elapsed();

        // set up block assignment container based on algorithm configuration
        if(config.rle_length == 0) {
            block_assignments = std::make_shared<RunLengthCompressionVector<PartitionID>>();
        }
        else if (config.rle_length > 0) {
            block_assignments = std::make_shared<BatchRunLengthCompression<PartitionID>>((config.total_nodes /
                                                                                          config.rle_length) + 1);
        }

        onepass_partitioner->instantiate_blocks(config.remaining_stream_nodes,
                                                config.remaining_stream_edges,
                                                config.k, config.imbalance);

        for (int i = 0; i < config.parallel_nodes; i++) {
            config.all_blocks_to_keys[i].resize(config.k);
            for (auto &b : config.all_blocks_to_keys[i]) {
                b = INVALID_PARTITION;
            }
            config.neighbor_blocks[i].resize(config.k);
            config.next_key[i] = 0;
        }

        processing_t.restart();

        for (LongNodeID curr_node = 0; curr_node < config.nodes_to_generate; curr_node++) {
            int my_thread = 0;
            io_t.restart();
            std::vector<NodeID> generated_neighbors; 
            if ((config.one_pass_algorithm != ONEPASS_HASHING) &&
                (config.one_pass_algorithm != ONEPASS_HASHING_CRC32)) {
		   	streamGenerator.streamVertex(curr_node + 1, MPI_COMM_WORLD, generated_neighbors);
                graph_io_stream::generateNeighborhood(config, generated_neighbors, input);
            }
            buffer_io_time += io_t.elapsed();
            // ***************************** perform partitioning
            // ***************************************
            t.restart();
            graph_io_stream::readNodeOnePass(config, curr_node, my_thread, input, block_assignments,
                                             onepass_partitioner);
            PartitionID block = onepass_partitioner->solve_node(
                    curr_node, 1, config.previous_assignment, config.kappa, my_thread);

	    if (((config.one_pass_algorithm == ONEPASS_HASHING) ||
                 (config.one_pass_algorithm == ONEPASS_HASHING_CRC32)) && !config.evaluate) {
            } else {
                if (config.set_part_zero) {
                    block = 0;
                }
                if (config.rle_length == -1) {
                    (*config.stream_nodes_assign)[curr_node] = block;
                } else if (config.rle_length == 0) {
                    block_assignments->Append(block);
                } else if (config.rle_length == -2) {
                    // add blocks to PQ
                    auto &read_ew = config.read_ew;
                    auto &read_nw = config.read_nw;
                    NodeWeight weight;
                    LongNodeID cursor = (config.ram_stream) ? curr_node : 0;

                    std::vector <LongNodeID> &line_numbers = (*input)[cursor];
                    LongNodeID col_counter = 0;
                    weight = (read_nw) ? line_numbers[col_counter++] : 1;

                    while (col_counter < line_numbers.size()) {
                        LongNodeID target = line_numbers[col_counter++];
                        EdgeWeight edge_weight = (read_ew) ? line_numbers[col_counter++] : 1;
                        if (curr_node < (target - 1)) {
                            (*config.external_pq_partition_assign)
                                    .appendToExtPQ(block, target - 1);
                        }
                    }
                    if (config.evaluate) {
                        (*config.stream_nodes_assign)[curr_node] = block;
                    }
                } else {
                    block_assignments->BatchAppend(curr_node / config.rle_length, block);
                }
            }
            if (!config.ram_stream) {
                delete input;
            }

            config.previous_assignment = block;
            (*config.stream_blocks_weight)[block] += 1;

            global_mapping_time += t.elapsed();
        }
        total_time += processing_t.elapsed();

        if (config.ram_stream) {
            delete input;
            /* delete lines; */
        }
    }

    long overall_max_RSS = getMaxRSS();
    std::string baseFilename = extractBaseFilename(graph_filename);
    if (!config.suppress_output) {
        if (((config.one_pass_algorithm == ONEPASS_HASHING) ||
             (config.one_pass_algorithm == ONEPASS_HASHING_CRC32)) && !config.evaluate) {
            total_edge_cut = 0;
        } else {
            if (config.rle_length != -2 || config.evaluate) {
                streamEvaluatePartitionGenerated(config, graph_filename,
                                                 total_edge_cut, qap, streamGenerator,
                                                 block_assignments);
            }
        }
        balance = qm.balance_full_stream(*config.stream_blocks_weight);
    }

    CapturedValues capturedValues;
    //if (config.rle_length == -1) {
    //    std::streambuf *original_cout_buffer = std::cout.rdbuf();
    //    std::cout.rdbuf(redirected_cout.rdbuf());
    //    std::cout << "Performing RLE compression test..." << std::endl;
    //    cpi::RunLengthCompression rlc(*config.stream_nodes_assign);
    //    auto p_id = rlc[42];    // access partition ids;
    //    rlc.print_statistics(); // print statistics
    //    std::string output_str = redirected_cout.str();
    //    capturedValues = parseCapturedValues(output_str);
    //    std::cout.rdbuf(original_cout_buffer);
    //} else {
        capturedValues.space_in_bytes = 0;
        capturedValues.uncompressed_space_in_bytes = 0;
        capturedValues.space_in_mib = 0;
        capturedValues.relative = 0;
    //}

    // write the partition to the disc
    std::stringstream filename;
    if (!config.filename_output.compare("")) {
        filename << baseFilename << "_" << config.k;
    } else {
        filename << config.filename_output;
    }

    if (!config.suppress_file_output) {
        if (config.rle_length != -2 || config.evaluate) {
            graph_io_stream::writePartitionStream(config, filename.str(), block_assignments);
        }
    } else {
        std::cout << "No partition will be written as output." << std::endl;
    }

    if (config.write_results) {
        FlatBufferWriter fb_writer;
        fb_writer.updateResourceConsumption(buffer_io_time, global_mapping_time,
                                            total_time, overall_max_RSS);
        fb_writer.updatePartitionMetrics(total_edge_cut, balance);
        fb_writer.updateCompressionStatistics(
                capturedValues.space_in_bytes,
                capturedValues.uncompressed_space_in_bytes, capturedValues.space_in_mib,
                capturedValues.relative);
        fb_writer.write(baseFilename, config);
    }

    return 0;
}

void initialize_onepass_partitioner(PartitionConfig &config,
                                    vertex_partitioning *&onepass_partitioner) {
    switch (config.one_pass_algorithm) {
        case ONEPASS_HASHING:
        case ONEPASS_HASHING_CRC32:
            onepass_partitioner = new vertex_partitioning(
                    0, config.k - 1, config.stream_rec_bisection_base,
                    config.parallel_nodes, true);
            break;
        case ONEPASS_LDG:
            onepass_partitioner =
                    new onepass_ldg(0, config.k - 1, config.stream_rec_bisection_base,
                                    config.parallel_nodes, false);
            break;
        case ONEPASS_FENNEL:
            onepass_partitioner =
                    new onepass_fennel(0, config.k - 1, config.stream_rec_bisection_base,
                                       config.parallel_nodes, false, config.fennel_gamma);
            break;
        case ONEPASS_FENNEL_APPROX_SQRT:
        default:
            onepass_partitioner = new onepass_fennel_approx_sqrt(
                    0, config.k - 1, config.stream_rec_bisection_base,
                    config.parallel_nodes, false, config.fennel_gamma);
            break;
    }
}

long getMaxRSS() {
    struct rusage usage;

    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        // The maximum resident set size is in kilobytes
        return usage.ru_maxrss;
    } else {
        std::cerr << "Error getting resource usage information." << std::endl;
        // Return a sentinel value or handle the error in an appropriate way
        return -1;
    }
}

// Redirect cout to the stringstream
std::ostream &cout_redirect() {
    static std::ostream cout_redirector(redirected_cout.rdbuf());
    return cout_redirector;
}

// Function to parse the captured values from the redirected output
CapturedValues parseCapturedValues(const std::string &output_str) {
    CapturedValues values;

    // Extract values using stream extraction
    std::istringstream stream(output_str);

    // Ignore text up to '=' and then extract values
    stream.ignore(std::numeric_limits<std::streamsize>::max(), '=');
    stream >> values.space_in_bytes;

    stream.ignore(std::numeric_limits<std::streamsize>::max(), '=');
    stream >> values.uncompressed_space_in_bytes;

    stream.ignore(std::numeric_limits<std::streamsize>::max(), '=');
    stream >> values.space_in_mib;

    stream.ignore(std::numeric_limits<std::streamsize>::max(), '=');
    stream >> values.relative;

    return values;
}

// Function to extract the base filename without path and extension
std::string extractBaseFilename(const std::string &fullPath) {
    size_t lastSlash = fullPath.find_last_of('/');
    size_t lastDot = fullPath.find_last_of('.');

    if (lastSlash != std::string::npos) {
        // Found a slash, extract the substring after the last slash
        return fullPath.substr(lastSlash + 1, lastDot - lastSlash - 1);
    } else {
        // No slash found, just extract the substring before the last dot
        return fullPath.substr(0, lastDot);
    }
}

void streamEvaluatePartitionGenerated(PartitionConfig &config,
                                                       const std::string &filename,
                                                       EdgeWeight &edgeCut,
                                                       EdgeWeight &qap,
                                                       kagen::StreamingGenerator &streamGenerator2,
                                                       const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments) {
    long nmbNodes = config.total_nodes;
    long nmbEdges = config.total_edges;
    int ew = 0;
    std::cout << "Running eval..." << std::endl;
    //MPI_Init(NULL, NULL);
    std::cout << "Nodes to generate: " << config.nodes_to_generate << std::endl;
    kagen::StreamingGenerator streamGenerator(MPI_COMM_WORLD, config.kagen_chunk_count);
    if(config.rgg2d) {
        std::cout << "KaGen radius: " << config.kagen_r << std::endl;
        streamGenerator.setupConfig_RGG2D(config.nodes_to_generate, 0, config.kagen_r, false);
    } else if (config.ba) {
        std::cout << "KaGen average degree: " << config.kagen_d_ba << std::endl;
        streamGenerator.setupConfig_BA(config.nodes_to_generate, 0, config.kagen_d_ba, false, false);
    } else if (config.rgg3d) {
        std::cout << "KaGen radius: " << config.kagen_r << std::endl;
        streamGenerator.setupConfig_RGG3D(config.nodes_to_generate, 0, config.kagen_r, false);
    } else { // (config.rhg) {
        std::cout << "KaGen average degree : " << config.kagen_r << ", KaGen gamma: " << config.kagen_gamma << std::endl;
        streamGenerator.setupConfig_RHG(config.nodes_to_generate, 0, config.kagen_d_rhg, config.kagen_gamma, false);
    }
    streamGenerator.setRandomSeed(config.seed);
    streamGenerator.setupChunkGeneration(MPI_COMM_WORLD);
    pid_t processID = getpid();
    //std::cout << "----------" << std::endl;
    NodeID target;
    NodeWeight total_nodeweight = 0;
    EdgeWeight total_edgeweight = 0;
    edgeCut = 0;
    qap = 0;
    LongNodeID real_edge_count = 0;

    for(NodeID node = 0; node < config.nodes_to_generate; ++node) {
        //std::cout << "Node " << node << ": ";
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

        LongNodeID col_counter = 0;

        NodeWeight weight = 1;
        std::vector<NodeID> generated_neighbors;
	   	streamGenerator.streamVertex(node + 1, MPI_COMM_WORLD, generated_neighbors);
        for(auto & neighbor: generated_neighbors) {
            real_edge_count++;
		target = neighbor - 1;
            EdgeWeight edge_weight = 1;
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
    }
    edgeCut = edgeCut / 2; // Since every edge is counted twice
    qap = edgeCut;
    config.total_edges=real_edge_count;
}
