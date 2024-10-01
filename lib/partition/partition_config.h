/******************************************************************************
 * partition_config.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PARTITION_CONFIG_DI1ES4T0
#define PARTITION_CONFIG_DI1ES4T0

#include "definitions.h"
#include "data_structure/buffered_map.h"
#include "data_structure/single_adj_list.h"
#include "data_structure/RunLengthEncodedVector.h"
#include "cpi/run_length_compression.hpp"
#include "data_structure/ExternalPQ.h"

#include <omp.h>
#include <memory>

typedef struct {
    PartitionID block;
    int gain;
    int degree;
} DELTA;

class matrix;

// Configuration for the partitioning.
struct PartitionConfig {
    PartitionConfig() {}

    PermutationQuality permutation_quality;
    ImbalanceType imbalance;

    double time_limit;
    double epsilon;

    std::string input_partition;
    int seed;
    bool balance_edges;
    // number of blocks the graph should be partitioned in
    PartitionID k;
    bool compute_vertex_separator;
    bool only_first_level;
    bool use_balance_singletons;
    int amg_iterations;
    std::string graph_filename;
    std::string filename_output;
    std::string output_path;
    double balance_factor;

    //=======================================
    //========== Stream Partition ===========
    //=======================================

    bool stream_input;
    LongNodeID stream_buffer_len;
    LongNodeID rle_length;
    double kappa;
    PartitionID previous_assignment;
    LongNodeID remaining_stream_nodes;
    LongEdgeID remaining_stream_edges;
    LongEdgeID total_edges;
    LongNodeID total_nodes;
    LongNodeID max_degree;
    bool set_part_zero;
    bool write_results;
    int remaining_stream_ew;
    LongNodeID total_stream_nodeweight;
    LongNodeID total_stream_nodecounter;
    LongNodeID stream_assigned_nodes;
    LongNodeID stream_n_nodes;
    std::ifstream *stream_in;
    LongNodeID lower_global_node;
    LongNodeID upper_global_node;
    std::size_t uncompressed_runs = 64;
    std::vector <PartitionID> *stream_nodes_assign;
    ExternalPQ *external_pq_partition_assign;
    std::vector <NodeWeight> *stream_blocks_weight;
    LongNodeID nmbNodes;
    std::vector <std::vector<EdgeWeight>> *degree_nodeBlock;
    std::vector <std::vector<EdgeWeight>> *ghostDegree_nodeBlock;
    int one_pass_algorithm;
    bool full_stream_mode;
    LongNodeID stream_total_upperbound;
    double fennel_gamma;
    double fennel_alpha;
    double fennel_alpha_gamma;
    bool use_fennel_objective; // maps global blocks to current stream blocks
    int fennel_dynamics;
    bool ram_stream;
    bool evaluate;
    bool fennel_contraction;
    int fennel_batch_order;
    int quotient_nodes;
    int lhs_nodes;
    bool stream_initial_bisections;
    LongNodeID n_batches;
    int curr_batch;
    double stream_global_epsilon;
    bool stream_output_progress;
    double batch_inbalance;
    bool skip_outer_ls;
    bool use_fennel_edgecut_objectives;
    std::vector <PartitionID> one_pass_neighbor_blocks;

    // KaGen Streaming
    bool streaming_graph_generation;
    LongNodeID nodes_to_generate;
    NodeID kagen_chunk_count;
    int kagen_d_ba;
    bool rgg2d;
    double kagen_r;
	bool rgg3d; 
	bool rdg2d; 
	bool rdg3d; 
	bool ba; 
	bool rhg; 
	double kagen_d_rhg; 
	double kagen_gamma; 

    // Initial partition via growing multiple BFS trees
    bool initial_part_multi_bfs;
    int multibfs_tries;

    // Initial partitioning via Fennel on the coarsest level
    int initial_part_fennel_tries;

    // Restreaming and partial restreaming
    int num_streams_passes;
    int restream_number;
    bool restream_vcycle;

    int xxx;
    double *t1;
    double *t2;
    double *t3;

    //=======================================
    //============= Stream Map ==============
    //=======================================

    double specify_alpha_multiplier;
    bool stream_multisection;
    std::vector <std::vector<NodeWeight>> *stream_modules_weight;
    std::vector <std::vector<std::vector < EdgeWeight>>> *
    degree_nodeLayerModule;
    std::vector <std::vector<std::vector < EdgeWeight>>> *
    ghostDegree_nodeLayerModule;
    int stream_weighted_msec_type;
    bool onepass_pipelined_input;
    bool onepass_simplified_input;
    bool multicore_pipeline;
    pipelist_nodes *nodes_pipeline;
    PartitionID pipeline_stages;
    int parallel_nodes;
    PartitionID hashify_layers;
    int fast_alg;
    NodeWeight one_pass_my_weight;
    std::vector <std::vector<std::pair < PartitionID, EdgeWeight>>>
    neighbor_blocks;
    std::vector <std::vector<PartitionID>> all_blocks_to_keys;
    std::vector <PartitionID> next_key;
    bool stream_rec_bisection;
    PartitionID stream_rec_bisection_base;
    bool stream_rec_biss_orig_alpha;
    PartitionID non_hashified_layers;
    float percent_non_hashified_layers;

    LongEdgeID edges;
    bool read_ew;
    bool read_nw;
    bool suppress_output;
    bool suppress_file_output;

    void LogDump(FILE *out) const {
    }
};


#endif /* end of include guard: PARTITION_CONFIG_DI1ES4T0 */
