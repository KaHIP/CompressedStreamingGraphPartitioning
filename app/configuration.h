/******************************************************************************
 * configuration.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/


#ifndef CONFIGURATION_3APG5V7Z
#define CONFIGURATION_3APG5V7Z

#include "partition/partition_config.h"

class configuration {
public:
    configuration() {};

    virtual ~configuration() {};

    void strong(PartitionConfig &config);

    void eco(PartitionConfig &config);

    void fast(PartitionConfig &config);

    void strong_separator(PartitionConfig &config);

    void eco_separator(PartitionConfig &config);

    void fast_separator(PartitionConfig &config);

    void standard(PartitionConfig &config);

    void fastsocial(PartitionConfig &config);

    void ecosocial(PartitionConfig &config);

    void integrated_mapping(PartitionConfig &partition_config);

    void fsocialmap(PartitionConfig &partition_config);

    void stream_partition(PartitionConfig &partition_config);

    void stream_map(PartitionConfig &partition_config);

};

inline void configuration::strong(PartitionConfig &partition_config) {
    standard(partition_config);
}

inline void configuration::eco(PartitionConfig &partition_config) {
    standard(partition_config);
}

inline void configuration::fast(PartitionConfig &partition_config) {
    standard(partition_config);
}

inline void configuration::standard(PartitionConfig &partition_config) {
    partition_config.filename_output = "";
    partition_config.seed = 0;

    partition_config.epsilon                                = 3;
    partition_config.imbalance                              = 3;

    partition_config.balance_edges = false;

    partition_config.time_limit = 0;

    // Stream Partition
    partition_config.stream_input = false;
    partition_config.max_degree = 0;
    partition_config.set_part_zero = false;
    partition_config.write_results = false;
    partition_config.stream_buffer_len = 32768;
    partition_config.rle_length = -1;
    partition_config.uncompressed_runs = 64;
    partition_config.kappa = 1;
    partition_config.previous_assignment = -1;
    partition_config.remaining_stream_nodes = UNDEFINED_LONGNODE;
    partition_config.remaining_stream_edges = UNDEFINED_LONGEDGE;
    partition_config.remaining_stream_ew = 0;
    partition_config.total_stream_nodeweight = 0;
    partition_config.total_stream_nodecounter = 0;
    partition_config.stream_assigned_nodes = 0;
    partition_config.stream_n_nodes = 0;
    partition_config.stream_in = NULL;
    partition_config.lower_global_node = 1;
    partition_config.upper_global_node = partition_config.k;
    partition_config.stream_nodes_assign = NULL;
    partition_config.stream_blocks_weight = NULL;
    partition_config.nmbNodes = 0;
    partition_config.degree_nodeBlock = NULL;
    partition_config.one_pass_algorithm = ONEPASS_FENNEL_APPROX_SQRT;
    partition_config.full_stream_mode = false;
    partition_config.stream_total_upperbound = 0;
    partition_config.fennel_gamma = 1.5;
    partition_config.fennel_alpha = 1;
    partition_config.fennel_alpha_gamma = partition_config.fennel_alpha * partition_config.fennel_gamma;
    partition_config.use_fennel_objective = false;
    partition_config.fennel_dynamics = FENNELADP_ORIGINAL;
    partition_config.ram_stream = false;
    partition_config.evaluate = false;
    partition_config.fennel_contraction = false;
    partition_config.fennel_batch_order = FENN_ORDER_UNCHANGED;
    partition_config.quotient_nodes = 0;
    partition_config.lhs_nodes = 0;
    partition_config.stream_initial_bisections = false;
    partition_config.n_batches = 1;
    partition_config.curr_batch = 1;
    partition_config.stream_global_epsilon = (partition_config.imbalance) / 100.;
    partition_config.stream_output_progress = false;
    partition_config.num_streams_passes = 1;
    partition_config.restream_number = 0;
    partition_config.restream_vcycle = false;
    partition_config.batch_inbalance = 20;
    partition_config.skip_outer_ls = false;
    partition_config.use_fennel_edgecut_objectives = false;
    partition_config.xxx = 4;
    partition_config.suppress_output = false;
    partition_config.suppress_file_output = false;

    // KaGen Streaming
    partition_config.streaming_graph_generation=false;
    partition_config.nodes_to_generate = 0;
    partition_config.kagen_chunk_count = 0;
    partition_config.kagen_d_ba = 0;
    partition_config.kagen_r = 0;
    partition_config.rgg2d = false;
	partition_config.rgg3d = false; 
	partition_config.rdg2d = false; 
	partition_config.rdg3d = false; 
	partition_config.ba = false; 
	partition_config.rhg = false; 
	partition_config.kagen_d_rhg = 0;
	partition_config.kagen_gamma = 0; 

    // Stream Map
    partition_config.specify_alpha_multiplier = 1;
    partition_config.stream_multisection = false;
    partition_config.stream_modules_weight = NULL;
    partition_config.degree_nodeLayerModule = NULL;
    partition_config.stream_weighted_msec_type = WMSEC_COMMUNIC_COST;
    partition_config.onepass_pipelined_input = false;
    partition_config.onepass_simplified_input = false;
    partition_config.multicore_pipeline = false;
    partition_config.nodes_pipeline = NULL;
    partition_config.pipeline_stages = 0;
    partition_config.parallel_nodes = 1;
    partition_config.hashify_layers = 0;
    partition_config.fast_alg = ONEPASS_HASHING;
    partition_config.stream_rec_bisection = false;
    partition_config.stream_rec_bisection_base = 2;
    partition_config.stream_rec_biss_orig_alpha = false;
    partition_config.non_hashified_layers = std::numeric_limits<PartitionID>::max();
    partition_config.percent_non_hashified_layers = 1.0;
    partition_config.neighbor_blocks.resize(1);
    partition_config.all_blocks_to_keys.resize(1);
    partition_config.next_key.resize(1);
    partition_config.read_ew = false;
    partition_config.read_nw = false;

    partition_config.t1 = new double[1];
    partition_config.t2 = new double[1];
    partition_config.t3 = new double[1];
    (*partition_config.t1) = 0;
    (*partition_config.t2) = 0;
    (*partition_config.t3) = 0;


    // Initial partitioning via growning multiple BFS trees around artificial nodes
    partition_config.initial_part_multi_bfs = false;
    partition_config.multibfs_tries = 1;


    // Initial partitioning via Fennel on the coarsest level
    partition_config.initial_part_fennel_tries = 1;

}

inline void configuration::fastsocial(PartitionConfig &partition_config) {
    eco(partition_config);
    partition_config.balance_factor = 0;
}

inline void configuration::ecosocial(PartitionConfig &partition_config) {
    eco(partition_config);
    partition_config.balance_factor = 0.016;
}



inline void configuration::stream_map(PartitionConfig &partition_config) {
    fastsocial(partition_config);
    partition_config.stream_input = true;
    partition_config.stream_initial_bisections = true;
    partition_config.fennel_dynamics = FENNELADP_ORIGINAL;
    partition_config.batch_inbalance = 20;
    partition_config.use_fennel_objective = true;
}

inline void configuration::stream_partition(PartitionConfig &partition_config) {
    fastsocial(partition_config);
    partition_config.stream_input = true;
    partition_config.stream_initial_bisections = true;
    partition_config.fennel_dynamics = FENNELADP_ORIGINAL;
    partition_config.batch_inbalance = 20;
    partition_config.use_fennel_objective = true;
}

#endif /* end of include guard: CONFIGURATION_3APG5V7Z */
