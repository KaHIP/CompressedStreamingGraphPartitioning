/******************************************************************************
 * fennel.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj
 *****************************************************************************/

#ifndef ONEPASS_FENNEL_7I4IR31Y
#define ONEPASS_FENNEL_7I4IR31Y

#include <algorithm>

#include "definitions.h"
#include "random_functions.h"
#include "timer.h"

#include "vertex_partitioning.h"

class onepass_fennel : public vertex_partitioning {
public:
  onepass_fennel(PartitionID k0, PartitionID kf, PartitionID max_blocks,
                 NodeID n_threads, bool hashing = false, float gamma = 1.5);
  virtual ~onepass_fennel();
  void instantiate_blocks(LongNodeID n, LongEdgeID m, PartitionID k,
                          ImbalanceType epsilon);

protected:
  float compute_score(floating_block &block, int my_thread);
  float gamma;
};

#endif /* end of include guard: ONEPASS_FENNEL_7I4IR31Y */
