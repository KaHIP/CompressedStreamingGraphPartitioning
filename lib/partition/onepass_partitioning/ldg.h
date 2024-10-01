/******************************************************************************
 * ldg.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj
 *****************************************************************************/

#ifndef ONEPASS_LDG_7I4IR31Y
#define ONEPASS_LDG_7I4IR31Y

#include <algorithm>

#include "definitions.h"
#include "random_functions.h"
#include "timer.h"

#include "vertex_partitioning.h"

class onepass_ldg : public vertex_partitioning {
public:
  onepass_ldg(PartitionID k0, PartitionID kf, PartitionID max_blocks,
              NodeID n_threads, bool hashing = false);
  virtual ~onepass_ldg();
  void instantiate_blocks(LongNodeID n, LongEdgeID m, PartitionID k,
                          ImbalanceType epsilon);

protected:
  float compute_score(floating_block &block, int my_thread);
};

#endif /* end of include guard: ONEPASS_LDG_7I4IR31Y */
