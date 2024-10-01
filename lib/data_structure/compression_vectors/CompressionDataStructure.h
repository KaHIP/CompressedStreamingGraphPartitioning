//
// Created by adilchhabra on 16.05.24.
//

#pragma once
#include "definitions.h"

template <typename T, bool enable_buffer = false>
class CompressionDataStructure {
public:
  virtual void Append(T value) = 0;
  virtual void BatchAppend(NodeID batch_index, T value) = 0;
  virtual T GetValueByIndex(NodeID index) const = 0;
  virtual T GetValueByBatchIndex(NodeID batch_index, NodeID data_index) const = 0;
};
