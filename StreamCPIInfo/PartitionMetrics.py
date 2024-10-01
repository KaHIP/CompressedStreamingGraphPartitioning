# automatically generated by the FlatBuffers compiler, do not modify

# namespace: StreamCPIInfo

import flatbuffers
from flatbuffers.compat import import_numpy
np = import_numpy()

class PartitionMetrics(object):
    __slots__ = ['_tab']

    @classmethod
    def GetRootAs(cls, buf, offset=0):
        n = flatbuffers.encode.Get(flatbuffers.packer.uoffset, buf, offset)
        x = PartitionMetrics()
        x.Init(buf, n + offset)
        return x

    @classmethod
    def GetRootAsPartitionMetrics(cls, buf, offset=0):
        """This method is deprecated. Please switch to GetRootAs."""
        return cls.GetRootAs(buf, offset)
    # PartitionMetrics
    def Init(self, buf, pos):
        self._tab = flatbuffers.table.Table(buf, pos)

    # PartitionMetrics
    def EdgeCut(self):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(4))
        if o != 0:
            return self._tab.Get(flatbuffers.number_types.Uint64Flags, o + self._tab.Pos)
        return 0

    # PartitionMetrics
    def Balance(self):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(6))
        if o != 0:
            return self._tab.Get(flatbuffers.number_types.Float64Flags, o + self._tab.Pos)
        return 0.0

def PartitionMetricsStart(builder):
    builder.StartObject(2)

def Start(builder):
    PartitionMetricsStart(builder)

def PartitionMetricsAddEdgeCut(builder, edgeCut):
    builder.PrependUint64Slot(0, edgeCut, 0)

def AddEdgeCut(builder, edgeCut):
    PartitionMetricsAddEdgeCut(builder, edgeCut)

def PartitionMetricsAddBalance(builder, balance):
    builder.PrependFloat64Slot(1, balance, 0.0)

def AddBalance(builder, balance):
    PartitionMetricsAddBalance(builder, balance)

def PartitionMetricsEnd(builder):
    return builder.EndObject()

def End(builder):
    return PartitionMetricsEnd(builder)