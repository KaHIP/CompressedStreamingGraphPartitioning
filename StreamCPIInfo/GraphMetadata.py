# automatically generated by the FlatBuffers compiler, do not modify

# namespace: StreamCPIInfo

import flatbuffers
from flatbuffers.compat import import_numpy
np = import_numpy()

class GraphMetadata(object):
    __slots__ = ['_tab']

    @classmethod
    def GetRootAs(cls, buf, offset=0):
        n = flatbuffers.encode.Get(flatbuffers.packer.uoffset, buf, offset)
        x = GraphMetadata()
        x.Init(buf, n + offset)
        return x

    @classmethod
    def GetRootAsGraphMetadata(cls, buf, offset=0):
        """This method is deprecated. Please switch to GetRootAs."""
        return cls.GetRootAs(buf, offset)
    # GraphMetadata
    def Init(self, buf, pos):
        self._tab = flatbuffers.table.Table(buf, pos)

    # GraphMetadata
    def Filename(self):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(4))
        if o != 0:
            return self._tab.String(o + self._tab.Pos)
        return None

    # GraphMetadata
    def NumNodes(self):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(6))
        if o != 0:
            return self._tab.Get(flatbuffers.number_types.Uint64Flags, o + self._tab.Pos)
        return 0

    # GraphMetadata
    def NumEdges(self):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(8))
        if o != 0:
            return self._tab.Get(flatbuffers.number_types.Uint64Flags, o + self._tab.Pos)
        return 0

    # GraphMetadata
    def MaxDegree(self):
        o = flatbuffers.number_types.UOffsetTFlags.py_type(self._tab.Offset(10))
        if o != 0:
            return self._tab.Get(flatbuffers.number_types.Uint64Flags, o + self._tab.Pos)
        return 0

def GraphMetadataStart(builder):
    builder.StartObject(4)

def Start(builder):
    GraphMetadataStart(builder)

def GraphMetadataAddFilename(builder, filename):
    builder.PrependUOffsetTRelativeSlot(0, flatbuffers.number_types.UOffsetTFlags.py_type(filename), 0)

def AddFilename(builder, filename):
    GraphMetadataAddFilename(builder, filename)

def GraphMetadataAddNumNodes(builder, numNodes):
    builder.PrependUint64Slot(1, numNodes, 0)

def AddNumNodes(builder, numNodes):
    GraphMetadataAddNumNodes(builder, numNodes)

def GraphMetadataAddNumEdges(builder, numEdges):
    builder.PrependUint64Slot(2, numEdges, 0)

def AddNumEdges(builder, numEdges):
    GraphMetadataAddNumEdges(builder, numEdges)

def GraphMetadataAddMaxDegree(builder, maxDegree):
    builder.PrependUint64Slot(3, maxDegree, 0)

def AddMaxDegree(builder, maxDegree):
    GraphMetadataAddMaxDegree(builder, maxDegree)

def GraphMetadataEnd(builder):
    return builder.EndObject()

def End(builder):
    return GraphMetadataEnd(builder)