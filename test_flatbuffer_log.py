import sys
import flatbuffers
from StreamCPIInfo import Partition, GraphMetadata, PartitionConfiguration, RunTime, PartitionMetrics, MemoryConsumption

def read_partition_from_file(file_path):
    with open(file_path, 'rb') as file:
        # Read the FlatBuffer binary file
        buf = file.read()

        # Create a FlatBuffer object
        partition = Partition.Partition.GetRootAsPartition(buf, 0)

        # Access and print the fields
        graph_metadata = partition.GraphMetadata()
        print("Graph Metadata:")
        print(f"  Filename: {graph_metadata.Filename().decode('utf-8')}")
        print(f"  Num Nodes: {graph_metadata.NumNodes()}")
        print(f"  Num Edges: {graph_metadata.NumEdges()}")
        print(f"  Max Degree: {graph_metadata.MaxDegree()}")

        partition_config = partition.PartitionConfiguration()
        print("\nPartition Configuration:")
        print(f"  k: {partition_config.K()}")
        print(f"  Seed: {partition_config.Seed()}")
        print(f"  Input Balance: {partition_config.InputBalance()}")
        if(partition_config.SetPartZero() == 1) :
            print(f"  Part ID: {0}")
        else :
            print(f"  Part ID: {'Fennel'}")
        if(partition_config.RleLength() == 18446744073709551615) :
            print(f"  RLE Length: {None}")
        else :
            print(f"  Part ID: {partition_config.RleLength()}")

        runtime = partition.Runtime()
        print("\nRun Time:")
        print(f"  IO Time: {runtime.IoTime()}")
        print(f"  Mapping Time: {runtime.MappingTime()}")
        print(f"  Total Time: {runtime.TotalTime()}")

        metrics = partition.Metrics()
        print("\nPartition Metrics:")
        print(f"  Edge Cut: {metrics.EdgeCut()}")
        print(f"  Balance: {metrics.Balance()}")

        memory_consumption = partition.MemoryConsumption()
        print("\nMemory Consumption:")
        print(f"  Uncompressed Vec Bytes: {memory_consumption.UncompressedVecBytes()}")
        print(f"  RLE Vector Bytes: {memory_consumption.RleVectorBytes()}")
        print(f"  RLE Vector MB: {memory_consumption.RleVectorMb()}")
        print(f"  Relative Compression: {memory_consumption.RelativeCompression()}")
        print(f"  Overall Max RSS: {memory_consumption.OverallMaxRss()}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <binary_file_path>")
        sys.exit(1)

    file_path = sys.argv[1]
    read_partition_from_file(file_path)
