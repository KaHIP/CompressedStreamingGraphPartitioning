StreamCPI 1.00
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++](https://img.shields.io/badge/C++-20-blue.svg)](https://isocpp.org/)
[![CMake](https://img.shields.io/badge/CMake-3.24+-064F8C.svg)](https://cmake.org/)
[![Linux](https://img.shields.io/badge/Linux-supported-success.svg)](https://github.com/KaHIP/CompressedStreamingGraphPartitioning)
[![macOS](https://img.shields.io/badge/macOS-supported-success.svg)](https://github.com/KaHIP/CompressedStreamingGraphPartitioning)
[![GitHub Stars](https://img.shields.io/github/stars/KaHIP/CompressedStreamingGraphPartitioning)](https://github.com/KaHIP/CompressedStreamingGraphPartitioning/stargazers)
[![GitHub Issues](https://img.shields.io/github/issues/KaHIP/CompressedStreamingGraphPartitioning)](https://github.com/KaHIP/CompressedStreamingGraphPartitioning/issues)
[![Last Commit](https://img.shields.io/github/last-commit/KaHIP/CompressedStreamingGraphPartitioning)](https://github.com/KaHIP/CompressedStreamingGraphPartitioning/commits)
[![Homebrew](https://img.shields.io/badge/Homebrew-available-orange)](https://github.com/KaHIP/homebrew-kahip)
[![arXiv](https://img.shields.io/badge/arXiv-2410.07732-b31b1b.svg)](https://arxiv.org/abs/2410.07732)
[![ACDA'25](https://img.shields.io/badge/ACDA'25-published-blue)](https://arxiv.org/abs/2410.07732)
[![Heidelberg University](https://img.shields.io/badge/Heidelberg-University-c1002a)](https://www.uni-heidelberg.de)
=====

<p align="center">
  <img src="https://raw.githubusercontent.com/KaHIP/CompressedStreamingGraphPartitioning/main/logo/streamcpi-banner.png" alt="StreamCPI Banner" width="900"/>
</p>

**StreamCPI** is a framework for reducing the memory consumption of streaming graph partitioners by compressing block assignments using run-length encoding. Part of the [KaHIP](https://github.com/KaHIP) organization.

| | |
|:--|:--|
| **What it solves** | Memory-efficient streaming graph partitioning for trillion-edge graphs on edge devices |
| **Techniques** | Run-length compressed partition indices (CPI), modified Fennel scoring, batch-wise compression, external memory PQ |
| **Interfaces** | CLI |
| **Requires** | C++20, CMake 3.24+, MPI, OpenMP, Argtable |

## Quick Start

### Install via Homebrew

```bash
brew install KaHIP/kahip/streamcpi
```

### Or build from source

```bash
git clone --recursive https://github.com/KaHIP/CompressedStreamingGraphPartitioning.git
cd CompressedStreamingGraphPartitioning
./compile.sh
```

Alternatively, use the standard CMake build process:

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

The resulting binaries are `deploy/stream_cpi` and `deploy/stream_cpi_generated`.

### Run

```bash
# Partition a METIS graph into k blocks
./deploy/stream_cpi <graph> --k=<number of blocks>

# With full run-length compression (recommended)
./deploy/stream_cpi <graph> --k=<number of blocks> --rle_length=0

# With kappa scaling for further memory reduction
./deploy/stream_cpi <graph> --k=<number of blocks> --rle_length=0 --kappa=20

# Full parameter list
./deploy/stream_cpi --help
```

---

## Compression Modes

The `--rle_length` flag selects the compression mode:

| rle_length  | Mode                                                                                     |
|-------------|------------------------------------------------------------------------------------------|
| 0           | Complete run-length compression (recommended)                                            |
| -1          | std::vector (fastest, no compression)                                                    |
| -2          | External memory PQ (using STXXL, configurable in `lib/data_structure/ExternalPQ.h`)      |
| 100+        | Batch-wise compression: each compression vector handles `rle_length` nodes               |

## CPI Compression Vector

The (semi-)dynamic compression vector can be used as a drop-in replacement for `std::vector` in any streaming algorithm that stores arrays with repeating values. The standalone library is available at [kurpicz/cpi](https://github.com/kurpicz/cpi).

## Streaming Graph Generator

The included `stream_cpi_generated` binary partitions graphs generated on-the-fly using a streaming graph generator:

```bash
# Barabasi-Albert graph
./deploy/stream_cpi_generated <output> --k=<blocks> --rle_length=0 --kappa=20 \
    --ba --nodes_to_generate=<n> --kagen_d_ba=<avg_degree> --kagen_chunk_count=<chunks>

# RGG2D graph
./deploy/stream_cpi_generated <output> --k=<blocks> --rle_length=0 --kappa=20 \
    --rgg2d --nodes_to_generate=<n> --kagen_r=<radius> --kagen_chunk_count=<chunks>
```

See [adilchhabra/KaGen](https://github.com/adilchhabra/KaGen) for graph generation models and parameters.

---

## Notes

- Results are stored as [FlatBuffer](https://github.com/google/flatbuffers) `.bin` files when passing `--write_results`.
- 64-bit vertex IDs are enabled by default. To disable, set `64BITVERTEXMODE` to `OFF` in `CMakeLists.txt`.
- For the METIS graph format, refer to the [KaHIP manual](https://github.com/KaHIP/KaHIP/raw/master/manual/kahip.pdf).

## Data References

Graphs used in our experiments were sourced from:
- [SNAP Dataset](https://snap.stanford.edu/data/)
- [10th DIMACS Challenge](https://sites.cc.gatech.edu/dimacs10/downloads.shtml)
- [Laboratory for Web Algorithmics](https://law.di.unimi.it/datasets.php)
- [Network Repository](https://networkrepository.com/)

---

## Citing

If you use StreamCPI in your research, please cite:

```bibtex
@inproceedings{chhabra2025streamcpi,
    title     = {Partitioning Trillion Edge Graphs on Edge Devices},
    author    = {Adil Chhabra and Florian Kurpicz and Christian Schulz and Dominik Schweisgut and Daniel Seemaier},
    booktitle = {SIAM Conference on Applied and Computational Discrete Algorithms (ACDA)},
    year      = {2025}
}
```

## Licensing

StreamCPI is distributed under the MIT License. See [LICENSE](LICENSE) for details.
