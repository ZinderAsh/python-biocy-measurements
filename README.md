# KIVS: Graph **K**-mer **I**ndexer and **V**ariant **S**ignature Finder

Benchmarking for KIVS, a master's thesis project about discovering and indexing k-mers from genome graphs.

## Description

This project aims to create a high-performance cython-based python module for creating genome graphs and discovering/indexing the k-mers within the graph.

## Requirements

For requirements and setup, see the [KIVS repository](https://github.com/ZinderAsh/python-kivs).
Some tests require the setup if [kage-indexing](https://github.com/ivargr/kage-indexing) as well.
The `time` command is required for some tests.

## Reproducing Thesis Benchmarks

### Python Prototype Performance

For runtime measurements of the Python prototypes, run:
```bash
python test_python_speed.py
```
Results are in milliseconds.

### Compare Single Element Operations in Python Lists and NumPy Arrays

To get runtimes for setting the value of 100 million individual elements with a standard for-loop in Python lists and NumPy arrays, run:
```bash
python python_vs_numpy.py
```
Results are in milliseconds.

### Compare Map Encoding to Bit-Operation Encoding

Encoding is here referred to as "hashing". `cd` into the `hashing_speed` directory for these tests.

To get runtime measurements of encoding 100 million 31-mers with map-based encoding and bit-operation encoding, run:
```bash
make test-bit-hash
```

### Compare Encoding on Demand to Encoding Ahead of Time

These tests are also in the `hashing_speed` directory.

To get runtime measurements run:
```bash
make test-graph-hash
```

### K-mer Finding Performance

These tests include the performance of both [vg](https://github.com/vgteam/vg) and [odgi](https://github.com/pangenome/odgi). Executables for these need to be downloaded and put into the root directory of this project, as files named `vg` and `odgi`. To run the tests do:
```bash
# All tests
make test

# Single tests
make test-vg
make test-odgi
make test-kivs
make test-kivs-stdout
make test-kivs-full
```

### Variant Signature Performance and Accuracy

These tests require [kage-indexing](https://github.com/ivargr/kage-indexing) to be properly set up, and for KIVS to be installed in its conda environment. Further, the yeast dataset needs to be put into the `local_data` directory in kage-indexing. These can be found at [zenodo](https://zenodo.org/record/7929047). Mind that the yeast tests require a large amount of available memory and take a good while to complete.
For runtime tests, run these commands in the kage-indexing directory:
```bash
# KAGE implementation
snakemake --use-conda --config max_variant_nodes=3 k_threads_kmer_index=1 --cores 1 --forceall test_yeast_full
snakemake --use-conda --config max_variant_nodes=3 k_threads_kmer_index=16 --cores 1 --forceall test_yeast_full
snakemake --use-conda --config max_variant_nodes=3 k_threads_kmer_index=1 --cores 1 --forceall test_yeast_full
snakemake --use-conda --config max_variant_nodes=10 k_threads_kmer_index=16 --cores 1 --forceall test_yeast_full

# KIVS implementation
snakemake --use-conda --config max_variant_nodes=3 use_kivs=True kivs_minimize_overlaps=True kivs_align_windows=True --cores 1 --forceall test_yeast_full
snakemake --use-conda --config max_variant_nodes=10 use_kivs=True kivs_minimize_overlaps=True kivs_align_windows=True --cores 1 --forceall test_yeast_full
```
After each command, the relevant runtime can be found in `data/yeast_whole_genome/benchmarks/get_variant_kmers.tsv`

For accuracy tests, first copy `kivs_accuracy_test.sh` and `kivs_accuracy_yeast.sh` from python-kivs-benchmarking to the kage-indexing root directory. These automate the process for accuracy tests for all rows included in the accuracy tables. To run the tests, then do:
```bash
# For a subset of human chromosome 1
./kivs_accuracy_test.sh

# For the whole yeast genome
./kivs_accuracy_yeast.sh
```
