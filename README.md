# KIVS: Graph **K**-mer **I**ndexer and **V**ariant **S**ignature Finder

Benchmarking for KIVS, a master's thesis project about discovering and indexing k-mers from genome graphs.

## Description

This project aims to create a high-performance cython-based python module for creating genome graphs and discovering/indexing the k-mers within the graph.

For requirements and setup, see the [KIVS repository](https://github.com/ZinderAsh/python-kivs).
Some tests require the setup if [kage-indexing](https://github.com/ivargr/kage-indexing) as well.

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

```

### Variant Signature Performance and Accuracy


