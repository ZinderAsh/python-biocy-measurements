from kivs import Graph, KmerFinder
from obgraph import Graph as OBGraph
import sys

if len(sys.argv) < 2:
    print("Usage:")
    print("python test_kmer_index_speed.py filepath.kivs")
    exit(0)

graph = Graph.from_file(sys.argv[1])
kf = KmerFinder(graph, 31)

stdout = False
include_spanning_nodes = False

if len(sys.argv) == 3:
    if sys.argv[2] == "stdout":
        stdout = True
    elif sys.argv[2] == "full":
        include_spanning_nodes = True

kf.find(stdout=stdout, include_spanning_nodes=include_spanning_nodes)
