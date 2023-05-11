from kivs import Graph, KmerFinder
from obgraph import Graph as OBGraph
import sys

if len(sys.argv) != 2:
    print("Usage:")
    print("python test_kmer_index_speed.py filepath.gfa")
    exit(0)

graph = Graph.from_file(sys.argv[1])
kf = KmerFinder(graph, 31)
kf.find(stdout=True, include_spanning_nodes=False)
