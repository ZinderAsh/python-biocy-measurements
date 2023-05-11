from biocy import Graph, KmerFinder
import sys

if len(sys.argv) != 3:
    print("Usage:")
    print("python gfa_to_bcg.py filepath.gfa filepath.bcg")
    exit(1)

graph = Graph.from_gfa(sys.argv[1])
graph.to_file(sys.argv[2])

exit(0)
