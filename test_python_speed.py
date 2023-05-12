#from kivs import Graph, KmerFinder, hash_kmer
from python_v1.Graph import Graph as Graph_v1
from python_v2.Graph import Graph as Graph_v2
from python_v3.Graph import Graph as Graph_v3
from python_v4.Graph import Graph as Graph_v4
from python_v5.Graph import Graph as Graph_v5
from python_v6.Graph import Graph as Graph_v6
from python_v7.Graph import Graph as Graph_v7
from obgraph import Graph as OBGraph
#from graph_kmer_index import kmer_hash_to_sequence
from graph_kmer_index.kmer_finder import DenseKmerFinder
from timeit import timeit
#import numpy as np
#import npstructures as nps

def measure_time_for_class(obgraph, GraphClass, k):
    graph = GraphClass.from_obgraph(obgraph)
    #graph.save_results = False
    durations = []
    runs = 7
    for i in range(runs):
        print("Run", i + 1, "of", runs)
        duration = timeit(lambda: graph.create_kmer_index(k), number=1)
        durations.append(round(duration * 1000))
    return durations

graph_classes = [
        (Graph_v1, "Initial Prototype"),
        (Graph_v2, "Numpy Arrays"),
        #(Graph_v3, "Replace Break with Return"),
        #(Graph_v4, "Added Reference Property"),
        (Graph_v5, "Numpy Array Slicing"),
        #(Graph_v6, "Empty Node Support"),
        (Graph_v7, "2-Bit Encode Graph")
]

obgraph = OBGraph.from_file("data/smaller_graph.npz")
for i in range(len(graph_classes)):
    durations = measure_time_for_class(obgraph, graph_classes[i][0], 31)
    print(graph_classes[i][1], sorted(durations))

exit(0)

def measure_time(obgraph, k, max_variant_nodes):
        finder = DenseKmerFinder(obgraph, k=k, max_variant_nodes=max_variant_nodes)
        finder.find()

obgraph = OBGraph.from_file("data/example_graph2.npz")
print(len(obgraph.nodes))
finder = DenseKmerFinder(obgraph, k=31, max_variant_nodes=31)
finder.find()

for k in [3, 15, 31]:
    duration = timeit(lambda: measure_time(obgraph, k, k), number=5)
    print("KAGE solution used", duration, "for k =", k)

