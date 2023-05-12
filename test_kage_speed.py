#from kivs import Graph, KmerFinder, hash_kmer
from obgraph import Graph as OBGraph
#from graph_kmer_index import kmer_hash_to_sequence
from graph_kmer_index.kmer_finder import DenseKmerFinder
from timeit import timeit
#import numpy as np
#import npstructures as nps

def measure_time(obgraph, k, max_variant_nodes):
        finder = DenseKmerFinder(obgraph, k=k, max_variant_nodes=max_variant_nodes)
        finder.find()

obgraph = OBGraph.from_file("data/small_graph.npz")
print(len(obgraph.nodes))
finder = DenseKmerFinder(obgraph, k=31, max_variant_nodes=31)
finder.find()

for k in [3, 15, 31]:
    duration = timeit(lambda: measure_time(obgraph, k, k), number=5)
    print("KAGE solution used", duration, "for k =", k)

