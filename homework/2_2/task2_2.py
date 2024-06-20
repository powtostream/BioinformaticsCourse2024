import subprocess
from matplotlib import pyplot as plt

from Bio import SeqIO
import networkx as nx
from collections import defaultdict
import sys

test_fasta = "reads.fastq"


# class Node:
#
#     def __init__(self, seq=None, start=None, read=None):
#         self.seq = seq
#         self.start = start
#         self.read = read


def preprocess(fasta_path, trimmomatics=False, coral=False):
    if trimmomatics:
        subprocess.run(
            f"java -jar /usr/share/java/trimmomatic.jar SE -phred33 {fasta_path} filtered_reads.fastq SLIDINGWINDOW:4:15 LEADING:5 TRAILING:5 MINLEN:36")
        fasta_path = "filtered_reads.fastq"
    if coral:
        subprocess.run(
            f"/home/sebas/Documents/HSE/bioinf/coral-1.4.1/coral -fq {fasta_path} -o corrected.fastq -e 0.5")
        fasta_path = "corrected.fastq"
    reads = list(SeqIO.parse(fasta_path, "fastq"))
    return reads


def plot_graph(graph):
    pos = nx.spring_layout(graph)
    nx.draw(graph, pos)
    node_labels = {k:k for k, v in graph.nodes.items()}

    nx.draw_networkx_labels(graph, pos, labels=node_labels)
    # edge_seq = nx.get_edge_attributes(graph, 'seq')
    # edge_cov = nx.get_edge_attributes(graph, 'coverage')
    # edge_labels = {key: f"{edge_seq[key]}_{edge_cov[key]}" for key in edge_cov.keys()}
    edge_labels = nx.get_edge_attributes(graph, 'seq')
    nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels)
    plt.show()


def build_de_bruijn(fasta_path, k, trimmomatic=False, coral=False):
    reads = preprocess(fasta_path, trimmomatic, coral)
    graph = nx.DiGraph()
    edge_coverage = defaultdict(int)

    for num, read in enumerate(reads):
        print(num)
        for i in range(len(read) - k):
            kmer = str(read.seq[i:i + k])
            next_kmer = str(read.seq[i + 1:i + k + 1])
            graph.add_edge(kmer, next_kmer, seq=str(read.seq[i:i + k + 1]))
            edge_coverage[(kmer, next_kmer)] += 1

    nx.set_edge_attributes(graph, edge_coverage, 'coverage')
    nx.set_edge_attributes(graph, 1, 'condensed')
    return graph


def compress_graph(graph):
    # nx.set_edge_attributes(graph, 1, 'condensed')
    iter_list = list(graph.nodes())
    for node in iter_list:
        if graph.in_degree(node) == 1 and graph.out_degree(node) == 1:
            prev_node = next(graph.predecessors(node))
            next_node = next(graph.successors(node))
            print(prev_node, node, next_node)
            prev_edge = graph.edges[(prev_node, node)]
            next_edge = graph.edges[(node, next_node)]
            print(prev_edge, next_edge)
            condensed = prev_edge["condensed"]+next_edge["condensed"]
            coverage = (
                prev_edge["coverage"]*prev_edge["condensed"] +
                next_edge["coverage"]*next_edge["condensed"]
                ) / condensed
            seq = prev_edge["seq"]+next_edge["seq"][len(node):]
            graph.add_edge(
                prev_node, next_node, coverage=coverage, condensed=condensed,
                seq=seq
            )
            graph.remove_node(node)
    return graph

# gr = build_de_bruijn("corrected.fastq", k=80, trimmomatic=False, coral=False)
# nx.write_graphml(gr, "de_brujn")
# gr2 = nx.read_graphml("de_brujn")
# compressed_graph = compress_graph(gr2)
# nx.draw(compressed_graph)
# plt.show()


gr = build_de_bruijn("reads_light.fastq", k=6, trimmomatic=False, coral=False)
plot_graph(gr)
gr = compress_graph(gr)
plot_graph(gr)