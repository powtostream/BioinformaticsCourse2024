import subprocess
from datetime import datetime

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
    print(datetime.utcnow())
    # reads = list(SeqIO.parse(fasta_path, "fastq"))
    reads = []
    with open("corrected.fasta", "r") as file:
        for i, line in enumerate(file.readlines()):
            if i % 4 == 1:
                reads.append(line.strip().upper())
    print(reads[-1])
    print(datetime.utcnow())
    return reads


def plot_graph(graph, coverage=False, nodes=False):
    pos = nx.spring_layout(graph)
    nx.draw(graph, pos)
    if nodes:
        node_labels = {k:k for k, v in graph.nodes.items()}
        nx.draw_networkx_labels(graph, pos, labels=node_labels)
    if coverage:
        edge_seq = nx.get_edge_attributes(graph, 'seq')
        edge_cov = nx.get_edge_attributes(graph, 'coverage')
        edge_labels = {key: f"{edge_seq[key]}_{edge_cov[key]:.1f}" for key in edge_cov.keys()}
    else:
        edge_labels = nx.get_edge_attributes(graph, 'seq')
    nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels)
    plt.show()


def build_de_bruijn(fasta_path, k, trimmomatic=False, coral=False):
    reads = preprocess(fasta_path, trimmomatic, coral)
    graph = nx.DiGraph()
    edge_coverage = defaultdict(int)

    start = datetime.utcnow()
    for num, read in enumerate(reads):
        if num % 100 == 1:
            now = datetime.utcnow()
            print(f"{num}/{len(reads)}", now, "left:", (now-start).seconds/60/num*(len(reads)-num))
        for i in range(len(read) - k):
            # kmer = str(read.seq[i:i + k])
            # next_kmer = str(read.seq[i + 1:i + k + 1])
            # graph.add_edge(kmer, next_kmer, seq=str(read.seq[i:i + k + 1]))
            # edge_coverage[(kmer, next_kmer)] += 1
            kmer = str(read[i:i + k])
            next_kmer = str(read[i + 1:i + k + 1])
            graph.add_edge(kmer, next_kmer, seq=str(read[i:i + k + 1]))
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
            # break
    return graph


def cut_tails(graph):
    iter_list = list(graph.nodes())
    tails = []
    for node in iter_list:
        if graph.in_degree(node) == 0 and graph.out_degree(node) == 1:
            tail_start = node
            tail_list = [tail_start]
            next_node = next(graph.successors(node))
            tail_len = 1
            tail_cov = graph.edges[(node, next_node)]["coverage"]
            while graph.in_degree(next_node) == 1 and graph.out_degree(next_node) == 1:
                tail_list.append(next_node)
                node, next_node = next_node, next(graph.successors(node))
                tail_len += 1
                tail_cov *= graph.edges[(node, next_node)]["coverage"]
            tails.append((tail_cov/tail_len, tail_start, "start", tail_list))

        elif graph.in_degree(node) == 1 and graph.out_degree(node) == 0:
            tail_end = node
            tail_list = [tail_end]
            next_node = next(graph.predecessors(node))
            tail_len = 1
            tail_cov = graph.edges[(next_node, node)]["coverage"]
            while graph.in_degree(next_node) == 1 and graph.out_degree(next_node) == 1:
                tail_list.append(next_node)
                node, next_node = next_node, next(graph.predecessors(node))
                tail_len += 1
                tail_cov *= graph.edges[(next_node, node)]["coverage"]
            tails.append((tail_cov/tail_len, tail_end, "end", tail_list))

    tails_sorted = sorted(tails, key=lambda x: x[0])

    threshold = int(len(tails_sorted)*0.7)+1
    for _, tail_node, tail_type, lll in tails_sorted[:int(len(tails_sorted)*0.3)+1]:
    # tail_num = 0
    # while tail_num <= threshold:
    #     _, tail_node, tail_type, lll = tails_sorted[tail_num]
    #     tail_num += 1
        try:
            if tail_type == "start":
                node = next(graph.successors(tail_node))
                graph.remove_edge(tail_node, node)
                while graph.out_degree(node) == 1 and graph.out_degree(node) == 1:
                    next_node = next(graph.successors(node))
                    # graph.remove_node(node)
                    graph.remove_edge(node, next_node)
                    node = next_node

            if tail_type == "end":
                node = next(graph.predecessors(tail_node))
                # graph.remove_node(tail_node)
                graph.remove_edge(node, tail_node)
                while graph.out_degree(node) == 1 and graph.out_degree(node) == 1:
                    next_node = next(graph.predecessors(node))
                    # graph.remove_node(node)
                    graph.remove_edge(next_node, node)
                    node = next_node
        except KeyError:
            # graph.remove_node(tail_node)
            # threshold += 1
            pass
        except StopIteration:
            # graph.remove_node(tail_node)
            # threshold += 1
            pass
    iter_list = list(graph.nodes())
    for node in iter_list:
        if graph.in_degree(node) == 0 and graph.out_degree(node) == 0:
            graph.remove_node(node)
    return graph

# gr = build_de_bruijn("corrected.fastq", k=80, trimmomatic=False, coral=False)
# nx.write_graphml(gr, "de_brujn")
# gr2 = nx.read_graphml("de_brujn")
# compressed_graph = compress_graph(gr2)
# nx.draw(compressed_graph)
# plt.show()


gr = build_de_bruijn("corrected.fasta", k=48, trimmomatic=False, coral=False)
nx.write_graphml(gr, "de_brujn_init")
# plot_graph(gr, True, True)
# for i in range(20):
#     gr = compress_graph(gr)
#     plot_graph(gr, True, True)
gr = compress_graph(gr)
# plot_graph(gr, True, True)
gr = cut_tails(gr)
nx.write_graphml(gr, "de_brujn_compressed_cut")
# plot_graph(gr, True, True)