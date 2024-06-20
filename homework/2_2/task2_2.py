import subprocess

from Bio import SeqIO
import networkx as nx
from collections import defaultdict
import sys

test_fasta = "reads.fastq"


class Node:

    def __init__(self, seq=None, start=None, read=None):
        self.seq = seq
        self.start = start
        self.read = read


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


def build_de_bruijn(fasta_path, k, trimmomatic=False, coral=False):
    reads = preprocess(fasta_path, trimmomatic, coral)
    nodes = dict()

    edge_coverage = defaultdict(int)

    for num, read in enumerate(reads):
        print(num)
        for i in range(len(read) - k - 2):
            kmer = str(read.seq[i:i + k])
            if nodes.get(kmer) is None:
                kmer_node = Node(seq=kmer, start=i, read=num)
                nodes[kmer] = kmer_node
            next_kmer = str(read.seq[i + 1:i + k + 1])
            if nodes.get(next_kmer) is None:
                next_kmer_node = Node(seq=next_kmer, start=i+1, read=num)
                nodes[next_kmer] = next_kmer_node
            edge_coverage[(kmer, next_kmer)] += 1
    return nodes, edge_coverage, reads


def compress_graph(graph):
    for node in graph.nodes():
        pass


nodes, cov, reads = build_de_bruijn("corrected.fastq", k=80, trimmomatic=False, coral=False)
print("END")
