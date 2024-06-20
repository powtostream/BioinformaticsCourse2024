from datetime import datetime

from Bio import SeqIO

print(datetime.utcnow())
# reads = list(SeqIO.parse("corrected.fasta", "fastq"))
reads = []
with open("corrected.fasta", "r") as file:
    for i, line in enumerate(file.readlines()):
        if i % 4 == 1:
            reads.append(line.strip().upper())
print(len(reads))
print(reads[7], reads[9], reads[0], reads[123123])
print(datetime.utcnow())
