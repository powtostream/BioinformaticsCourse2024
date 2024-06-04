import random
from datetime import datetime

import pandas as pd
import json
import subprocess


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

with open("data.json", "r") as file:
   data = json.load(file)

subprocess.run("java -jar /usr/share/java/trimmomatic.jar SE -phred33 reads.fastq filtered_reads.fastq SLIDINGWINDOW:4:15 LEADING:5 TRAILING:5 MINLEN:36")
"""
TrimmomaticSE: Started with arguments:
 -phred33 reads.fastq filtered_reads.fastq SLIDINGWINDOW:4:15 LEADING:5 TRAILING:5 MINLEN:36 MAXINFO:40:0.5
Automatically using 4 threads
Input Reads: 50000 Surviving: 0 (0.00%) Dropped: 50000 (100.00%)
TrimmomaticSE: Completed successfully
"""
# trimmomatics удалил все мои риды
"""
sebas@pc1:~/Documents/HSE/bioinf/coral-1.4.1$ ./coral -fq ../BioinformaticsCourse2024/homework/2_1/reads.fastq -o ../BioinformaticsCourse2024/homework/2_1/corrected.fastq -e 0.5
Counting reads in file ../BioinformaticsCourse2024/homework/2_1/reads.fastq
Found 50000 reads
Reading reads from file ../BioinformaticsCourse2024/homework/2_1/reads.fastq
Found 50000 reads
Length of k-mers: 21
Counting k-mers
0 reads, 0 different k-mers
4957826 different k-mers, 11475312 total k-mers
Constructing read lists
0 reads
Q-gram index built.
Correcting reads
Thread 0/8 starting
Thread 6/8 starting
Thread 6/8 finishing
Thread 1/8 starting
Thread 5/8 starting
Thread 2/8 starting
Thread 7/8 starting
Thread 3/8 starting
Thread 4/8 starting
Thread 5/8 finishing
Thread 7/8 finishing
Thread 3: 39999/50000 Perfect alignments: 238945, Good alignments: 1427, Bad alignments: 31951
Thread 3/8 finishing
Thread 4: 49999/50000 Perfect alignments: 240021, Good alignments: 1429, Bad alignments: 32087
Thread 4/8 finishing
Thread 1: 19999/50000 Perfect alignments: 240279, Good alignments: 1433, Bad alignments: 32112
Thread 1/8 finishing
Thread 2: 29999/50000 Perfect alignments: 240280, Good alignments: 1433, Bad alignments: 32112
Thread 2/8 finishing
Thread 0: 9999/50000 Perfect alignments: 240603, Good alignments: 1433, Bad alignments: 32158
Thread 0/8 finishing
Perfect alignments: 240603, Good alignments: 1433, Bad alignments: 32158
50000 reads were aligned, 50000 reads were aligned in a good alignment
Quick alignments: 60723285, Banded alignments: 58049, Full alignments: 0
Correction finished. Outputting corrected reads.
"""
true_seq = data["true"]
corrected_reads = list(SeqIO.parse("corrected.fastq", "fastq"))
init_reads = list(SeqIO.parse("reads.fastq", "fastq"))

errors = data["errors"]
errors = {key:val for key, val in enumerate(errors)}
correction_df = pd.DateaFrame(columns=["read#", "pos", "init", "true", "corrected"])
# for read_err in errors:
#    for err in read_err:
#       try:
#          true_read = true_seq[err[2]]
#          init_read = init_reads[err[0]][err[1]]
#          corrected_read = corrected_reads[err[0]][err[1]]
#          correction_df.loc[len(correction_df)] = {
#             "read#": err[0],
#             "pos": err[1],
#             "init": init_read,
#             "true": true_read,
#             "corrected": corrected_read
#          }
#       except IndexError:
#          continue

tp = 0
tn = 0
fp = 0
fn = 0
start = datetime.utcnow()
for i, raw_read in enumerate(init_reads):
    print(i, "left:", (datetime.utcnow() - start).seconds/(i+1)*(50000-i))
    corrected_read = corrected_reads[i]
    read_errors = {read_ind: whole_ind for _, read_ind, whole_ind in errors[i]}
    for pos, raw_nucl in enumerate(raw_read):
        try:
            er_pos = read_errors.get(pos)
            corrected = corrected_read[pos]
            if er_pos:
                true = true_seq[er_pos]
                if true == corrected:
                    tp += 1
                else:
                    fn += 1
            else:
                if raw_nucl == corrected:
                    tn += 1
                else:
                    fp += 1

        except IndexError:
            continue
print("tp:", tp, "fp:", fp, "tn:", tn, "fn:", fn)
acc = (tp+tn)/(tp+tn+fp+fn)
recall = tp/(tp+fn)
prec = tp/(tp+fp)
print("acc:", acc, "recall:", recall, "precision:", prec)
"""
tp: 525752 fp: 143515 tn: 11798695 fn: 7249
acc: 0.988 recall: 0.986 precision: 0.786
Вывод у инструмента coral высокая полнота, большинство ошибок исправлено
Однако при этом не слишком высокая точность, порядка 20 процентов из 
исправленных позиций исправлены ошибочно. Возможно подбор параметров 
может улучшить качество исправления.
"""