import random
import numpy as np
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

alphabet = ["A", "T", "C", "G"]
alph_num = list(range(4))
al_dict = {nucl: i for i, nucl in enumerate(alphabet)}
alph_set = set(alphabet)
true_len = 20
true_string = "".join(random.choices(alphabet, k=true_len))
# print(true_string)


def get_quality(occ):
    q = -10*np.log10(1-occ/100)
    # c = chr(33+int(q))
    return int(q)



read_len_mean = 10
read_len_std = 0
read_num = 60
reads = []
reads_quality = []
errors = []
err_probs = [1/3 for _ in range(3)]
reads_start_pos = dict()
for i in range(read_num):
    read_len = int(random.gauss(read_len_mean, read_len_std))
    start = random.randint(0, true_len-read_len)
    read = []
    read_quality = []
    read_errors = []
    reads_start_pos[i] = start
    for j in range(start, start+read_len):
        # err = random.randint(0, 75)
        err = 0
        # true_prob = (100 - err)/100
        # err_prob = (err / 3)/100

        true = true_string[j]
        if err == 0:
            read.append(true)
            read_quality.append(93)
            continue
        err_list = [al_dict[n] for n in alphabet if n != true]
        # probs = [true_prob if nucl == true else err_prob for nucl in alphabet]
        choices = np.random.choice(err_list, size=err, p=err_probs)
        counts = np.bincount(choices, minlength=3)
        # cur_read_counts.append(counts)
        max_count = np.argmax(counts)
        if counts[max_count] > 100 - err:
            read_errors.append((i, j-start, j))
            read.append(alphabet[max_count])
            q = get_quality(counts[max_count])
            read_quality.append(q)
        else:
            read.append(true)
            q = get_quality(100-err)
            read_quality.append(q)
    reads.append("".join(read))
    reads_quality.append(read_quality)
    errors.append(read_errors)
    print(i)




records = []
for i, seq in enumerate(reads):
    record = SeqRecord(
        Seq(seq),
        id=f"read{i+1}",
        description="",
        letter_annotations={"phred_quality": reads_quality[i]}
    )
    records.append(record)


with open("reads_light.fastq", "w") as output_handle:
    SeqIO.write(records, output_handle, "fastq")

data_json = {"errors": errors, "true": true_string, "reads_start_pos": reads_start_pos}
with open("data_light.json", "w") as file:
    json.dump(data_json, file, indent=4)
