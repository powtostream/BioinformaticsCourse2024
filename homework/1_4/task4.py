import random

import numpy as np
from Bio import SeqIO
# from Bio.Align import substitution_matrices, PairwiseAligner

is_path = "data/islands.fasta"
nonis_path = "data/nonIslands.fasta"
islands = []
islands_fasta = SeqIO.parse(is_path, "fasta")
for record in islands_fasta.records:
    islands.append(str(record.seq))

nonislands_fasta = SeqIO.parse(nonis_path, "fasta")
nonislands = []
for record in nonislands_fasta.records:
    nonislands.append(str(record.seq))

nucleotides = ["A", "C", "G", "T"]
islands_one = {k: 0 for k in nucleotides}
non_islands_one = {k: 0 for k in nucleotides}

pairs = [a+b for a in nucleotides for b in nucleotides]
islands_pairs = {k: 0 for k in pairs}
non_islands_pairs = {k: 0 for k in pairs}


for seq in islands:
    for i in range(len(seq)-1):
        islands_pairs[seq[i]+seq[i+1]] += 1
        islands_one[seq[i]] += 1
    islands_one[seq[-1]] += 1


for seq in nonislands:
    for i in range(len(seq)-1):
        non_islands_pairs[seq[i]+seq[i+1]] += 1
        non_islands_one[seq[i]] += 1
    non_islands_one[seq[-1]] += 1

total_islands_ones = sum(islands_one.values())
total_islands_pairs = total_islands_ones - len(islands)
total_non_islands_ones = sum(non_islands_one.values())
total_non_islands_pairs = total_non_islands_ones - len(nonislands)


for n in nucleotides:
    islands_one[n] /= total_islands_ones
    non_islands_one[n] /= total_non_islands_ones

for p in pairs:
    islands_pairs[p] /= total_islands_pairs
    non_islands_pairs[p] /= total_non_islands_pairs

random.seed(333)


class Nucleotides_generator:
    cur_state = "-"
    # states = ["-", "+"]
    stay_probs = {"-": (("-", "+"), (0.995, 0.005)),
                  "+": (("+", "-"), (0.95, 0.05))}
    nucl_probs = {"-": list(non_islands_one.values()),
                  "+": list(islands_one.values())}
    nucleotides = list(non_islands_one.keys())
    print(nucleotides)
    print(nucl_probs)

    def __next__(self):
        cur_data = self.stay_probs.get(self.cur_state)
        self.cur_state = random.choices(cur_data[0], weights=cur_data[1], k=1)[
            0]

        nucleotid = random.choices(self.nucleotides,
                                   weights=self.nucl_probs.get(self.cur_state),
                                   k=1)[0]
        return self.cur_state, nucleotid

seqs = list()
states = list()
# for length in range(1000, 100001, 5000):
for length in range(1000, 1001):
    new_seq = list()
    new_states = list()
    generator = Nucleotides_generator()
    for i in range(length):
        state, nucl = next(generator)
        new_states.append(state)
        new_seq.append(nucl)
    seqs.append(new_seq)
    states.append(new_states)

transfer_probs = {0: (0.995, 0.05), 1: (0.005, 0.95)}
nucl_probs = {0: non_islands_one, 1: islands_one}

def viterbi(seq):
    s = 2  # number of states
    st_dict = {0: "-", 1: "+"}
    v = np.zeros(shape=(s, len(seq)))
    prev = np.zeros(shape=(s, len(seq))).astype(int)
    init_state = 0
    # v[0, 0] = np.log( transfer_probs[init_state][0]) + np.log( nucl_probs[0][seq[0]])
    # v[1, 0] = np.log( transfer_probs[init_state][1]) + np.log( nucl_probs[1][seq[0]])
    for i in range(s):
        v[i, 0] = np.log(transfer_probs[init_state][i])+np.log(nucl_probs[i][seq[0]])
    # v[:,0] /= v[:,0].sum()
    # v[1, 0] = transfer_probs[init_state][1]*nucl_probs[1][seq[0]]

    for i in range(1, len(seq)):
        for st in range(s):
            # probs = [
            #     np.log( nucl_probs[st][seq[i]]) +
            #     np.log( transfer_probs[st][j]) +
            #     v[j][i-1] for j in range(s)
            # ]
            probs = np.array([
                np.log(nucl_probs[st][seq[i]])+np.log(transfer_probs[st][j]) +
                v[j][i-1] for j in range(s)
            ])
            # probs /= probs.sum()
            max_ind = np.argmax(probs)
            prev[st, i] = max_ind
            v[st, i] = probs[max_ind]

    result = [0] * len(seq)
    final_state = np.argmax(v[:, -1])

    for i in range(len(seq)-1, -1, -1):
        result[i] = st_dict.get(final_state)
        final_state = prev[final_state, i]

    return result

def calc_metrics(true, pred):
    metrics = {"TP": 0, "TN": 0, "FP": 0, "FN": 0}
    for i, state in enumerate(true):
        equal = state == pred[i]
        if state == "-":
            if equal:
                metrics["TN"] += 1
            else:
                metrics["FP"] += 1
        else:
            if equal:
                metrics["TP"] += 1
            else:
                metrics["FN"] += 1

    metrics["Recall"] = f'{metrics["TP"]/(metrics["TP"]+ metrics["FN"]):.3f}' if (metrics["TP"]+ metrics["FN"]) != 0 else 0
    metrics["Precision"] = f'{metrics["TP"]/(metrics["TP"]+metrics["FP"]):.3f}' if (metrics["TP"]+metrics["FP"]) != 0 else 1
    metrics["Accuracy"] = f'{(metrics["TP"]+metrics["TN"])/len(true):.3f}'
    return metrics

# res = viterbi(seqs[0])
# metrics = calc_metrics(states[0], res)
# print(metrics)

# probs to get to current state from + and -
back_probs = {0: (0.995, 0.005), 1: (0.05, 0.95)}
# nucl_probs = {0: non_islands_one, 1: islands_one}


def forward_backward(seq):
    s = 2  # number of states
    st_dict = {0: "-", 1: "+"}
    v = np.zeros(shape=(s, len(seq)))
    back = np.zeros(shape=(s, len(seq)))
    # prev = np.zeros(shape=(s, len(seq))).astype(int)
    init_state = 0
    for i in range(s):
        v[i, 0] = np.log(transfer_probs[init_state][i]) + np.log(
            nucl_probs[i][seq[0]])
        # v[i, 0] = transfer_probs[init_state][i])*nucl_probs[i][seq[0]]

    for i in range(1, len(seq)):
        for st in range(s):
            # probs = np.array([
            #     np.log(nucl_probs[st][seq[i]])+np.log(transfer_probs[st][j]) +
            #     v[j][i-1] for j in range(s)
            # ])
            log_probs = np.array([
                np.log(transfer_probs[st][j]) +
                v[j][i - 1] for j in range(s)
            ])
            max_log = max(log_probs)
            sum_exp_probs = np.exp(log_probs - max_log).sum()

            # max_ind = np.argmax(probs)
            # prev[st, i] = max_ind
            v[st, i] = np.log(sum_exp_probs) + max_log + np.log(
                nucl_probs[st][seq[i]])
            print(v[:, i])

    back[0, -1] = transfer_probs[0][0]
    back[1, -1] = transfer_probs[1][1]
    for i in range(len(seq)-2, -1, -1):
        for st in range(s):
            # probs = np.array([
            #     np.log(nucl_probs[st][seq[i]])+np.log(transfer_probs[st][j]) +
            #     v[j][i-1] for j in range(s)
            # ])
            log_probs = np.array([
                np.log(back_probs[st][j]) + np.log(nucl_probs[j][seq[i+1]]) +
                back[j][i + 1] for j in range(s)
            ])
            max_log = max(log_probs)
            sum_exp_probs = np.exp(log_probs - max_log).sum()

            # max_ind = np.argmax(probs)
            # prev[st, i] = max_ind
            back[st, i] = np.log(sum_exp_probs) + max_log
        print(back[:, i])

    result = [0] * len(seq)
    final_state = np.argmax(v[:, -1])
    # final_state = random.choices([0, 1], weights = v[:, -1]/v[:, -1]), k=1)[0]
    # for i in range(len(seq)-1, -1, -1):
    # result[i] = st_dict.get(final_state)
    # final_state = prev[final_state, i]

    return result

res = forward_backward(seqs[0])

print("END")
