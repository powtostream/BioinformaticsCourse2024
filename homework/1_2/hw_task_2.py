import numpy as np
from Bio import SeqIO
from Bio.Align import substitution_matrices, PairwiseAligner

f8_path = "data/f8.fasta"
gattaca_path = "data/gattaca.fasta"
f8_seq = []
f8 = SeqIO.parse(f8_path, "fasta")
for record in f8.records:
    f8_seq.append(str(record.seq))

gattaca = SeqIO.parse(gattaca_path, "fasta")
gattaca_seq = []
for record in gattaca.records:
    gattaca_seq.append(str(record.seq))
substitution_matrix = substitution_matrices.load('BLOSUM62')

def needleman_wunsch(seq1, seq2, substitution_matrix, gap):
    n = len(seq1)
    m = len(seq2)
    score_matrix = np.zeros((n + 1, m + 1))
    score_matrix[:, 0] = [gap * i for i in range(n+1)]
    score_matrix[0, :] = [gap * i for i in range(m + 1)]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            subst = substitution_matrix[seq1[i - 1]][seq2[j - 1]]
            score_matrix[i, j] = max(
                score_matrix[i - 1, j - 1] + subst,
                score_matrix[i - 1, j] + gap,
                score_matrix[i, j - 1] + gap
            )
    res_1 = ''
    res_2 = ''
    while n > 0 or m > 0:
        if n > 0 and m > 0 and (
                score_matrix[n, m] == score_matrix[n - 1, m - 1] +
                substitution_matrix[seq1[n - 1]][seq2[m - 1]]
        ):
            res_1 = seq1[n - 1] + res_1
            res_2 = seq2[m - 1] + res_2
            n -= 1
            m -= 1
        elif n > 0 and score_matrix[n, m] == score_matrix[n - 1, m] + gap:
            res_1 = seq1[n - 1] + res_1
            res_2 = '-' + res_2
            n -= 1
        else:
            res_1 = '-' + res_1
            res_2 = seq2[m - 1] + res_2
            m -= 1
    print(res_1)
    for i in range(len(res_1)):
        if res_1[i] == res_2[i]:
            print("|", end="")
        elif res_1[i] == "-" or res_2[i] == "-":
            print(" ", end="")
        else:
            print(".", end="")
    print("")
    print(res_2)
    return res_1, res_2, score_matrix[-1, -1]


def affine_gap(seq1, seq2, substitution_matrix, alpha, beta):
    n = len(seq1)
    m = len(seq2)
    score_matrix = np.zeros((n + 1, m + 1))
    score_matrix[1:, 0] = [alpha + beta * (i-1) for i in range(1, n + 1)]
    score_matrix[0, 1:] = [alpha + beta * (i-1) for i in range(1, m + 1)]
    # score_matrix[0,0] = 0
    # seq_1_gap = np.full(shape=(n + 1, m + 1), fill_value=-np.infty)
    # seq_2_gap = np.full(shape=(n + 1, m + 1), fill_value=-np.infty)
    seq_1_gap = np.zeros((n + 1, m + 1))
    seq_2_gap = np.zeros((n + 1, m + 1))
    seq_1_gap[1:, 0] = [alpha + beta*(i-1) for i in range(1, n + 1)]
    seq_2_gap[0, 1:] = [alpha + beta*(i-1) for i in range(1, m + 1)]
    # seq_2_gap[0,0] = alpha
    # seq_1_gap[0, 0] = alpha

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            seq_1_gap[i, j] = max(
                score_matrix[i,j-1] + alpha,
                seq_1_gap[i,j-1] + beta
            )
            seq_2_gap[i, j] = max(
                score_matrix[i-1, j] + alpha,
                seq_2_gap[i-1, j] + beta
            )
            subst = substitution_matrix[seq1[i - 1]][seq2[j - 1]]
            score_matrix[i, j] = max(
                score_matrix[i - 1, j - 1] + subst,
                seq_1_gap[i, j],
                seq_2_gap[i, j]
                # seq_1_gap[i-1, j-1] + subst,
                # seq_2_gap[i-1, j-1] + subst
            )
            # seq_1_gap[i, j] = max(
            #     score_matrix[i - 1, j - 1] + alpha + beta,
            #     seq_1_gap[i-1, j] + beta
            # )
            # seq_2_gap[i, j] = max(
            #     score_matrix[i - 1, j - 1] + alpha + beta,
            #     seq_2_gap[i, j-1] + beta
            # )
    res_1 = ''
    res_2 = ''
    while n > 0 or m > 0:

        if n > 0 and score_matrix[n, m] == seq_1_gap[n, m]:
            res_1 = seq1[n - 1] + res_1
            res_2 = '-' + res_2
            n -= 1
        elif n > 0 and m > 0 and (
                score_matrix[n, m] == score_matrix[n - 1, m - 1] +
                substitution_matrix[seq1[n - 1]][seq2[m - 1]]
        ):
            res_1 = seq1[n - 1] + res_1
            res_2 = seq2[m - 1] + res_2
            n -= 1
            m -= 1
        else:
        # elif m > 0 and score_matrix[n, m] == seq_2_gap[n, m]:
            res_1 = '-' + res_1
            res_2 = seq2[m - 1] + res_2
            m -= 1
        # else:
        #     res_1 = seq1[n - 1] + res_1
        #     res_2 = seq2[m - 1] + res_2
        #     n -= 1
        #     m -= 1
    # while n > 0 or m > 0:
    #     if n > 0 and m > 0 and (
    #             score_matrix[n, m] == score_matrix[n - 1, m - 1] +
    #             substitution_matrix[seq1[n - 1]][seq2[m - 1]]
    #     ):
    #         res_1 = seq1[n - 1] + res_1
    #         res_2 = seq2[m - 1] + res_2
    #         n -= 1
    #         m -= 1
    #     elif n > 0 and score_matrix[n, m] == seq_1_gap[n, m]:
    #         res_1 = seq1[n - 1] + res_1
    #         res_2 = '-' + res_2
    #         n -= 1
    #     else:
    #         res_1 = '-' + res_1
    #         res_2 = seq2[m - 1] + res_2
    #         m -= 1
    print(seq1)
    print(seq2)

    print(res_1)
    for i in range(len(res_1)):
        if res_1[i] == res_2[i]:
            print("|", end="")
        elif res_1[i] == "-" or res_2[i] == "-":
            print(" ", end="")
        else:
            print(".", end="")
    print("")
    print(res_2)
    print(score_matrix)
    print("----------")
    print(seq_1_gap)
    print("----------")
    print(seq_2_gap)
    return res_1, res_2, score_matrix[-1, -1]


gap = -3
seq1 = "ACGTAASDFDDD"
seq2 = "ACTTAADFQQ"
# seq1 = "MMMFRERRRY"
# seq2 = "MNFRY"
alignment, _, score = affine_gap(f8_seq[0][20:40], f8_seq[1][8:20], substitution_matrix, -2, -1)


#
# alignment, _, score = affine_gap(seq1, seq2, substitution_matrix, -10, -2)
print("Optimal Alignment:")
print(alignment, _)
print("Alignment Score:", score)
