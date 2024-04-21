from collections import defaultdict

import numpy as np


def levenshtein_dist(str_a, str_b, match, mismatch, delete, insert):
    # assume str_a is shorter, swap them if not
    if len(str_a) > len(str_b):
        str_a, str_b = str_b, str_a

    # create a list with length of the shortest string and fill it with seq
    # from 0 to length
    dp_line = [i*delete for i in range(len(str_a)+1)]

    # outer loop over the longest string
    for i in range(1, len(str_b)+1):
        # create current value as we need to compare 3 elements,
        # 2 of them are in the same line in dp_line list.
        # Remaining is in the cur_val, for the beginning of the string
        # it equals the insert
        cur_val = i*insert
        for j in range(1, len(str_a)+1):

            # if elements of the strings are the same set param to 0, else 1
            if str_a[j - 1] == str_b[i - 1]:
                param = match
            else:
                param = mismatch

            # choose the maximum and assign it to new_val
            new_val = max(
                dp_line[j] + insert,
                cur_val + delete,
                dp_line[j-1] + param
            )

            # we can update the previous element in dp_line as we don't
            # need it anymore for next calculations
            dp_line[j-1] = cur_val
            # update cur_val
            cur_val = new_val

        # finally, after inner loop update the last element of dp_line
        dp_line[-1] = cur_val

    return dp_line[-1]


def needleman_wunsch(seq1, seq2, match=2, mismatch=-1, delete=-2, insert=-2):
    n = len(seq1)
    m = len(seq2)
    score_matrix = np.zeros((n + 1, m + 1))
    score_matrix[:, 0] = [delete * i for i in range(n+1)]
    score_matrix[0, :] = [insert * i for i in range(m + 1)]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq2[j - 1] == seq1[i - 1]:
                param = match
            else:
                param = mismatch
            score_matrix[i, j] = max(
                score_matrix[i - 1, j - 1] + param,
                score_matrix[i - 1, j] + delete,
                score_matrix[i, j - 1] + insert
            )

    res_1 = ''
    res_2 = ''

    while n > 0 or m > 0:
        diag, up, left = -float("inf"), -float("inf"), -float("inf")
        if m > 0 and n > 0:
            if seq2[m - 1] == seq1[n - 1]:
                param = match
            else:
                param = mismatch
            diag = score_matrix[n-1, m-1] + param
            left = score_matrix[n, m-1] + insert
            up = score_matrix[n-1, m] + delete
        elif m > 0:
            left = score_matrix[n, m - 1] + insert
        else:
            up = score_matrix[n - 1, m] + delete

        best = max(diag, left, up)
        if best == diag:
            res_1 = seq1[n - 1] + res_1
            res_2 = seq2[m - 1] + res_2
            n -= 1
            m -= 1
        elif best == left:
            res_1 = '-' + res_1
            res_2 = seq2[m - 1] + res_2
            m -= 1
        else:
            res_1 = seq1[n - 1] + res_1
            res_2 = '-' + res_2
            n -= 1

    return res_1, res_2


def mult_align(sequences, match=2, mismatch=-1, delete=-2, insert=-2):
    seqs = sequences.copy()
    profiles = [-1] * len(seqs)
    cons_dict = {i:[i] for i in range(len(seqs))}
    pairwise = np.full(shape=(len(seqs), len(seqs)), fill_value=-float("inf"))
    used = set()
    final_consensus = None
    for i in range(len(seqs)-1):
        for j in range(i+1, len(seqs)):
            pairwise[i][j] = levenshtein_dist(
                seqs[i], seqs[j], match=match, mismatch=mismatch,
                delete=delete, insert=insert
            )
    for stage in range(len(seqs)-1):
        i, j = np.unravel_index(np.argmax(pairwise), pairwise.shape)
        prof_i, prof_j = needleman_wunsch(
            seqs[i], seqs[j], match=match, mismatch=mismatch,
            delete=delete, insert=insert
        )
        # if profiles[i] == -1:
        #     profiles[i] = prof_i
        # if profiles[j] == -1:
        #     profiles[j] = prof_j
        for pr in cons_dict[i]:
            profiles[pr] = prof_i
        for pr in cons_dict[j]:
            profiles[pr] = prof_j
        cons_dict[i] += cons_dict[j]
        del cons_dict[j]
        consensus = ""
        for letter in range(len(prof_i)):
            cur_letters = defaultdict(int)
            for seq in cons_dict[i]:
                cur_letters[profiles[seq][letter]] += 1
            consensus += max(cur_letters.items(), key=lambda x: (x[1], x[0]))[0]
        seqs[i] = consensus
        final_consensus = consensus
        # used.add(i)
        used.add(j)
        recount = np.full(len(seqs), fill_value=-float("inf"))
        for s in range(len(seqs)):
            if s in used or s == i:
                continue
            else:
                recount[s] = levenshtein_dist(
                    consensus, seqs[s], match=match, mismatch=mismatch,
                    delete=delete, insert=insert
                )
        pairwise[:i, i] = recount[:i]
        pairwise[i, i:] = recount[i:]
        pairwise[j, :] = -float("inf")
        pairwise[:, j] = -float("inf")
    final = list()
    print(final_consensus)
    for seq in sequences:
        result, _ = needleman_wunsch(
            seq, final_consensus, match=match, mismatch=mismatch,
            delete=delete, insert=insert
        )
        final.append(result)
    return final

delete = -2
insert = -2
match = 2
mismatch = -1

seqs = ["ATTC", "GAGTC", "AGATTC", "GATTC", "GGATTC", "CGATTCA"]

res = mult_align(seqs, match=match, mismatch=mismatch, delete=delete, insert=insert)
print(*seqs)
for align in res:
    print(align)
