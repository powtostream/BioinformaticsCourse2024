import numpy as np
from collections import defaultdict

def levenshtein_dist(str_a, str_b, match, mismatch, delete, insert):
    # assume str_a is shorter, swap them if not
    # if len(str_a) > len(str_b):
    #     str_a, str_b = str_b, str_a

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

    return np.array(dp_line)

calls = defaultdict(list)

def call_rec(str_a, str_b, match, mismatch, delete, insert, lvl=0):
    calls[lvl].append((str_a, str_b))
    print(str_a, str_b)
    if len(str_b) <= 1:
        return
    div_ind = len(str_a)//2
    a_left = str_a[:div_ind]
    a_right = str_a[div_ind:]
    left = levenshtein_dist(
        str_b, a_left,  match, mismatch, delete, insert)
    right = levenshtein_dist(
        str_b[::-1], a_right[::-1], match, mismatch, delete, insert)
    best_b_div = np.argmax((left+right[::-1]))
    b_left = str_b[:best_b_div]
    b_right = str_b[best_b_div:]
    call_rec(a_left, b_left, match, mismatch, delete, insert, lvl+1)
    call_rec(a_right, b_right, match, mismatch, delete, insert, lvl+1)

a = "AGTACGCA"
b = "TATGC"
delete = -2
insert = -2
match = 2
mismatch = -1
call_rec(a, b, match, mismatch, delete, insert)
