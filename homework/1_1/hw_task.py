def hamming_dist(str_a, str_b):
    # iterate over strings and count if the elements are different
    return sum(str_a[i] != str_b[i] for i in range(len(str_a)))


def levenshtein_dist(str_a, str_b):
    # assume str_a is shorter, swap them if not
    if len(str_a) > len(str_b):
        str_a, str_b = str_b, str_a

    # create a list with length of the shortest string and fill it with seq
    # from 0 to length
    dp_line = [i for i in range(len(str_a)+1)]

    # outer loop over the longest string
    for i in range(1, len(str_b)+1):
        # create current value as we need to compare 3 elements,
        # 2 of them are in the same line in dp_line list.
        # Remaining is in the cur_val, for the beginning of the string
        # it equals the index
        cur_val = i
        for j in range(1, len(str_a)+1):

            # if elements of the strings are the same set param to 0, else 1
            if str_a[j - 1] == str_b[i - 1]:
                param = 0
            else:
                param = 1

            # choose the minimum and assign it to new_val
            new_val = min(dp_line[j] + 1, cur_val, dp_line[j-1] + param)

            # we can update the previous element in dp_line as we don't
            # need it anymore for next calculations
            dp_line[j-1] = cur_val
            # update cur_val
            cur_val = new_val

        # finally, after inner loop update the last element of dp_line
        dp_line[-1] = cur_val

    return dp_line[-1]



def closest_match(key_str, search_str):
    result = len(key_str)
    pos = None
    # iterate over string and count hamming dist for each position
    for i in range(len(search_str)-len(key_str)+1):
        h_dist = hamming_dist(key_str, search_str[i:i+len(key_str)])
        # update result and position if h_dist is smaller
        if h_dist < result:
            result = h_dist
            pos = i
    if pos is not None:
        return pos, search_str[pos:pos+len(key_str)], result
    else:
        return None, None, None
