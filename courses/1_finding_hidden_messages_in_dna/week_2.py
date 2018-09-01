from week_1 import pattern_to_number, number_to_pattern, reverse_complement

def skew(s):
    """ Returns list of differences G-C in the first i-th element """
    diff_list = [0]
    for i in range(len(s)):
        if s[i] == 'C':
            diff_list.append( diff_list[-1] - 1  )
        elif s[i] == 'G':
            diff_list.append( diff_list[-1] + 1 )
        else:
            diff_list.append( diff_list[-1] )
    return diff_list


def minimum_skew(s):
    """ Return indices at which skew function achieves its minimum """
    diff_values = skew(s)
    min_skew = min(diff_values)
    inds = []
    for i in range(len(s)):
        if diff_values[i] == min_skew:
            inds.append(i)
    return inds


def maximum_skew(s):
    """ Return indices at which skew function achieves its maximum """
    diff_values = skew(s)
    max_skew = max(diff_values)
    inds = []
    for i in range(len(s)):
        if diff_values[i] == max_skew:
            inds.append(i)
    return inds


def hamming_distance(s1, s2):
    """ Return number of mismatching symbols in the strings s1 and s2 """
    return sum(map(lambda i: s1[i] != s2[i], range(len(s1))))


def approx_pattern_matching(s, pattern, d):
    """
    Find indices of pattern begging in the string s.
    Substring consider to be a matching if its hamming distance with pattern
    is less of equal with d
    """
    inds = []
    for i in range(len(s) - len(pattern) + 1):
        if hamming_distance(pattern, s[i:i+len(pattern)]) <= d:
            inds.append(i)
    return inds


def approx_pattern_count(s, pattern, d):
    """
    Count number of pattern matching substrings in given string.
    """
    return len(approx_pattern_matching(s, pattern, d))


### Approximate Frequent Worlds problem section ###

def frequent_words_with_mismatches(s, k, d):
    """ Compute frequent k-mers with mismatches in string d """
    frequent_patterns = []
    freq_arr = compute_frequences_with_mismatches(s, k, d)
    max_count = max(freq_arr)
    for i in range(4**k):
        if freq_arr[i] == max_count:
            pattern = number_to_pattern(i, k)
            frequent_patterns.append(pattern)
    return frequent_patterns


def compute_frequences_with_mismatches(s, k, d):
    freq_arr = [0] * 4**k
    for i in range(len(s) - k + 1):
        pattern = s[i:i+k]
        neighbourhood = neighbours(pattern, d)
        for approx_p in neighbourhood:
            j = pattern_to_number(approx_p)
            freq_arr[j] += 1
    return freq_arr


def frequent_words_with_mismatches_and_reverse_complements(s, k, d):
    """ Function with the long name. """
    frequent_patterns = []
    freq_arr = compute_frequences_with_mismatches_and_reverse_complements(s, k, d)
    max_count = max(freq_arr)
    for i in range(4**k):
        if freq_arr[i] == max_count:
            pattern = number_to_pattern(i, k)
            frequent_patterns.append(pattern)
    return frequent_patterns


def compute_frequences_with_mismatches_and_reverse_complements(s, k, d):
    freq_arr = [0] * 4**k
    for i in range(len(s) - k + 1):
        for pattern in [s[i:i+k], reverse_complement(s[i:i+k])]:
            neighbourhood = neighbours(pattern, d)
            for approx_p in neighbourhood:
                freq_arr[pattern_to_number(approx_p)] += 1
    return freq_arr


def neighbours(pattern, d):
    """
    Return set of neighbours for given patternself.
    Substring is considered as neighbour if its hamming distance id less then d
    """
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return {'A', 'C', 'G', 'T'}
    neighbourhood = set()
    suffix_neighbours = neighbours(pattern[1:], d)
    for s in suffix_neighbours:
        if hamming_distance(pattern[1:], s) < d:
            for nuc in ['A', 'C', 'G', 'T']:
                neighbourhood.add(nuc + s)
        else:
            neighbourhood.add(pattern[0] + s)
    return neighbourhood


### ... ###

if __name__ == "__main__":
    # (1) Test Skew function
    s = "GAGCCACCGCGATA"
    #print(' '.join(map(str, skew(s))))

    # (2) Test minimum_skew function
    s = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"
    #print(' '.join(map(str, minimum_skew(s))))

    # (3) Test hamming distance
    s1 = "GGGCCGTTGGT"
    s2 = "GGACCGTTGAC"
    #print(hamming_distance(s1, s2))

    # (4) Test approximate pattern matching problem
    s = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"
    pattern = "ATTCTGGA"
    d = 3
    #print(approx_pattern_matching(s, pattern, d))

    # (5) Test approximate pattern counting problem
    s = "TTTAGAGCCTTCAGAGG"
    pattern = "GAGG"
    d = 2
    #print(approx_pattern_count(s, pattern, d))

    # (6) Test neighbour
    pattern = "GTTAAGGTCCT"
    n = 2
    with open('neigbours.txt', 'w') as f:
        f.write('\n'.join(neighbours(pattern, n)))

    # (7) Test frequent words with mismatches
    s = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
    k = 4
    d = 1
    #print(' '.join(frequent_words_with_mismatches(s, k, d)))

    # (8) Test frequent words with mismatches and reverse complements
    s = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
    k = 4
    d = 1
    #print(' '.join(frequent_words_with_mismatches_and_reverse_complements(s, k, d)))

    # Quiz 2
    s1 = "CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA"
    s2 = "CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG"
    print(hamming_distance(s1, s2))

    print(maximum_skew("GCATACACTTCCCAGTAGGTACTG"))

    print(approx_pattern_count("CGTGACAGTGTATGGGCATCTTT", "TGT", 1))

    print(len(neighbours("CCAGTCAATG", 1)))
