def pattern_count(source_text, pattern):
    """
    Count number of occurences of pattern in source_text.
    """
    return sum((source_text[i:i+len(pattern)] == pattern
        for i in range(len(source_text) - len(pattern) + 1)))

def frequent_words(st, k):
    """
    Count most frequent words of length k in the given string.
    Complexity: O(n^2)
    """
    counts = [pattern_count(st, st[i:i+k]) for i in range(len(st) - k + 1)]
    max_count = max(counts)
    print("Max count: ", max_count)
    frequent_patterns = []
    for i in range(len(counts)):
        if counts[i] == max_count and st[i:i+k] not in frequent_patterns:
            frequent_patterns.append(st[i:i+k])
    return frequent_patterns

### Faster Frequent words section ###

def symbol_to_number(s):
    """ Convert given symbol to corresponding nucleotide index. """
    d = {"A": 0, "C": 1, "G": 2, "T": 3}
    return d[s]

def number_to_symbol(n):
    """ Convert given index to the corresponding nucleotide. """
    d = {0: "A", 1: "C", 2: "G", 3: "T"}
    return d[n]

def pattern_to_number(pattern):
    """ Convert given string pattern to number. """
    n = 0   # accumulator
    m = 1   # multiplier
    while len(pattern):
        n += m * symbol_to_number(pattern[-1])
        m *= 4
        pattern = pattern[:-1]
    return n

def number_to_pattern(n, l):
    """ Convert given number to the string with alphabet of AGCT. """
    pattern = ""
    while l:
        pattern = number_to_symbol(n % 4) + pattern
        n //= 4
        l -= 1
    return pattern

def compute_frequencies(st, k):
    """
    Count most frequent words of length k in the given string.
    """
    freq_arr = [0] * 4**k
    for i in range(len(st) - k + 1):
        pattern = st[i:i+k]
        freq_arr[pattern_to_number(pattern)] += 1
    return freq_arr

def faster_frequent_words(s, k):
    """ Compute frequent words by building accumulator length of 4**k  """
    frequent_patterns = []
    freq_arr = compute_frequencies(s, k)
    max_count = max(freq_arr)
    for i in range(4**k):
        if freq_arr[i] == max_count:
            pattern = number_to_pattern(i, k)
            frequent_patterns.append(pattern)
    return frequent_patterns

### ... ###

### Finding Frequent Words By Sorting Section

def frequent_words_by_sorting(s, k):
    """
    Compute most frequent words by counting k-mers, sorting it and
    finding most repeatable sequence of k-merself.
    """
    frequent_patterns = []
    inds = []
    count = []
    for i in range(0, len(s) - k + 1):
        inds.append( pattern_to_number( s[i:i+k] ) )
        count.append(1)
    sorted_index = sorted(inds)
    for i in range(1, len(s) - k + 1):
        if sorted_index[i] == sorted_index[i - 1]:
            count[i] = count[i - 1] + 1
    max_count = max(count)
    for i in range(len(s) - k + 1):
        if count[i] == max_count:
            pattern = number_to_pattern(sorted_index[i], k)
            if pattern not in frequent_patterns:
                frequent_patterns.append(pattern)
    return frequent_patterns

### ... ###

def reverse_complement(s):
    """ Get reverse complement of DNA area """
    d = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(map(lambda x: d[x], reversed(s)))

def pattern_matching(s, pattern):
    """ Return indices of pattern beggining in the string s. """
    inds = []
    for i in range(len(s) - len(pattern) + 1):
        if s[i:i+len(pattern)] == pattern:
            inds.append(i)
    return inds

def find_clumps(g, k, l, t):
    """
    Retrun founded clumps in string g.
    Clump is the t occurences of the same k-mer in the substring of length l.
    """
    clumps = []
    kmers = dict()

    def collect_kmer(kmers, t, acc):
        for kmer in kmers.keys():
            if kmers[kmer] >= t and kmer not in clumps:
                clumps.append(kmer)

    # Form dictionary for the first l symbols
    for i in range(l - k + 1):
        if g[i:i+k] not in kmers.keys():
            kmers[g[i:i+k]] = 1
        else:
            kmers[g[i:i+k]] += 1

    collect_kmer(kmers, t, clumps)

    for i in range(1, len(g) - l + 1, 1):
        if g[i+l-k:i+l] not in kmers.keys():
            kmers[g[i+l-k:i+l]] = 1
        else:
            kmers[g[i+l-k:i+l]] += 1
        collect_kmer(kmers, t, clumps)
        kmers[g[i:i+k]] -= 1
    return clumps


if __name__ == "__main__":
    # (1) Test pattern_count
    source_text = "ACTGTACGATGATGTGTGTCAAAG"
    pattern = "TGT"
    #print(pattern_count(source_text, pattern))

    # (2) Test frequent words
    source_text = "TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT"
    k = 3
    print(frequent_words(source_text, k))

    # (3) Reverse compelement Test
    #print(reverse_complement("GCTAGCT"))

    # (4) Pattern counts
    #print(pattern_matching("AAACATAGGATCAAC", "AA"))

    # (4) Clump problem
    g = "..."
    k = 11
    l = 519
    t = 17
    #print(find_clumps(g, k, l, t))

    # CHARGE STATION
    # (1) Test pattern to number
    #print(pattern_to_number("AGT"))

    # (2) Test number to pattern
    #print(number_to_pattern(5217, 8))

    # (3) Test compute_frequencies
    #inds = #compute_frequencies("ATAGTAGTTAACAGACACTGATATATTCTTCTGTTCGGGTCGCAACGTAAACCAAGTTCGCGGGGTGCTGGCCCCAGTACAATTGAGTCACACCGGTAGTCTGGAGGCAACCGAAGTTCCCCCGCATTGCAGACCGGCGTGAACTTGTCGAGGCGCTCCTCGAGTCGGACACGCCAGGAGCTACCTAGGTGGTGTCAGCGGCCTATGCGTGTGCGCTCTTGGCTTAGCATTGCTCTGATTGAGTGACGAATACACAAGGTGGGCGTTCCAATTGCCCTAGGCCCAATCCCGTTGGGATTTCGCGGTGTGCCTGTTCGGAGCCAACAGATACAGTTAGTTACTTCTAGCTAACCCATGAGGGTGTGAACCACAAAAAATTGAACCCGCGTCGGTCTAGTGGACGCTAATGGTTTGCTGCACCCACCCTACTAAGAGGAGGGCATCCGATGCAGTAGTCCAGTGAGACCAGGTCATGGAGCTATGCAATAACTGCATACCGACATAACCAGTAACTCATACTGAGCTATGCAGGCGTGCACCCGCCCACTCTCGGCGATCCCACCAATGCAATTTCGCGCCAAACGGTAGAGCCGAAGGCGCTCTCTCTCTTGCCCCTCGTGACTGTTGGGGCAATTCCATAACAAACCACGGTTTGCTTGGATCCTATGCGG", 5)
    #print(' '.join(map(str, inds)))

    # (4) Test faster faster_frequent_words
    source_text = "TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT"
    k = 3
    print(faster_frequent_words(source_text, k))

    # (5) Test frequent words by sorting
    source_text = "TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT"
    k = 3
    print(frequent_words_by_sorting(source_text, k))
