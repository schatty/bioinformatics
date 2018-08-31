def pattern_count(source_text, pattern):
    """
    Count number of occurences of pattern in source_text.
    """
    return sum((source_text[i:i+len(pattern)] == pattern
        for i in range(len(source_text) - len(pattern) + 1)))

def frequent_words(st, k):
    """
    Fount most frequent words of length k in the given string.
    """
    counts = [pattern_count(st, st[i:i+k]) for i in range(len(st) - k + 1)]
    max_count = max(counts)
    frequent_patterns = []
    for i in range(len(counts)):
        if counts[i] == max_count and st[i:i+k] not in frequent_patterns:
            frequent_patterns.append(st[i:i+k])
    return frequent_patterns

def reverse_complement(s):
    d = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(map(lambda x: d[x], reversed(s)))

def pattern_matching(s, pattern):
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
    #print(frequent_words(source_text, k))

    # (3) Reverse compelement Test
    #print(reverse_complement("GCTAGCT"))

    # (4) Pattern counts
    print(pattern_matching("AAACATAGGATCAAC", "AA"))

    # (4) Clump problem
    g = "..."
    k = 11
    l = 519
    t = 17
    #print(find_clumps(g, k, l, t))
