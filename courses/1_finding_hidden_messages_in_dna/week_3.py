import numpy as np
from week_2 import neighbours, approx_pattern_count, hamming_distance


def motif_enumeration(dna, k, d):
    """
    Find motifs length k that occures in each of dna strings with
    at most d mismatches.
    """
    dna_glued = ''.join(dna)
    kmers = set(map(lambda i: dna_glued[i:i+k], range(len(dna_glued)-k+1)))
    patterns = set()

    for kmer in kmers:
        neigbs = neighbours(kmer, d)
        for nb in neigbs:
            appear_in_all_strings = True
            for s in dna:
                if approx_pattern_count(s, nb, d) == 0:
                    appear_in_all_strings = False
                    break
            if appear_in_all_strings:
                patterns.add(nb)
    return patterns


def count_motifs(motifs, laplace_rule=True):
    """ Compute count matrix with ACGT rows and len(motifs[0]) columns. """
    count = {"A": [], "C": [], "G": [], "T": []}
    lpl = [0, 1][laplace_rule]
    for j in range(len(motifs[0])):
        for nuc in "ACGT":
            count[nuc].append( sum((s[j] == nuc for s in motifs)) + lpl )
    return count


def score_motifs(motifs):
    """ Return score from list of motifs strings. """
    counts = count_motifs(motifs)
    score = 0
    for j in range(len(counts["A"])):
        max_val = max([counts[n][j] for n in "ACGT"])
        score += sum((counts[n][j] for n in "ACGT" if counts[n][j] != max_val))
    return score


def profile_motifs(motifs, laplace_rule=True):
    """ Compute profile matrix from list of motifs strings. """
    profile = count_motifs(motifs, laplace_rule=laplace_rule)
    for j in range(len(profile["A"])):
        sum_nucl = float( sum([profile[n][j] for n in "ACGT"]) )
        for n in "ACGT":
            profile[n][j] /= sum_nucl
    return profile


def entropy(motifs):
    """ Compute entropy for the list of motif strings """
    pr = profile_motifs(motifs)
    entr = 0
    for j in range(len(pr["A"])):
        e_sum = 0
        for n in "ACGT":
            if pr[n][j] == 0:
                e_sum += 0
            else:
                e_sum += pr[n][j] * np.log2(pr[n][j])
        entr += -e_sum
    return entr


def median_string(dna, k):
    """ Return one median kmer from dna strings """
    def d(pattern, dna):
        k = len(pattern)
        res = 0
        for i in range(len(dna)):
            res += min((hamming_distance(dna[i][j:j+k], pattern)
            for j in range(len(dna[i])-k)))
        return res

    distance = np.inf
    dna_glued = ''.join(dna)
    kmers = set(map(lambda i: dna_glued[i:i+k], range(len(dna_glued)-k+1)))
    for kmer in kmers:
        if distance > d(kmer, dna):
            distance = d(kmer, dna)
            median = kmer
    return median


def median_strings(dna, k):
    """ Return multiple kmers from list of dna and kmer """
    def d(pattern, dna):
        k = len(pattern)
        res = 0
        for i in range(len(dna)):
            res += min((hamming_distance(dna[i][j:j+k], pattern)
            for j in range(len(dna[i])-k)))
        return res

    distance = np.inf
    dna_glued = ''.join(dna)
    kmers = set(map(lambda i: dna_glued[i:i+k], range(len(dna_glued)-k+1)))
    medians = []
    for kmer in kmers:
        if distance == d(kmer, dna):
            medians.append(kmer)
        if distance > d(kmer, dna):
            distance = d(kmer, dna)
            medians = []
            medians.append(kmer)
    return medians


def proba_profile(s, profile):
    """ Compute probability of the given string within profile """
    p = 1
    for i, c in enumerate(s):
        p *= profile[c][i]
    return p


def profile_most_probable_kmer(s, profile, k):
    """ Find most probable k-mer in string s according to profile """
    max_pr = -1
    kmer = ""
    for i in range(len(s) - k + 1):
        pr = proba_profile(s[i:i+k], profile)
        if pr > max_pr:
            max_pr = pr
            kmer = s[i:i+k]
    return kmer


def greedy_motif_search(dna, k, t):
    best_motifs = [s[:k] for s in dna]
    for i in range(len(dna[0]) - k + 1):
        motifs = [dna[0][i:i+k]] + [s[:k] for s in dna[1:]]
        for j in range(1, t):
            prfl = profile_motifs(motifs[:j])
            motifs[j] = profile_most_probable_kmer(dna[j], prfl, k)
        if score_motifs(motifs) < score_motifs(best_motifs):
            best_motifs = motifs
    return best_motifs

if __name__ == "__main__":
    # (1) Test Motif Enumeration
    dna = ["ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"]
    k = 3
    d = 1
    #print(motif_enumeration(dna, k, d))

    # (2) Test entropy
    motifs = [
    "TCGGGGGTTTTT",
    "CCGGTGACTTAC",
    "ACGGGGATTTTC",
    "TTGGGGACTTTT",
    "AAGGGGACTTCC",
    "TTGGGGACTTCC",
    "TCGGGGATTCAT",
    "TCGGGGATTCCT",
    "TAGGGGAACTAC",
    "TCGGGTATAACC"
    ]
    #print(entropy(motifs))

    # (3) Test Median String
    dna = [
    'GCGGTGCAAGGACAGTGGCCAGTAAACCAGACTGCTCCTCGA',
    'TACCAGCCAGGCTGAGCTTCTTGACGACACGGAGCTAAACAA',
    'GTTTTAGACTTAATCAGTTACCAGCCGCCTGTTACATTCATA',
    'ATACCCGCTTGATACCAGGTCGACCAAGCATGGGACGTTCGG',
    'TACTTGTAGCGCATTACAATGGAGGATTTGAACCAGGTCCGG',
    'TCGCAAGACCAGACCGCAGGCTGCCTTGGTAAGCCCGTGCCA',
    'CTTGACTGCGCTCACCAGATGATTAGGGCATGGCGATCCTAT',
    'CGCACGGGCTGCGACCAGACGGGTAGCTTGGCGTGGTCACAA',
    'GAACATGGTAAAGAACCCGACCAGCCGTCACCTCTTGCGAGG',
    'ACTTACTACCAGGTTGAGGCCCGTTACTTTACCTGACAGTAT']
    #print(median_string(dna, 6))

    # (4) Test Probability of the Profile
    prof = {"A": [.2,.2,0,0,0,0,.9,.1,.1,.1,.3, 0],
    "C": [.1, .6, .0, .0, .0, .0, .0, .4, .1, .2, .4, .6],
    "G": [0, 0, 1, 1, .9, .9, .1, 0, 0, 0, 0, 0],
    "T": [0.7, 0.2, 0.0, 0.0, .1, .1, 0, .5, .8, .7, .3, .4]}
    #print(proba_profile("TCGTGGATTTCC", prof))

    # (5) Test Profile Most Probable K-mer Problem
    s = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
    profile = {
    "A": [0.2, 0.2, 0.3, 0.2, 0.3],
    "C": [0.4, 0.3, 0.1, 0.5, 0.1],
    "G": [0.3, 0.3, 0.5, 0.2, 0.4],
    "T": [0.1, 0.2, 0.1, 0.1, 0.2]
    }
    k = 5
    #print(profile_most_probable_kmer(s, profile, k))

    # (6) Test Greedy Motif Search
    k = 3
    t = 5
    dna = [
    'GGCGTTCAGGCA',
    'AAGAATCAGTCA',
    'CAAGGAGTTCGC',
    'CACGTCAATCAC',
    'CAATAATATTCG',
    ]
    #print(greedy_motif_search(dna, k, t))

    #### Quiz Problems ###

    def entr(x):
        r = 0
        for p in x:
            r += [p * np.log2(p), 0][p == 0]
        return -r
    print("Entropy A: ", entr([0.5, 0.0, 0.0, 0.5]))
    print("Entropy B: ", entr([0.25, 0.25, 0.25, 0.25]))
    print("Entropy C: ", entr([0, 0, 0, 1]))
    print("Entropy D: ", entr([0.25, 0, 0.5, 0.25]))

    dna = [
    'CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC',
    'GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC',
    'GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG'
    ]
    print(median_strings(dna, 7))

    prfl = {
    "A": [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
    "C": [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
    "G": [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
    "T": [0.3, 0.1, 0.0, 0.4, 0.5, 0.0]
    }
    print(proba_profile("AAGTTC", prfl))
