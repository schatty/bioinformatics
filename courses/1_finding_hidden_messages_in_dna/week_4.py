import random
import numpy as np
from week_3 import profile_most_probable_kmer, profile_motifs, score_motifs
from week_3 import proba_profile, entropy


def select_motifs(profile, dna):
    """ Compute most probable kmers for each of dna string from profile """
    k = len(profile["A"])
    return [profile_most_probable_kmer(s, profile, k) for s in dna]


def select_random_motifs(dna, k, t):
    """ Select random k-mer motif per string in dna """
    motifs = []
    for t in range(len(dna)):
        i = random.randint(0, len(dna[0]) - k)
        motifs.append(dna[t][i:i+k])
    return motifs


def randomized_motif_search(dna, k, t, iter=1000):
    """ Perform randomized motif search algorithm on dna strings """
    motifs = select_random_motifs(dna, k, t)
    best_motifs = motifs
    for i in range(iter):
        motifs = select_random_motifs(dna, k, t)
        while True:
            prfl = profile_motifs(motifs)
            motifs = select_motifs(prfl, dna)
            if score_motifs(motifs, laplace_rule=False) < score_motifs(best_motifs, laplace_rule=False):
                best_motifs = motifs[:]
            else:
                break
    return best_motifs


def generate_rand_kmer(s, profile, k):
    """ Generate random kmer from string s using weighted dice technique """
    kmers = [s[i:i+k] for i in range(len(s) - k + 1)]
    probas = [proba_profile(s[i:i+k], profile) for i in range(len(s) - k + 1)]
    sum_probas = float(sum(probas))
    probas = [p / sum_probas for p in probas]
    probas_lined = [probas[0]]
    for p in probas[1:]:
        probas_lined.append(p + probas_lined[-1])
    rand_n = np.random.random()
    for i in range(len(probas_lined)):
        if rand_n < probas_lined[i]:
            return kmers[i]


def gibbs_sampler(dna, k, t, N):
    """ Find most promising motifs with gibbs sampler algorithm """
    global_motifs = []
    global_score = np.inf
    for global_step in range(20):
        motifs = select_random_motifs(dna, k, t)
        best_motifs = motifs
        for j in range(N):
            i = np.random.randint(0, t-1)
            profile = profile_motifs(motifs[:i] + motifs[i+1:])
            motif_i = generate_rand_kmer(dna[i], profile, k)
            motifs = motifs[:i] + [motif_i] + motifs[i+1:]
            if score_motifs(motifs, laplace_rule=False) < score_motifs(best_motifs, laplace_rule=False):
                best_motifs = motifs[:]

        if score_motifs(best_motifs, laplace_rule=False) < global_score:
            global_score = score_motifs(best_motifs, laplace_rule=False)
            global_motifs = best_motifs[:]
    return global_motifs


if __name__ == "__main__":
    # (1) Test Motifs from Profile and Dna
    prfl = {
    "A": [4./5, 0, 0, 1./5],
    "C": [0, 3./5, 1./5, 0],
    "G": [1./5, 1./5, 4./5, 0],
    "T": [0, 1./5, 0, 4./5]
    }
    dna = [
    "ttaccttaac".upper(),
    "gatgtctgtc".upper(),
    "acggcgttag".upper(),
    "ccctaacgag".upper(),
    "cgtcagaggt".upper()]
    #print(select_motifs(prfl, dna))

    # (2) Test Randomized Motif Search
    dna = [
    'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
    'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
    'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
    'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
    'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'
    ]
    k = 8
    t = 5
    #print(randomized_motif_search(dna, k, t))

    # (3) Test Gibbs Sampling
    k = 8
    t = 5
    N = 100
    dna = [
    'CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA',
    'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
    'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
    'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
    'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'
    ]
    #print('\n'.join(gibbs_sampler(dna, k, t, N)))

    # Quiz Questions

    dna = ["ATGAGGTC",
    "GCCCTAGA",
    "AAATAGAT",
    "TTGTGCTA"]
    motifs = ["GTC", "CCC", "ATA", "GCT"]
    print(select_motifs(profile_motifs(motifs), dna))
