import random

from week_3 import Consensus, Score, Count, Pr, ProfileMostProbablePattern

def CountWithPseudocounts(Motifs):
    """
    Input:  A set of kmers Motifs
    Output: CountWithPseudocounts(Motifs)
    """
    t = len(Motifs)
    k = len(Motifs[0])

    count = {
        "A": [1] * k,
        "C": [1] * k,
        "G": [1] * k,
        "T": [1] * k
    }

    for sym in "ACGT":
        for i in range(t):
            for j in range(k):
                if Motifs[i][j] == sym:
                    count[sym][j] += 1
    return count


def ProfileWithPseudocounts(Motifs):
    # Input:  A set of kmers Motifs
    # Output: ProfileWithPseudocounts(Motifs)
    k = len(Motifs[0])
    t = len(Motifs)
    profile = CountWithPseudocounts(Motifs)
    num_nuc = [float(profile["A"][i] + profile["C"][i] + profile["G"][i] + profile["T"][i]) for i in range(k)]
    for sym in "ACGT":
        for i in range(k):
            profile[sym][i] /= float(t+4)
    return profile


def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


def Motifs(Profile, Dna):
    k = len(Profile['A'])
    motifs = []
    for i in range(len(Dna)):
        motifs.append(ProfileMostProbablePattern(Dna[i], k, Profile))
    return motifs


def RandomMotifs(Dna, k, t):
    motifs = []
    for i in range(len(Dna)):
        start_ind = random.randint(0, len(Dna[i])-k)
        motifs.append(Dna[i][start_ind:start_ind+k])
    return motifs


def RandomizedMotifSearch(Dna, k, t, N):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    for i in range(N):
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
    return BestMotifs


def Normalize(Probabilities):
    sm = sum(Probabilities.values())
    for k in Probabilities.keys():
        Probabilities[k] /= sm
    return Probabilities


def WeightedDie(Probabilities):
    kmer = '' # output variable
    probas_1d = {}
    point = 0
    for k in Probabilities.keys():
        probas_1d[k] = (point, point+Probabilities[k])
        point += Probabilities[k]
    p_pick = random.uniform(0, 1)

    for k in probas_1d.keys():
        if probas_1d[k][0] <= p_pick < probas_1d[k][1]:
            return k


def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}

    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)

    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)


def GibbsSampler(Dna, k, t, N):
    BestMotifs = RandomMotifs(Dna, k, t)
    for j in range(N):
        i = random.randint(0, t-1)
        Profile = ProfileWithPseudocounts(BestMotifs[:i] + BestMotifs[i+1:])
        motif_i = ProfileGeneratedString(Dna[i], Profile, k)
        Motifs = BestMotifs[:i] + [motif_i] + BestMotifs[i+1:]
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


if __name__ == "__main__":
    '''
    # Test CountWithPseudocounts
    Motifs = ["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"]
    print(CountWithPseudocounts(Motifs))

    # Test ProfileWithPseudocounts
    print(ProfileWithPseudocounts(Motifs))

    # Test GreedyMotifSearch
    k = 3
    t = 5
    Dna = [
        "GGCGTTCAGGCA",
        "AAGAATCAGTCA",
        "CAAGGAGTTCGC",
        "CACGTCAATCAC",
        "CAATAATATTCG"
    ]
    print(GreedyMotifSearchWithPseudocounts(Dna, k, t))

    Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC",
       "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG",
       "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC",
       "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC",
       "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG",
       "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA",
       "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA",
       "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG",
       "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG",
       "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]

    # set t equal to the number of strings in Dna and k equal to 15
    t = len(Dna)
    k = 15

    # Call GreedyMotifSearchWithPseudocounts(Dna, k, t) and store the output in a variable called Motifs
    Motifs = GreedyMotifSearchWithPseudocounts(Dna, k, t)

    # Print the Motifs variable
    print(Motifs)

    # Print Score(Motifs)
    print(Score(Motifs))
    '''

    # Question 3

    Dna = ["TGACGTTC",
    "TAAGAGTT",
    "GGACGAAA",
    "CTGTTCGC"]
    k = 3
    t = 4
    N = 1

    M = ["TGA", "GTT", "GAA", "TGT"]
    BestMotifs = M

    Profile = ProfileWithPseudocounts(M)
    M = Motifs(Profile, Dna)
    if Score(M) < Score(BestMotifs):
        BestMotifs = M
    print("Best motifs: ", BestMotifs)

    print(Normalize([0.15, 0.6, 0.225, 0.225, 0.3]))
