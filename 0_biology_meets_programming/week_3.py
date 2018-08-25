def Count(Motifs):
    """
    Input: A set of kmers Motifs
    Output: Count(Motifs)
    """
    count = {}

    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)

    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return count


def Profile(Motifs):
    """
    Input: A list of kmers Motifs
    Output: the profile matrix of Motifs, as a dictionary of lists.
    """
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}

    profile = Count(Motifs)
    for nuc in ['A', 'C', 'G', 'T']:
        for i in range(k):
            profile[nuc][i] = profile[nuc][i] / t
    return profile


def Consensus(Motifs):
    """
    Input: A set of kmers Motifs
    Output: A consensus string of Motifs
    """
    count = Count(Motifs)
    k = len(Motifs[0])
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


def Score(Motifs):
    """
    Input: A set of k-mers Motifs
    Output: The score of these k-mers
    """
    consensus_str = Consensus(Motifs)
    score = 0
    for i in range(len(consensus_str)):
        for j in range(len(Motifs)):
            score += Motifs[j][i] != consensus_str[i]
    return score


def Pr(Text, Profile):
    """
    Input: String Text and profile matrix Profile
    Output: Pr(Text, Profile)
    """
    p = 1
    for i in range(len(Text)):
        p *= Profile[Text[i]][i]
    return p

def ProfileMostProbableKmer(text, k, profile):
    profile = dict(zip(["A", "C", "G", "T"], profile))
    best_p = -1
    best_str = None
    for i in range(len(text) - k + 1):
        pb = Pr(text[i:i+k], profile)
        if pb > best_p:
            best_p = pb
            best_str = text[i:i+k]
    return best_str


def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


"""
sample_input = ["AACGTA",
"CCCGTT",
"CACCTT",
"GGATTA",
"TTCCGG"]
print(Score(sample_input))
"""


if __name__ == "__main__":
    # Profile Most Probable k-mer
    Text = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
    k = 5
    profile = [[0.2, 0.2, 0.3, 0.2, 0.3],
    [0.4, 0.3, 0.1, 0.5, 0.1],
    [0.3, 0.3, 0.5, 0.2, 0.4],
    [0.1, 0.2, 0.1, 0.1, 0.2]]

    #print(ProfileMostProbableKmer(Text, k, profile))

    profile = {
    "A": [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
    "C": [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
    "G": [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
    "T": [0.3, 0.1, 0.0, 0.4, 0.5, 0.0]
    }
    print(Pr("GAGCTA", profile))
