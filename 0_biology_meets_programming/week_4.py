
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
    profile = CountWithPseudocounts(Motifs)
    num_nuc = [float(profile["A"][i] + profile["C"][i] + profile["G"][i] + profile["T"][i]) for i in range(k)]
    for sym in "ACGT":
        for i in range(k):
            profile[sym][i] /= num_nuc[i]
    return profile


if __name__ == "__main__":

    # Test CountWithPseudocounts
    Motifs = ["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"]
    print(CountWithPseudocounts(Motifs))

    # Test ProfileWithPseudocounts
    print(ProfileWithPseudocounts(Motifs))
