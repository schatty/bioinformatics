def SymbolArray(Genome, symbol):
    # Input:  Strings Genome and symbol
    # Output: SymbolArray(Genome, symbol)
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array


# Reproduce the PatternCount function here.
def PatternCount(Pattern, Text):
    n = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            n += 1
    return n

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look ath the fist half of Genome to compute first array value
    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i - 1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i] - 1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i] + 1
    return array


def SkewArray(Genome):
    Skew = [0] * (len(Genome)+1)
    for i in range(len(Genome)):
        if Genome[i] == 'A' or Genome[i] == 'T':
            Skew[i+1] = Skew[i]
        elif Genome[i] == 'G':
            Skew[i+1] = Skew[i] + 1
        else:
            Skew[i+1] = Skew[i] - 1
    return Skew


def MinimumSkew(Genome):
    skew_array = SkewArray(Genome)
    min_val = min(skew_array)
    min_inds = []
    for i, skew in enumerate(skew_array):
        if skew == min_val:
            min_inds.append(i)
    return min_inds


def HammingDistance(p, q):
    return sum([1 for (i, j) in zip(p, q) if i != j])


def ApproximatePatternMatching(Text, Pattern, d):
    positions = [] # initializing list of positions
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions


def ApproximatePatternCount(Pattern, Text, d):
    return len(ApproximatePatternMatching(Text, Pattern, d))
