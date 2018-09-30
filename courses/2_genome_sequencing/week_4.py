from collections import Counter


def cyclic_spectrum(peptide, amino_mass):
    """ Return cyclic spectrum of given peptide """
    prefix_mass = [0]
    for i in range(len(peptide)):
        for k in amino_mass.keys():
            if peptide[i] == k:
                prefix_mass.append(prefix_mass[-1] + amino_mass[k])
                break
    peptide_mass = prefix_mass[-1]
    cyclic_spectrum = [0]
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide)+1):
            print(i, j, len(prefix_mass))
            cyclic_spectrum.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < len(peptide):
                cyclic_spectrum.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))
    return sorted(cyclic_spectrum)


def linear_spectrum(peptide, amino_mass):
    """ Return linear spectrum of given peptide """
    prefix_mass = [0]
    for i in range(len(peptide)):
        for k in amino_mass.keys():
            if peptide[i] == k:
                prefix_mass.append(prefix_mass[-1] + amino_mass[k])
                break
    linear_spectrum = [0]
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide)+1):
            linear_spectrum.append(prefix_mass[j] - prefix_mass[i])
    return sorted(linear_spectrum)


def score_cyclopeptide(pep, spec, amino_mass):
    """ Count matching subpeptides between given peptide and spectrum """
    count_orig = Counter(spec)
    count_pep = Counter(cyclic_spectrum(pep, amino_mass))
    n_matching = 0
    for c in count_orig.keys():
        if count_orig[c] <= count_pep[c]:
            n_matching += count_orig[c]
        else:
            n_matching += count_pep[c]
    return n_matching


def score_linear_peptide(pep, spec, amino_mass):
    """ Count matching linear subpeptides between given peptide and spectrum """
    count_orig = Counter(spec)
    count_pep = Counter(linear_spectrum(pep, amino_mass))
    n_matching = 0
    for c in count_orig.keys():
        if count_orig[c] <= count_pep[c]:
            n_matching += count_orig[c]
        else:
            n_matching += count_pep[c]
    return n_matching


def trim(leaderboard, spec, N, amino_mass):
    """ Return all peptides which scores are in top-N """
    lin_scores = []
    for pep in leaderboard:
        lin_scores.append( score_linear_peptide(pep, spec, amino_mass) )

    # Sort leaderboard according to decreasing order of scores
    lin_scores_sorted = sorted(lin_scores, reverse=True)
    leaderboard_sorted = []
    for sc in lin_scores_sorted:
        i_score = lin_scores.index(sc)
        lin_scores[i_score] = -1 # For excluding repeats
        leaderboard_sorted.append(leaderboard[i_score])

    for j in range(N, len(leaderboard)):
        if lin_scores_sorted[j] < lin_scores_sorted[N-1]:
            return leaderboard_sorted[:j]
    return leaderboard


def leaderboard_cyclopeptide_sequencing(spectrum, N):
    def expand(peptides):
        pep_symbols = "GASPVTCILNDKQEMHFRYW"
        expanded_peptides = []
        for pep in peptides:
            expanded_peptides += [pep+s for s in pep_symbols]
        return set(expanded_peptides)

    leaderboard = [""]
    leader_peptide = ""
    leader_score = -1
    parent_mass = spectrum[-1]
    while len(leaderboard):
        leaderboard = list(expand(leaderboard))
        for pep in leaderboard[:]:
            if calc_mass(pep) == parent_mass:
                score = score_cyclopeptide(pep, spectrum)
                if score > leader_score:
                    leader_peptide = pep
                    leader_score = score
            elif calc_mass(pep) > parent_mass:
                leaderboard.remove(pep)
        leaderboard = trim(leaderboard, spectrum, N)
    return leader_peptide


def convolution_cyclopeptide_sequencing(spectrum, M, N):
    """ Find most promising peptide with leaderboard and convolution. """

    spectrum = sorted(spectrum)
    leaderboard = [[]]
    leader_peptide = []
    leader_score = -1
    parent_mass = spectrum[-1]

    # For extended alphabet from convolution
    exts = spectral_convolution(spectrum)
    exts_filtered = [a for a in exts if 57 <= a <= 200]
    exts_counter = Counter(exts_filtered)
    spectrum_conv = []
    for v in sorted(set(exts_counter.values()), reverse=True)[:M]:
        for k in exts_counter:
            if exts_counter[k] == v and v > 1:
                spectrum_conv.append(k)

    spectrum_conv = sorted(spectrum_conv)
    amino_mass = {k:k for k in spectrum_conv}

    def expand(peptides, spec):
        expanded_peptides = []
        for pep in peptides:
            expanded_peptides += [pep+[s] for s in spec]
        return expanded_peptides

    def calc_mass(pep):
        return sum(pep)

    while len(leaderboard):
        leaderboard = list(expand(leaderboard, list(amino_mass.keys())))
        for pep in leaderboard[:]:
            pep_mass = calc_mass(pep)
            if pep_mass == parent_mass:
                score = score_cyclopeptide(pep, spectrum, amino_mass)
                if score > leader_score:
                    leader_peptide = pep
                    leader_score = score
            elif pep_mass > parent_mass:
                leaderboard.remove(pep)
        leaderboard = trim(leaderboard, spectrum, N, amino_mass)
    return leader_peptide


def spectral_convolution(spec, debug=False):
    """ Compute spectral convolution """
    counter_orig = Counter(spec)
    row = sorted(spec)
    if row[0] != 0:
        row = [0] + row
    counter = Counter()
    for i in range(1, len(row)):
        for j in range(i):
            if row[i] - row[j]:
                counter[row[i] - row[j]] += 1
    subs = [[k]*v for k, v in sorted(counter.items())]
    if debug:
        print(counter.items())
        print("Max Multiplicity: ", max(counter.values()))
    l = []
    for ll in subs:
        l += ll
    return l


if __name__ == "__main__":
    # (1) Test Cyclopeptide Scoring
    cycpep = "NQEL"
    spectrum = list(map(int, "0 99 113 114 128 227 257 299 355 356 370 371 484".split()))
    #print(score_cyclopeptide(cycpep, spectrum))

    # (2) Test Linear Cyclopeptide Scoring
    #pep = "NQEL"
    #spec = list(map(int, "0 99 113 114 128 227 257 299 355 356 370 371 484".split()))
    #print(score_linear_peptide(pep, spec))

    # (3) Test Trim
    leaderboard = "LAST ALST TLLT TQAS".split()
    spec = list(map(int, "0 71 87 101 113 158 184 188 259 271 372".split()))
    N = 2
    #print(" ".join(trim(leaderboard, spec, N)))

    # (4) Test Leaderboard Cyclopeptide Sequencing
    N = 10
    spec = list(map(int, "0 71 113 129 147 200 218 260 313 331 347 389 460".split()))
    #acid = leaderboard_cyclopeptide_sequencing(spec, N)
    #print('-'.join(map(str, [amino_mass[c] for c in acid])))

    # (5) Test spectral convolution
    #spec = list(map(int, "0 137 186 323".split()))
    #print(spectral_convolution(spec))

    # (6) Test convolution cyclopeptide sequencing
    M = 20
    N = 60
    #spec = list(map(int, "57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493".split()))
    #acid = convolution_cyclopeptide_sequencing(spec, M, N)
    #print('-'.join(map(str, acid)))

    # Quiz Questions

    # Problem 3
    amino_mass = {"G": 57, "A": 71, "S": 87, 'P': 97, "V": 99,
    "T": 101, "C": 103, 'I': 113, "L": 113, "N": 114, "D": 115,
    "K": 128, "Q": 128, 'E': 129, "M": 131, "H": 137, "F": 147,
    "R": 156, "Y": 163, "W": 186}
    pep = "MAMA"
    spec = list(map(int, "0 71 98 99 131 202 202 202 202 202 299 333 333 333 ".split()))
    print("Cyclic Score = ", score_cyclopeptide(pep, spec, amino_mass))

    # Problem 4
    pep = "PEEP"
    spec = list(map(int, "0 97 97 129 194 196 226 226 244 258 323 323 452".split()))
    print("Linear Score = ", score_linear_peptide(pep, spec, amino_mass))

    # Problem 5
    spec = list(map(int, "0 57 118 179 236 240 301".split()))
    spectral_convolution(spec, debug=True)
