def rna2dna(rna):
    return rna.replace('U', 'T')

def reverse_rna_complement(s):
    """ Get reverse complement of RNA area """
    d = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(map(lambda x: d[x], reversed(s)))


def read_acid_map(path):
    """ Return map with keys=codon, value=acid """
    acid_map = {}
    with open(path, 'r') as f:
        lines = f.readlines()
        for l in lines:
            if len(l.split()) > 1:
                acid = l.split()[1]
            else:
                acid = ""
            acid_map[l.split()[0]] = acid
    return acid_map


def rna_2_acid(rna, acid_map):
    """ Translate an RNA string into an amino acid string """
    return ''.join([acid_map[rna[i:i+3]] for i in range(0, len(rna)//3*3, 3)])


def find_dna_substrings_from_acid(dna, acid, acid_map):
    """ Return list of dna substrings coding given acid """
    rna = dna.replace('T', 'U')
    rna_substr_len = len(acid) * 3
    rna_sub_list = []
    for i in range(len(rna)):
        # Check original substring
        rna_sub_orig = rna[i:i+rna_substr_len]
        acid_str = ''.join([rna_2_acid(rna_sub_orig[j:j+3], acid_map) for j in range(0, len(rna_sub_orig)-len(rna_sub_orig)%3, 3)])
        if acid_str == acid:
            rna_sub_list.append(rna2dna(rna_sub_orig))
        # Check reverse complement
        rna_sub_rev = reverse_rna_complement(rna_sub_orig)
        acid_str = ''.join([rna_2_acid(rna_sub_rev[j:j+3], acid_map) for j in range(0, len(rna_sub_rev)-len(rna_sub_rev)%3, 3)])
        if acid_str == acid:
            rna_sub_list.append(rna2dna(rna_sub_orig))
    return rna_sub_list


def num_subpeptine(n):
    """ Calculate number of subpetides for peptide of length n """
    return n * (n-1)


def cyclic_spectrum(peptide):
    """ Return cyclic spectrum of given peptide """
    amino_mass = {"G": 57, "A": 71, "S": 87, 'P': 97, "V": 99,
    "T": 101, "C": 103, 'I': 113, "L": 113, "N": 114, "D": 115,
    "K": 128, "Q": 128, 'E': 129, "M": 131, "H": 137, "F": 147,
    "R": 156, "Y": 163, "W": 186}
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
            cyclic_spectrum.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < len(peptide):
                cyclic_spectrum.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))
    return sorted(cyclic_spectrum)

def linear_spectrum(peptide):
    """ Return linear spectrum of given peptide """
    amino_mass = {"G": 57, "A": 71, "S": 87, 'P': 97, "V": 99,
    "T": 101, "C": 103, 'I': 113, "L": 113, "N": 114, "D": 115,
    "K": 128, "Q": 128, 'E': 129, "M": 131, "H": 137, "F": 147,
    "R": 156, "Y": 163, "W": 186}
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

def calc_mass(peptide):
    amino_mass = {"G": 57, "A": 71, "S": 87, 'P': 97, "V": 99,
    "T": 101, "C": 103, 'I': 113, "L": 113, "N": 114, "D": 115,
    "K": 128, "Q": 128, 'E': 129, "M": 131, "H": 137, "F": 147,
    "R": 156, "Y": 163, "W": 186}
    return sum([amino_mass[a] for a in peptide])

def get_masses(acid):
    amino_mass = {"G": 57, "A": 71, "S": 87, 'P': 97, "V": 99,
    "T": 101, "C": 103, 'I': 113, "L": 113, "N": 114, "D": 115,
    "K": 128, "Q": 128, 'E': 129, "M": 131, "H": 137, "F": 147,
    "R": 156, "Y": 163, "W": 186}
    return '-'.join([str(amino_mass[a]) for a in acid])


def calc_mass(acid):
    amino_mass = {"G": 57, "A": 71, "S": 87, 'P': 97, "V": 99,
    "T": 101, "C": 103, 'I': 113, "L": 113, "N": 114, "D": 115,
    "K": 128, "Q": 128, 'E': 129, "M": 131, "H": 137, "F": 147,
    "R": 156, "Y": 163, "W": 186}
    return sum([amino_mass[a] for a in acid])


def cyclopeptide_sequencing(spectrum):

    def expand(peptides):
        pep_symbols = "GASPVTCILNDKQEMHFRYW"
        expanded_peptides = []
        for pep in peptides:
            expanded_peptides += [pep+s for s in pep_symbols]
        return set(expanded_peptides)

    def is_consistent(pep, spectrum):
        pep_spec = linear_spectrum(pep)
        for s in spectrum:
            if s in pep_spec:
                pep_spec.remove(s)
        return len(pep_spec) == 0

    seqs = []
    peptides = set([""])
    parent_mass = spectrum[-1]
    while len(peptides):
        peptides = expand(peptides)
        for pep in list(peptides):
            if calc_mass(pep) == parent_mass:
                if cyclic_spectrum(pep) == spectrum:
                    if get_masses(pep) not in seqs:
                        seqs.append(get_masses(pep))
                    peptides.remove(pep)
            elif not is_consistent(pep, spectrum):
                peptides.remove(pep)
    print(' '.join(seqs))


def num_linear_subpeptides(n):
    return sum([i for i in range(1, n+1)]) + 1

if __name__ == "__main__":
    # (1) Test protein translation
    acid_map = read_acid_map("RNA_codon_table_1.txt")
    rna_string = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
    acid_string = rna_2_acid(rna_string, acid_map)
    #print(acid_string)

    # (2) Test find_dna_substrings_from_acid
    dna_string = "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
    acid = "MA"
    #print('\n'.join(find_dna_substrings_from_acid(dna_string, acid, acid_map)))

    # (3) Test num_peptide
    #print( num_subpeptine(42067) )

    # (4) Test linear_spectrum
    #print(' '.join(map(str, cyclic_spectrum('YKPVKYYYAKDRNT') )))

    # (5) Test num_linear_subpeptides
    #print(num_linear_subpeptides(11379))

    # (6) Test cyclopeptide sequencing
    #print(cyclopeptide_sequencing([0, 113, 128, 186, 241, 299, 314, 427]))

    #with open('input.txt', 'r') as f:
    #    f_str = f.read()
    #    spectrum = list(map(int, f_str.split()))
    #print(cyclopeptide_sequencing(spectrum))

    # Quiz Questions

    # Problem 2
    print(rna_2_acid("CCAAGAACAGAUAUCAAU", acid_map) == "PRTEIN")
    print(rna_2_acid("CCACGUACUGAAAUUAAC", acid_map) == "PRTEIN")
    print(rna_2_acid("CCAAGUACAGAGAUUAAC", acid_map) == "PRTEIN")
    print(rna_2_acid("CCGAGGACCGAAAUCAAC", acid_map) == "PRTEIN")

    # Problem 3
    '''
    S -> 4
    Y -> 2
    N -> 2
    G -> 4
    E -> 2
    2^3*4*4=128
    '''

    # Problem 4
    # 'UGG' -> 186

    # Problem 5
    print("TMLA", cyclic_spectrum("TMLA"))
    print("IAMT", cyclic_spectrum("IAMT"))
    print("MLAT", cyclic_spectrum("MLAT"))
    print("TMIA", cyclic_spectrum("TMIA"))
    print("TAIM", cyclic_spectrum("TAIM"))
    print("MAIT", cyclic_spectrum("MAIT"))

    # Problem 6
    def is_consistent(pep, spectrum):
        pep_spec = linear_spectrum(pep)
        for s in spectrum:
            if s in pep_spec:
                pep_spec.remove(s)
        return len(pep_spec) == 0
    spectrum = list(map(int, "0 71 99 101 103 128 129 199 200 204 227 230 231 298 303 328 330 332 333".split()))
    print("QCV", is_consistent("QCV", spectrum))
    print("CTQ", is_consistent("CTQ", spectrum))
    print("ETC", is_consistent("ETC", spectrum))
    print("TCE", is_consistent("TCE", spectrum))
    print("CTV", is_consistent("CTV", spectrum))
    print("AQV", is_consistent("AQV", spectrum))
