def compl(s):
    d = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return d[s]


if __name__ == "__main__":
    with open('rosalind_dna.txt', 'r') as f:
        str_dna = f.read().strip()
    str_dna_compl = ''.join(reversed(list(map(compl, str_dna))))
    print(str_dna_compl)
