"""Given a DNA string transcribe it into RNA string """


if __name__ == "__main__":
    with open('rosalind_dna.txt', 'r') as f:
        str_dna = f.read()

    str_rna = ''.join(map(lambda x: 'U' if x == 'T' else x, str_dna))
    print(str_rna)
