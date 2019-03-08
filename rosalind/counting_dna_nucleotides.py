"""Given a DNA string count number of A, C, G and T symbols """


if __name__ == "__main__":
    with open('rosalind_dna.txt', 'r') as f:
        str = f.read()
    print(str.count('A'), str.count('C'), str.count('G'), str.count('T'))
