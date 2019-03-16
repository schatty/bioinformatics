"""
Given two DNA string compute their humming distance.
"""

def hamming_distance(str_1, str_2):
    return sum(map(lambda pair: pair[0] != pair[1], list(zip(str_1, str_2))))

if __name__ == "__main__":
    with open('rosalind_dna.txt', 'r') as f:
        dna_str_1, dna_str_2 = list(map(lambda x: x.strip(), f.readlines()))

    print(hamming_distance(dna_str_1, dna_str_2))
