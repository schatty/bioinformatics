"""
Given a collection of DNA strings find consensus sting and corresponding
profile matrix.
"""
import operator

def get_column_counts(dna_strings, pos):
    """
    dna_strings (list): list of dna strings
    pos (int): position in the string (column number)
    """
    d = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for dna in dna_strings:
        d[dna[pos]] += 1
    return d


if __name__ == "__main__":
    with open('rosalind_dna.txt', 'r') as f:
        lines = list(map(lambda x: x.strip(), f.readlines()))
    dna_strings = []
    current_dna = ''
    for line in lines[1:]:
        if line.startswith('>Rosalind'):
            dna_strings.append(current_dna)
            current_dna = ''
        else:
            current_dna += line
    dna_strings.append(current_dna)

    consensus_str = ''
    profile_columns = []
    for j in range(len(dna_strings[0])):
        column_nucl = get_column_counts(dna_strings, j)
        profile_columns.append(column_nucl)
        max_count_nucleotide = max(column_nucl.items(), key=operator.itemgetter(1))[0]
        consensus_str += max_count_nucleotide

    # Print consensus sting and profile matrix
    print(consensus_str)
    for i, nuc in enumerate(['A', 'C', 'G', 'T']):
        nuc_str = ' '.join(list(map(str, [col[nuc] for col in profile_columns])))
        print(f'{nuc}: {nuc_str}')
