"""
Given a collection of dna strings find the largest common substring
"""

def get_largest_commong_substr(strings):

    def substr_exist(strings, substr):
        for s in strings:
            if substr not in s:
                return False
        return True

    common_subs = ['A', 'C', 'G', 'T']
    for k in range(2, len(min(strings, key=len))):
        tmp_subs = []
        for nuc in 'ACGT':
            for common_sub in common_subs:
                if substr_exist(strings, common_sub + nuc):
                    tmp_subs.append(common_sub + nuc)
        if len(tmp_subs):
            common_subs = tmp_subs
        else:
            break
    return common_subs


if __name__ == "__main__":
    with open('input.txt', 'r') as f:
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

    common_subsrings = get_largest_commong_substr(dna_strings)
    print(common_subsrings[0])
