"""
Given a collection of DNA string construct graph of overlapping strings and
print all overlapping pairs.
"""

def build_dna_adj_list(dna_dict, k=3):
    """
    Construct dna adjacency list from given dna strings and suffix=prefix length

    dna_dict (dict): dictionary with keys=dna_id, values=dna string
    k (int): suffix length (=prefix length)

    Return (dict): adjacency list
    """
    dna_adj_list = {}
    for dna_id_1 in dna_dict:
        for dna_id_2 in dna_dict:
            if dna_id_1 == dna_id_2:
                continue
            if dna_dict[dna_id_1][-k:] == dna_dict[dna_id_2][:k]:
                if dna_id_1 not in dna_adj_list:
                    dna_adj_list[dna_id_1] = [dna_id_2]
                else:
                    dna_adj_list[dna_id_1].append(dna_id_2)
    return dna_adj_list


if __name__ == "__main__":
    with open('input.txt', 'r') as f:
        lines = list(map(lambda x: x.strip(), f.readlines()))
    dna_dict = {}
    current_key = lines[0]
    current_dna = ''
    for line in lines[1:]:
        if line.startswith('>Rosalind'):
            dna_dict[current_key] = current_dna
            current_dna = ''
            current_key = line
        else:
            current_dna += line
    dna_dict[current_key] = current_dna

    # Print each pair of adjacency list
    adj_list = build_dna_adj_list(dna_dict)
    for k in adj_list:
        for v in adj_list[k]:
            print(k[1:], v[1:])
