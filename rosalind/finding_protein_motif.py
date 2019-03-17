"""
Given a list of proteins' ids from unifor find indices with N glycosilation
for each of the protein.
"""
import urllib3


def check_N_glycosylation_hardcoded(protein_substr):
    """
    Return True if given srting satisfies hardcoded requirementsself.

    protein_substr (str): protein substring

    Return (bool): True if all requirements were satisfied
    """
    fatal_mistakes = [
        protein_substr[0] != 'N',
        protein_substr[1] == 'P',
        protein_substr[2] != 'S' and protein_substr[2] != 'T',
        protein_substr[3] == 'P'
    ]
    return sum(fatal_mistakes) == 0


def get_indices(protein_str, valid_func, motif_len=4):
    """
    Return indices of protein string containing required motif determined by
    validation function.

    protein_str (str): protein string
    valid_func (function): funtion for substring validation
    motif_len (int): length of motif to validate

    Return (list): list of indices from which motifs stared in a string
    """
    return [i+1 for i in range(len(protein_str) - motif_len + 1)
                if valid_func(protein_str[i:i+motif_len])]

def get_uniprot_protein(protein_id):
    """
    Return protein string from its id from uniprot database

    protein_id (str): id of the protein

    Return (str): protein string 
    """
    http = urllib3.PoolManager()
    request = http.request('GET',
        f'http://www.uniprot.org/uniprot/{protein_id}.fasta')
    str = request.data.decode('utf-8').strip()

    # Return string after first \n with no following newlines
    return str[str.find('\n')+1:].replace('\n', '')


if __name__ == "__main__":
    with open('input.txt', 'r') as f:
        lines = list(map(lambda x: x.strip(), f.readlines()))
    for prot_id in lines:
        inds = get_indices(get_uniprot_protein(prot_id), check_N_glycosylation_hardcoded)
        # If protein has N_glycosylation print it and corresponding indices
        if len(inds):
            print(prot_id)
            print(' '.join(list(map(str, inds))))
