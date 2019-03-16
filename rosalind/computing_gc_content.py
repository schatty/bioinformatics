"""
Given DNA string in FASTA format return ID of the string with the highest
GC content and the value GC-content value itself.
"""
import operator


def compute_gc_content(dna):
    return (dna.count('G') + dna.count('C')) / len(dna)


if __name__ == "__main__":
    # Read file by lines
    with open('rosalind_gc.txt', 'r') as f:
        lines = list(map(lambda x: x.strip(), f.readlines()))

    # Calculate GC content for each full dna string
    current_key = lines[0]
    current_dna = ''
    gc_dict = {}
    for line in lines[1:]:
        if line.startswith('>Rosalind'):
            gc_dict[current_key] = compute_gc_content(current_dna)
            current_dna = ''
            current_key = line
        else:
            current_dna += line
    gc_dict[current_key] = compute_gc_content(current_dna)

    # Get key with max value
    max_gc_content_str = max(gc_dict.items(), key=operator.itemgetter(1))[0]

    # Print ID with max GC content value and corrensponding GC value
    print(max_gc_content_str[1:])
    print(gc_dict[max_gc_content_str] * 100)
