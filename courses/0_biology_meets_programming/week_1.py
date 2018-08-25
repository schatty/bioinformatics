def PatternCount(Text, Pattern):
    """
    Function counts number of occurances of Pattern inside Text.
    Input:
        Text: source text to search
        Pattern: pattern to search in Text
    Output: count of occurences (int)
    """
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count


def PatternMatching(Pattern, Genome):
    """
    Return positions of the matching pattern in given DNA
    Input:
        Pattern: pattern to search for.
        Genome: source text to search in.
    """
    positions = [] # output variable
    for i in range(len(Genome) - len(Pattern) + 1):
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(i)
    return positions


def FrequencyMap(Text, k):
    """
    Construct frequency map of patterns in text
    Input:
        Text: source text
        k: pattern length
    Output: dict, key = pattern, value = number of occurences
    """
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        if Pattern not in freq:
            freq[Pattern] = 1
        else:
            freq[Pattern] += 1
    return freq


def FrequentWords(Text, k):
    """
    Input:
        Text: text to search frequent words
        k: length of word
    Ouput: list of most frequent words in given Text
    """
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words


def Reverse(Pattern):
    """
    Return reverse string.
    """
    rev = ''
    for c in Pattern:
        rev = c + rev
    return rev


def Complement(Pattern):
    """
    Construct DNA complement of given string
    Input:
        Pattern: DNA pattern
    Ouput: Complement of the given DNA pattern
    """
    comp = ''
    for c in Pattern:
        if c.lower() == 'a':
            comp += 'T'
        elif c.lower() == 'g':
            comp += 'C'
        elif c.lower() == 'c':
            comp += 'G'
        elif c.lower() == 't':
            comp += 'A'
    return comp


def ReverseComplement(Pattern):
    return Reverse(Complement(Pattern))
