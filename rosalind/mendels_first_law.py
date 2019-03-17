"""
Given three integers representing populations containing homozygous dominant,
heterozygous and homozygous recessive calculate probability that two randomly
selected mating organisms will produce an individual posessing a dominant
allelle.
"""

def calc_prob_dominant(k, m, n):
    """
    k (int): homozygous dominant
    m (int): heterozygous
    n (int): homozygous recessive

    Return (float): probability of producing individual containing dominant factor.
    """
    total_num = k + m + n
    p_aA_aA = m/total_num * (m-1)/(total_num-1) * 0.25
    p_aa_Aa = n/total_num * m/(total_num-1) * 0.5 * 2
    p_aa_aa = n/total_num * (n-1)/(total_num-1)
    return 1 - p_aA_aA - p_aa_Aa - p_aa_aa

if __name__ == "__main__":
    print(calc_prob_dominant(15, 25, 29))
