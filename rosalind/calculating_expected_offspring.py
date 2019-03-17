"""
Given six nonnegative integers corresponding to the folloging genotypes:
AA-AA, AA-Aa, AA-aa, Aa-Aa, Aa-aa, aa-aa return expected number of offspring
displaying the dominant fenotype in the next generation, under the assumption
that every couple has exactly two offspring.
"""

def calc_dominant_offspring_expectation(n_couples, n_offspring=2):
    """
    Return expected number of offspring displaying the dominant fenotype.
    """
    expected_AA_AA = n_offspring * n_couples[0]
    expected_AA_Aa = n_offspring * n_couples[1]
    expected_AA_aa = n_offspring * n_couples[2]
    expected_Aa_Aa = n_offspring * n_couples[3] * 0.75
    expected_Aa_aa = n_offspring * n_couples[4] * 0.5

    return expected_AA_AA + expected_AA_Aa + expected_AA_aa + expected_Aa_Aa + expected_Aa_aa


if __name__ == "__main__":
    print(calc_dominant_offspring_expectation([17442, 17921, 16086, 16493, 18649, 18889]))
