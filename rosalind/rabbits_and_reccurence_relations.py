"""
Given number of months and a litter of each pair return number of pairs
after n months
 """

def rabbits_pair_count(n_months, litter):
    """
    n_months (int): number of months to edure
    litter (int): number of pairs from one rabbit pair

    Return (int): number of pairs
    """
    if n_months <= 2:
        return 1
    n_pairs_cur = 1
    n_pairs_prev = 1
    for i in range(n_months-2):
        n_pairs_cur_tmp = n_pairs_cur
        n_pairs_cur = n_pairs_prev * litter + n_pairs_cur
        n_pairs_prev = n_pairs_cur_tmp
    return n_pairs_cur

if __name__ == "__main__":
    print(rabbits_pair_count(33, 2))
