"""
Calculate the total number of pairs of rabits that will remain after n-th month
if all rabits live for m month.
"""

def calc_mortal_rabbits(n, m):
    """
    n (int): number of months to endure
    m (int): number of months after which rabbits die

    Return (int): total number of rabbits pairs
    """
    pairs = [1, 1]
    cur_month = 2
    while cur_month < n:
        if cur_month < m:
            pairs.append(pairs[-2] + pairs[-1])
        elif cur_month == m:
            pairs.append(pairs[-2] + pairs[-1] - 1)
        else:
            pairs.append(pairs[-2] + pairs[-1] - pairs[-m-1])
        cur_month += 1
    return pairs[-1]


if __name__ == "__main__":
    print(calc_mortal_rabbits(94, 18))
