def dp_change(money, coins):
    min_num_coins = [0]
    for m in range(1, money+1):
        min_num_coins.append(1e10)
        for coin in coins:
            if m >= coin:
                if min_num_coins[m - coin] + 1 < min_num_coins[m]:
                    min_num_coins[m] = min_num_coins[m - coin] + 1
    return min_num_coins[-1]


def manhatten_tourist(n, m, down, right):
    s = [[0] * (m+1) for _ in range(n+1)]
    s[0][0] = 0
    for i in range(1, n):
        s[i][0] = s[i-1][0] + down[i-1][0]
    for j in range(1, m):
        s[0][j] = s[0][j-1] + right[0][j-1]
    for i in range(1, n+1):
        for j in range(1, m+1):
            s[i][j] = max(s[i-1][j] + down[i-1][j], s[i][j-1] + right[i][j-1])
    return s[-1][-1]


def lsc_backtrack(v, w):
    s = [[0] * (len(w)+1) for _ in range(len(v)+1)]
    backtrack = [[None] * (len(w)+1) for _ in range(len(v)+1)]
    for i in range(len(v) + 1):
        s[i][0] = 0
    for j in range(len(w) + 1):
        s[0][j] = 0
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            s[i][j] = max(s[i-1][j], s[i][j-1])
            if v[i-1] == w[j-1]:
                s[i][j] = max(s[i][j], s[i-1][j-1] + 1)
            if s[i][j] == s[i-1][j]:
                backtrack[i][j] = "↓"
            elif s[i][j] == s[i][j-1]:
                backtrack[i][j] = "→"
            elif s[i][j] == s[i-1][j-1]+1 and v[i-1] == w[j-1]:
                backtrack[i][j] = "↘"
    return backtrack


def output_lcs(backtrack, v, i, j):
    lcs = []
    while True:
        if i == 0 or j == 0:
            break
        if backtrack[i][j] == '↓':
            i -= 1
        elif backtrack[i][j] == '→':
            j -= 1
        else:
            lcs.append(v[i-1])
            i -= 1
            j -= 1
    return ''.join(reversed(lcs))


def dag_longest_path(start, end, dag):
    # Path weights
    dag_weights = {start: 0}
    # Save predecessors
    dag_paths = {start: None}
    for k in sorted(dag):
        v = dag[k]
        n = k[1]
        print("n = ", n)
        if n not in dag_weights:
            if k[0] in dag_weights:
                dag_weights[n] = dag_weights[k[0]] + v
                dag_paths[n] = k[0]
        else:
            if k[0] in dag_weights:
                if dag_weights[k[0]] + dag[k] > dag_weights[n]:
                    dag_weights[n] = dag_weights[k[0]] + dag[k]
                    dag_paths[n] = k[0]

    # Restoring path
    path = [end]
    n = dag_paths[end]
    while n != None:
        path.append(n)
        n = dag_paths[n]
    path = '->'.join(map(str, reversed(path)))

    return dag_weights[end], path

if __name__ == "__main__":
    # (1) Dynamic Change
    money = 18431
    coins = "16,15,12,10,5,3,1"

    coins = list(map(int, coins.split(',')))
    print(dp_change(money, coins))

    # (2) Manhatten longest path
    with open('input.txt', 'r') as f:
        firstline = f.readline()
        nrows = int(firstline.split(' ')[0])
        ncols = int(firstline.split(' ')[1])
        lines = [s.strip().split(' ') for s in f.readlines()]
        mdown = [list(map(int, l)) for l in lines[:lines.index(['-'])]]
        mright = [list(map(int, l)) for l in lines[lines.index(['-'])+1:]]

        print(manhatten_tourist(nrows, ncols, mdown, mright))

    # (3) Output LCS
    v = "AACCTTGG"
    w = "ACACTGTGA"
    s = output_lcs(lsc_backtrack(v, w), v, len(v), len(w))
    print(s)

    # (4) Longest path in DAG
    with open('input.txt', 'r') as f:
        start_node = int(f.readline())
        end_node = int(f.readline())
        lines = f.readlines()
        dag = {}
        for l in lines:
            input_node = int(l[:l.index('->')])
            output_node = int(l[l.index('->')+2:l.index(':')])
            weight = int(l[l.index(':')+1:])
            dag[(input_node, output_node)] = weight
        path_len, path = dag_longest_path(start_node, end_node, dag)
        print(path_len)
        print(path)
