def get_match_score(x, y):
    if x == y:
        return 1
    return -1


def global_alignment_score(v, w):
    sigma = 1
    s = [[0] * (len(w)+1) for _ in range(len(v)+1)]
    for j in range(1, len(w)+1):
        s[0][j]= s[0][j-1] - sigma
    for i in range(1, len(v)+1):
        s[i][0] = s[i-1][0] - sigma

    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            match_score = get_match_score(v[i-1], w[j-1])
            s[i][j] = max(s[i-1][j] - sigma, s[i][j-1] - sigma, s[i-1][j-1]+match_score)

    return s[-1][-1]


def global_alignment(v, w):
    sigma = 1
    s = [[0] * (len(w)+1) for _ in range(len(v)+1)]
    backtrack = [['~'] * (len(w)+1) for _ in range(len(v)+1)]
    for j in range(1, len(w)+1):
        s[0][j] = s[0][j-1] - sigma
    for i in range(1, len(v)+1):
        s[i][0] = s[i-1][0] - sigma

    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            match_score = get_match_score(v[i-1], w[j-1])
            s[i][j] = max(s[i-1][j] - sigma, s[i][j-1] - sigma, s[i-1][j-1]+match_score)

            if s[i][j] == s[i-1][j] - sigma:
                backtrack[i][j] = "↓"
            elif s[i][j] == s[i][j-1] - sigma:
                backtrack[i][j] = "→"
            else:
                backtrack[i][j] = "↘"

    v_res, w_res = output_lcs(backtrack, v, w, len(v), len(w))
    return v_res, w_res, s[-1][-1]


def output_lcs(backtrack, v, w, i, j):
    w_trace = []
    v_trace = []
    while True:
        if i == 0 or j == 0:
            break
        if backtrack[i][j] == '↓':
            v_trace.append(v[i-1])
            w_trace.append('-')
            i -= 1
        elif backtrack[i][j] == '→':
            w_trace.append(w[j-1])
            v_trace.append('-')
            j -= 1
        else:
            if v[i-1] != w[j-1]:
                v_trace.append(v[i - 1].lower())
                w_trace.append(w[j-1].lower())
            else:
                v_trace.append(v[i-1])
                w_trace.append(w[j-1])
            i -= 1
            j -= 1

    # Start of the strings
    if i == 0 and j != 0:
        w_trace.append(w[j-1])
        v_trace.append('-')
    elif j == 0 and i != 0:
        v_trace.append(v[i-1])
        w_trace.append('-')

    v_res = ''.join(reversed(v_trace))
    w_res = ''.join(reversed(w_trace))

    return v_res, w_res


if __name__ == "__main__":
    # Problem 1
    s1 = 'GGATTATCT'
    s2 = 'GGTTGTCT'
    #print(global_alignment_score(s1, s2))

    # Problem 2
    s1 = 'GGATTATCT'
    s2 = 'GGTTGTCT'
    v, w, score = global_alignment(s1, s2)
    print(score)
    print(v)
    print(w)