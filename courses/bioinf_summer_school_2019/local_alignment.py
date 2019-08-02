def get_match_score(x, y):
    if x == y:
        return 1
    return -1


def output_lcs_local(backtrack, s, v, w, i, j):
    w_trace = []
    v_trace = []
    while True:
        if i == 0 or j == 0:
            break
        if s[i][j] == 0:
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
    print("V: ", v_res)
    print("W: ", w_res)
    return v_res, w_res


def local_alignment(v, w):
    sigma = 1
    s = [[0] * (len(w)+1) for _ in range(len(v)+1)]
    backtrack = [['~'] * (len(w)+1) for _ in range(len(v)+1)]
    for j in range(1, len(w)+1):
        s[0][j]= s[0][j-1] - sigma
    for i in range(1, len(v)+1):
        s[i][0] = s[i-1][0] - sigma

    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            match_score = get_match_score(v[i-1], w[j-1])
            s[i][j] = max(s[i-1][j] - sigma, s[i][j-1] - sigma, s[i-1][j-1]+match_score, 0)
            if s[i][j] == s[i-1][j] - sigma:
                backtrack[i][j] = "↓"
            elif s[i][j] == s[i][j-1] - sigma:
                backtrack[i][j] = "→"
            elif s[i][j] == 0:
                backtrack[i][j] = 'taxi'
            else:
                backtrack[i][j] = "↘"

    max_score = max([max(a) for a in s])
    end_point = (0, 0)
    for i in range(len(s)):
        for j in range(len(s[0])):
            if s[i][j] == max_score:
                end_point = (i, j)
                break
    v_res, w_res = output_lcs_local(backtrack, s, v, w, end_point[0], end_point[1])
    return v_res, w_res, max_score


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
    print("V: ", v_res)
    print("W: ", w_res)
    return v_res, w_res


def fitting_alignment(v, w):

    def output_lcs_fitting(backtrack, v, w, i, j):
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
                v_trace.append(v[i-1])
                w_trace.append(w[j-1])
                i -= 1
                j -= 1

        if j != 0:
            w_trace.append(w[j-1])
            v_trace.append('-')

        v_res = ''.join(reversed(v_trace))
        w_res = ''.join(reversed(w_trace))
        print("V: ", v_res)
        print("W: ", w_res)
        return v_res, w_res

    sigma = 1
    s = [[0] * (len(w)+1) for _ in range(len(v)+1)]
    backtrack = [['~'] * (len(w)+1) for _ in range(len(v)+1)]
    for j in range(1, len(w)+1):
        s[0][j]= s[0][j-1] - sigma
    for i in range(1, len(v)+1):
        s[i][0] = s[i-1][0]# - sigma

    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            if v[i-1] == w[j-1]:
                match_score = 1
            else:
                match_score = -1
            if j == 1:
                s[i][j] = max(s[i-1][j] - sigma, s[i][j-1] - sigma, s[i-1][j-1]+match_score, 0)
            else:
                s[i][j] = max(s[i-1][j] - sigma, s[i][j-1] - sigma, s[i-1][j-1]+match_score)

            if s[i][j] == s[i-1][j] - sigma:
                backtrack[i][j] = "↓"
            elif s[i][j] == s[i][j-1] - sigma:
                backtrack[i][j] = "→"
            elif s[i][j] == s[i][j-1] and j <= 1:
                backtrack[i][j] = "→"
            elif s[i][j] == 0:
                backtrack[i][j] = '→'
            else:
                backtrack[i][j] = "↘"

    max_sc = -1e10
    max_ind = 0
    for i in range(len(s)):
        if s[i][-1] > max_sc:
            max_sc = s[i][-1]
            max_ind = i
        #print(s[i][-1])
    v_res, w_res = output_lcs_fitting(backtrack, v, w, max_ind, len(w))
    return v_res, w_res, max_sc


if __name__ == "__main__":
    s1 = "ACCTACAGGATTATCTTGCGGCGG"
    s2 = "GGTTGTCTT"

    v, w, score = fitting_alignment(s1, s2)
    print(score)
    print(v)
    print(w)