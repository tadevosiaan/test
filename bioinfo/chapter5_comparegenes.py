from tqdm import tqdm
import numpy as np
import random
from bio_constants import Down, Right, Diag, w_possible, w_ordinary, w_mu1
from bio_utils import print_matrix, HammingDistance, similar_parts


# CHAPTER 5. How Do We Compare Genes?
# Dynamic programming


def LongestCommonSubsequence(u: str, v: str) -> str:
    # Input: two strings
    # Output: a longest common subsequence of there strings
    pass


def ManhattanTourist(down: list, right: list, diag=None):
    # Input: a weighted n x m rectangular grid with n+1 rows and m+1 columns
    # Here we get a grid with weight arrays of size (n+1)x(m+1)
    # Output: a longest path from source (0,0) to sink (n,m) in the grid AND path
    # TODO: implement backtracking
    n = len(down) - 1
    m = len(down[0]) - 1
    s = [[0] * (m + 1) for i in range(n + 1)]
    s[0][0] = 0
    path = []
    for i in range(1, n + 1):
        s[i][0] = s[i - 1][0] + down[i][0]
    for j in range(1, m + 1):
        s[0][j] = s[0][j - 1] + right[0][j]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if diag is not None:
                s[i][j] = max(s[i - 1][j] + down[i][j], s[i][j - 1] + right[i][j], s[i - 1][j - 1] + diag[i][j])
                if s[i][j] == s[i - 1][j] + down[i][j]:
                    path.append((i - 1, j))
                elif s[i][j] == s[i][j - 1] + right[i][j]:
                    path.append((i, j - 1))
                elif s[i][j] == s[i - 1][j - 1] + diag[i][j]:
                    path.append((i - 1, j - 1))
            else:
                s[i][j] = max(s[i - 1][j] + down[i][j], s[i][j - 1] + right[i][j])
                if s[i][j] == s[i - 1][j] + down[i][j]:
                    path.append((i - 1, j))
                elif s[i][j] == s[i][j - 1] + right[i][j]:
                    path.append((i, j - 1))
    print_matrix(s)
    print(path)
    return s[n][m], path


def LongestPathInDirectedGraph(graph):
    # Input: an edge-weighted directed graph with source and sink nodes
    # Output: a longest path from source to sink in the directed graph
    pass


# 5A Code challenge
def DPChange(m: int, coins: list) -> int:
    # Input: int m, Coins = list of positive integers
    # Output: the minimum number of coins with denominators Coins that changes money
    # changes <=> money \in span(Coins) over R_+, so m=a_1 c_1 + ... + a_d c_d, a_i>=0
    # Complexity: O(m * |coins|) time, O(m) space -> can decrease
    # TODO: implement backtracking
    d = [0 for i in range(m + 1)]
    for i in tqdm(range(1, m + 1)):

        a = []
        for coin in coins:
            if i >= coin:
                a.append(d[i - coin] + 1)
        d[i] = min(a)
    print(d, '\n', d[m - 1], sep='')
    return d[m - 1]


def topological_sort(G):
    # Input: graph G (assuming DirectedAcylic)
    # Output: topological order of nodes
    # TODO: write it
    pass


def compute_score(b, G):
    # Input: node b, DAG G
    # Output: score s[b] for LongestPath
    pass


def LongestPathInDAG(G, source, sink):
    # Input: an edge-weighted directed acyclic graph with source and sink nodes
    # Output: a longest path from source to sink in the directed acyclic graph
    # TODO: rewrite it properly to work on arbitrary DA graphs
    score = dict()
    for b in G.nodes:  # assume G has nodes
        score[b] = float('-inf')
    score[source] = 0
    t = topological_sort(G)
    for b in t:
        score[b] = compute_score(b, G)
    return score[sink]


def global_alignment_matrix(str1, str2, weights, gap_weight):
    n, m = len(str1), len(str2)
    F = [[0] * (m + 1) for i in range(n + 1)]
    for i in range(n):
        F[i][0] = gap_weight * i
    for j in range(m):
        F[0][j] = gap_weight * j
    F[0][0] = 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = F[i - 1][j - 1] + weights[str1[i - 1]][str2[j - 1]]
            delete = F[i - 1][j] + gap_weight
            insert = F[i][j - 1] + gap_weight
            F[i][j] = max(match, delete, insert)
    return F


def local_alignment_matrix(str1, str2, weights, gap_weight):
    n, m = len(str1), len(str2)
    # F = [[0] * (m + 1) for i in range(n + 1)]
    F = np.zeros((n + 1, m + 1))
    for i in range(n):
        F[i, 0] = max(0, gap_weight * i)
    for j in range(m):
        F[0, j] = max(0, gap_weight * j)
    F[0, 0] = 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if i == n and j == m:
                F[i, j] = np.max(F)
            else:
                match = F[i - 1][j - 1] + weights[str1[i - 1]][str2[j - 1]]
                delete = F[i - 1][j] + gap_weight
                insert = F[i][j - 1] + gap_weight
                F[i][j] = max(0, match, delete, insert)
    return F


def LocalAlignment(a, b, w, d, F):
    # Input: strings a,b, matrix of weights w, weight of gap d, their matrix of local alignment F, computed by nw_matrix
    # Output: their local alignment
    # TODO: rewrite
    alignment_a = ''
    alignment_b = ''
    i, j = len(a), len(b)
    while i > 0 or j > 0:
        score = F[i][j]
        # print('1st while -> i: {}, j:{}'.format(i, j))
        score_source = 0
        score_diag = F[i - 1][j - 1]
        score_up = F[i][j - 1]
        score_left = F[i - 1][j]
        # print(score_up, score_left, score_diag)
        if score == score_diag + w[a[i - 1]][b[j - 1]]:
            alignment_a = a[i - 1] + alignment_a
            alignment_b = b[j - 1] + alignment_b
            i -= 1
            j -= 1
        elif score == score_left + d:
            alignment_a = a[i - 1] + alignment_a
            alignment_b = '-' + alignment_b
            i -= 1
        elif score == score_up + d:
            alignment_a = '-' + alignment_a
            alignment_b = b[j - 1] + alignment_b
            j -= 1
        elif score == 0:
            alignment_a = '-' * i + alignment_a
            alignment_b = '-' * j + alignment_b
            i -= 1
            j -= 1
    while i > 0:
        # print('2st while -> i: {}, j:{}'.format(i, j))
        alignment_a = a[i - 1] + alignment_a
        alignment_b = '-' + alignment_b
        i -= 1
    while j > 0:
        # print('3st while -> i: {}, j:{}'.format(i, j))
        alignment_a = '-' + alignment_a
        alignment_b = b[j - 1] + alignment_b
        j -= 1
    return alignment_a, alignment_b


def GlobalAlignment(a, b, w, d, F):
    # Input: strings a,b, matrix of weights w, weight of gap d, their matrix of global alignment F, computed by nw_matrix
    # Output: their global alignment
    alignment_a = ''
    alignment_b = ''
    i, j = len(a), len(b)
    while i > 0 or j > 0:
        score = F[i][j]
        score_diag = F[i - 1][j - 1]
        score_up = F[i][j - 1]
        score_left = F[i - 1][j]
        if score == score_diag + w[a[i - 1]][b[j - 1]]:
            alignment_a = a[i - 1] + alignment_a
            alignment_b = b[j - 1] + alignment_b
            i -= 1
            j -= 1
        elif score == score_left + d:
            alignment_a = a[i - 1] + alignment_a
            alignment_b = '-' + alignment_b
            i -= 1
        elif score == score_up + d:
            alignment_a = '-' + alignment_a
            alignment_b = b[j - 1] + alignment_b
            j -= 1
    while i > 0:
        alignment_a = a[i - 1] + alignment_a
        alignment_b = '-' + alignment_b
        i -= 1
    while j > 0:
        alignment_a = '-' + alignment_a
        alignment_b = b[j - 1] + alignment_b
        j -= 1
    return alignment_a, alignment_b


def align_info(str1, str2, w, d):
    # Input: strings str1, str2, matrix of weights w, weight of gap d
    # Output: info about alignment: alignment itself, strings themselves
    global_align = global_alignment_matrix(str1, str2, w, d)
    print('computed global alignment matrix')
    local_align = local_alignment_matrix(str1, str2, w, d)
    print('computed local alignment matrix')
    a = GlobalAlignment(str1, str2, w, d, global_align)
    print('constructed global alignment')
    b = LocalAlignment(str1, str2, w, d, local_align)
    print('constructed local alignment')
    score_global = global_align[len(str1)][len(str2)]
    score_local = local_align[len(str1)][len(str2)]
    global_alignment_string = a[0] + '\n' + similar_parts(a[0], a[1]) + '\n' + a[1]
    local_alignment_string = b[0] + '\n' + similar_parts(b[0], b[1]) + '\n' + b[1]
    print("Given sequences:\n{}\n{}\n".format(str1, str2))
    print('Global score is {}\n{}\n\nLocal score is {}\n{}'.format(score_global, global_alignment_string, score_local,
                                                                   local_alignment_string))


def LCSBacktrack(v: str, w: str) -> np.ndarray:
    # Input: string v,w
    # Output: compute their backtrack list in the alignment graph
    # Complexity: O(nm)
    n, m = len(v), len(w)
    s = [[0] * (m + 1) for i in range(n + 1)]
    backtrack = np.zeros((n + 1, m + 1, 3))
    # backtrack_ij contains boolean vector (a,b,c) which means existence of path (down, right, diag)
    for i in range(n + 1):
        s[i][0] = 0
        backtrack[i, 0, 0] = 1
    for j in range(m + 1):
        s[0][j] = 0
        backtrack[0, j, 1] = 1
    backtrack[0, 0, :] = [0, 0, 0]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = 1 if v[i - 1] == w[j - 1] else 0
            s[i][j] = max(s[i][j - 1], s[i - 1][j], s[i - 1][j - 1] + match)
            if s[i][j] == s[i - 1][j]:
                backtrack[i, j, 0] = 1

            if s[i][j] == s[i][j - 1]:
                backtrack[i, j, 1] = 1
            if s[i][j] == s[i - 1][j - 1] + match:
                backtrack[i, j, 2] = 1
    print_array(backtrack)
    return backtrack


# 5C Code challenge
def OutputLCS(backtrack: list, v: str, i: int, j: int) -> str:
    # Input: backtrack list, string v, positions i, j
    # Output: backtrack how we get there
    if i == 0 or j == 0:
        return ''
    if backtrack[i][j] == 'down':
        return OutputLCS(backtrack, v, i - 1, j)
    elif backtrack[i][j] == 'right':
        return OutputLCS(backtrack, v, i, j - 1)
    else:
        return OutputLCS(backtrack, v, i - 1, j - 1) + v[i - 1]


def Output_LCS(backtrack: np.ndarray, v: str, i: int, j: int) -> list:
    # TODO: rewrite it for np.ndarray
    if i == 0 or j == 0:
        return ['']
    s_down = []
    s_right = []
    s_diag = []
    if backtrack[i, j, 0] == 1:
        # print('i: {}, j: {},  k: {} -> {}'.format(i, j, 0, 1))
        s_down = Output_LCS(backtrack, v, i - 1, j)
    elif backtrack[i, j, 1] == 1:
        # print('i: {}, j: {},  k: {} -> {}'.format(i, j, 1, 1))
        s_right = Output_LCS(backtrack, v, i, j - 1)
    elif backtrack[i, j, 2] == 1:
        # print('i: {}, j: {},  k: {} -> {}'.format(i, j, 2, 1))
        s_diag = Output_LCS(backtrack, v, i - 1, j - 1)
        for h in range(len(s_diag)):
            s_diag[h] += v[i - 1]
    return s_down + s_right + s_diag


def LCS(v: str, w: str) -> str:
    backtrack = LCSBacktrack(v, w)
    a = Output_LCS(backtrack, v, len(v), len(w))
    for a_ in a:
        print(a)
