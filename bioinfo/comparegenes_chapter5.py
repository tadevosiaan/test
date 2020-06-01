from tqdm import tqdm
import numpy as np
import random


# CHAPTER 5. How Do We Compare Genes?
# ______________________________
# 1 Introduction to sequence alignment
# ______________________________

# blablabla insertions deletions matches mismatches, sum of matches=score
def LongestCommonSubsequence(u: str, v: str) -> str:
    # Input: two strings
    # Output: a longest common subsequence of there strings
    pass


# ______________________________
# 2 The Manhattan tourist problem
# ______________________________

Down = [[0, 0, 0, 0, 0],
        [1, 0, 2, 4, 3],
        [4, 6, 5, 2, 1],
        [4, 4, 5, 2, 1],
        [5, 6, 8, 5, 3]]
Right = [[0, 3, 2, 4, 0],
         [0, 3, 2, 4, 2],
         [0, 0, 7, 3, 4],
         [0, 3, 3, 0, 2],
         [0, 1, 3, 2, 2]]
Diag = [[0, 0, 0, 0, 0],
        [0, 5, 0, 2, 1],
        [0, 8, 4, 3, 0],
        [0, 10, 8, 9, 5],
        [0, 5, 6, 4, 7]]


def print_matrix(a):
    for list_ in a:
        for l in list_:
            print(l, end=' ')
        print()


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
    print_array(s)
    print(path)
    return s[n][m], path


def LongestPathInDirectedGraph(graph):
    # Input: an edge-weighted directed graph with source and sink nodes
    # Output: a longest path from source to sink in the directed graph
    pass


# ______________________________
# 3 Sequence alignment is Manhattan tourist in disguise
# ______________________________


# ______________________________
# 4 An intro to dynamic programming: the change problem
# ______________________________


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


# ______________________________
# 5 The Manhattan tourist problem revisited
# ______________________________
"""
class Node:
    def __init__(self, x, y, v=0):
        self.x = x
        self.y = y
        self.v = v

    def __str__(self):
        return 'node: ({}, {})'.format(self.x, self.y)


class Edge:
    def __init__(self, u: Node, v: Node, w: float):
        self.tail = u
        self.head = v
        self.w = w

    def __str__(self):
        return '({}, {}) -> ({}, {}); w = {}'.format(self.tail.x, self.tail.y, self.head.x, self.head.y,
                                                     self.w)


class Grid:
    def __init__(self, n, m):
        self.source = Node(0, 0)
        self.sink = Node(n, m)
        self.nodes = []
        self.edges = []
        for i in range(self.sink.x + 1):
            for j in range(self.sink.y + 1):
                if Node(i, j) not in self.nodes:
                    self._add_node(Node(i, j))
        self.set_randomly()

    def _add_node(self, n: Node):
        self.nodes.append(n)

    def add_edge(self, ux: int, uy: int, vx: int, vy: int, w: float):
        u = Node(ux, uy)
        v = Node(vx, vy)
        e = Edge(u, v, w)
        self.edges.append(e)

    def print_nodes(self):
        for n in self.nodes:
            print(n)

    def print_edges(self):
        for e in self.edges:
            print(e)

    def set_randomly(self):
        for i in range(self.sink.x):
            for j in range(self.sink.y + 1):
                w = np.random.randint(0, 11)
                self.add_edge(i, j, i + 1, j, w)
        for j in range(self.sink.y):
            for i in range(self.sink.x):
                w = np.random.randint(0, 11)
                self.add_edge(i, j, i, j + 1, w)
    __
    def SouthOrEast(self, i, j):
        #Find longest path to node (i,j) from (0,0)
        p = Node(i, j)
        if p.x == 0 and p.y == 0:
            return 0
        a = float('-inf')
        b = float('-inf')
        if p.x > 0:
            a = SouthOrEast(p.x - 1, p.y) + w[i,j]
        if p.y > 0:
            b = SouthOrEast(p.x, p.y - 1) + w[i, j]
        return max(a, b)
"""


# shitshitshit


def hamming_distance(u: str, v: str) -> int:
    assert len(u) == len(v)
    distance = 0
    for i in range(len(u)):
        if u[i] != v[i]:
            distance += 1
    return distance


# ______________________________
# 6 From Manhattan to an arbitrary DirectedACyclicGraph (=DAG)
# ______________________________
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
        s[b] = compute_score(b, G)
    return s[sink]


# ______________________________
# 7 Backtracking in the alignment graph
# ______________________________

# >>>>>>>>>>we are here<<<<<<<<<
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
    # TODO: check if it works
    alignment_a = ''
    alignment_b = ''
    i, j = len(a), len(b)
    while i > 0 or j > 0:
        score = F[i][j]
        print('1st while -> i: {}, j:{}'.format(i, j))
        score_source = 0
        score_diag = F[i - 1][j - 1]
        score_up = F[i][j - 1]
        score_left = F[i - 1][j]
        print(score_up, score_left, score_diag)
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
        print('2st while -> i: {}, j:{}'.format(i, j))
        alignment_a = a[i - 1] + alignment_a
        alignment_b = '-' + alignment_b
        i -= 1
    while j > 0:
        print('3st while -> i: {}, j:{}'.format(i, j))
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


def similar_parts(s1, s2):
    n = len(s1)
    s = ''
    for i in range(n):
        if s1[i] == s2[i]:
            s += s1[i]
        else:
            s += ' '
    return s


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


w_possible = {
    'A': {'A': 10, 'G': -1, 'C': -3, 'T': -4},
    'G': {'A': -1, 'G': 7, 'C': -5, 'T': -3},
    'C': {'A': -3, 'G': -5, 'C': 9, 'T': 0},
    'T': {'A': -4, 'G': -3, 'C': 0, 'T': 8}
}
w_ordinary = {
    'A': {'A': 1, 'G': 0, 'C': 0, 'T': 0},
    'G': {'A': 0, 'G': 1, 'C': 0, 'T': 0},
    'C': {'A': 0, 'G': 0, 'C': 1, 'T': 0},
    'T': {'A': 0, 'G': 0, 'C': 0, 'T': 1}
}
w_mu1 = {
    'A': {'A': 1, 'G': -1, 'C': -1, 'T': -1},
    'G': {'A': -1, 'G': 1, 'C': -1, 'T': -1},
    'C': {'A': -1, 'G': -1, 'C': 1, 'T': -1},
    'T': {'A': -1, 'G': -1, 'C': -1, 'T': 1}
}


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
def OutputLCS(backtrack: list, v, i, j) -> str:
    if i == 0 or j == 0:
        return ''
    if backtrack[i][j] == 'down':
        return OutputLCS(backtrack, v, i - 1, j)
    elif backtrack[i][j] == 'right':
        return OutputLCS(backtrack, v, i, j - 1)
    else:
        return OutputLCS(backtrack, v, i - 1, j - 1) + v[i - 1]


def Output_LCS(backtrack: np.ndarray, v: str, i: int, j: int) -> list:
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


v = 'GCC-C-AGTC-TATGT-CAGGGGGCACG--A-GCATGCACA'.replace('-', '')
w = 'GCCGCC-GTCGT-T-TTCAG----CA-GTTATGT-T-CAGAT'.replace('-', '')
align_info(v, w, w_mu1, -1)

# ______________________________
# 8 From global to local alignment
# ______________________________


# ______________________________
# 9 Penalizing insertions and deletions in sequence alignment
# ______________________________


# ______________________________
# 10 Multiple sequence alignment
# ______________________________
from skbio import alignment


