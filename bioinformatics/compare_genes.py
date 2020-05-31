from tqdm import tqdm
import numpy as np


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


def print_array(a):
    n = len(a)
    m = len(a[0])
    for i in range(n):
        for j in range(m):
            print(a[i][j], end=' ')
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


a = ManhattanTourist(Down, Right)


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
    #TODO: write it
    pass


def compute_score(b, G):
    # Input: node b, DAG G
    # Output: score s[b] for LongestPath
    pass


def LongestPathInDAG(G, source, sink):
    # Input: an edge-weighted directed acyclic graph with source and sink nodes
    # Output: a longest path from source to sink in the directed acyclic graph
    #TODO: rewrite it properly to work on arbitrary DA graphs
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

# ______________________________
# 8 From global to local alignment
# ______________________________


# ______________________________
# 9 Penalizing insertions and deletions in sequence alignment
# ______________________________


# ______________________________
# 10 Multiple sequence alignment
# ______________________________
