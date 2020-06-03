import numpy as np
from bio_utils import HammingDistance, gen_kMers
from chapter1_replication import Neighbors


# CHAPTER 2. Which DNA patterns play the role of molecular clocks?
# Randomized algorithms


def d(u: str, DNA: list) -> int:
    t = len(DNA)
    dist = 0
    for i in range(t):
        dist += HammingDistance(u, DNA[i])
    return dist


def Profile(M: list) -> list:
    # Input: list Motifs
    # Output: Profile matrix for these Motifs
    # TODO: add code
    pass


def Motifs(P: list, DNA: list) -> list:
    # Input: profile matrix P, list of strings DNA
    # Output: profle-most probable k-mers motifs
    # TODO: add code
    pass


def Score(M: list) -> int:
    # Input: list Motifs
    # Output: calculate score of Motifs
    pass


def iskdMotif(u: str, d: int, DNA: list) -> bool:
    # Input: string u of len k, int d, list of strings DNA
    # Output: check whether u is (k,d)-motif
    neighbors = Neighbors(u, d)
    for i in range(len(DNA)):
        flag = False
        for n in neighbors:
            if n in DNA[i]:
                flag = True
        if not flag:
            return False
    return True


def MotifEnumeration(DNA: list, k: int, d: int) -> set:
    motifs = set()
    for i in range(len(DNA)):
        for j in range(len(DNA[i]) - k + 1):
            s = DNA[i][j:j + k]
            neighbors = Neighbors(s, d)
            for n in neighbors:
                if iskdMotif(n, d, DNA):
                    motifs.add(n)
    return motifs


def MedianString(DNA: list, k: int) -> str:
    dist = float('-inf')
    kmers = gen_kMers(k)
    for Pattern in kmers:
        if dist > d(Pattern, DNA):
            dist = d(Pattern, DNA)
            median = Pattern
    return median


def RandomizedMotifSearch(DNA: list, k: int):
    # Input: list of strings DNA, integer k
    # Output: list of t k-mer motifs with best score
    t = len(DNA)
    m = []
    for i in range(t):
        j = np.random.randint(0, DNA[t] - k + 1)
        m.append(DNA[i][j:j + k])
    BestMotifs = m.copy()
    while True:
        p = Profile(m)
        m = Motifs(p, DNA)
        if Score(m) < Score(BestMotifs):
            BestMotifs = m
        else:
            return BestMotifs


def GreedyMotifSearch(DNA: list, k: int):
    # TODO: add code
    t = len(DNA)
    pass


def ProfileMatrix(cnt_matrix: np.ndarray) -> np.ndarray:
    # TODO: add code
    pass


def LaplaceRule(cnt_matrix: np.ndarray) -> np.ndarray:
    return (cnt_matrix + 1) / 8


def GibbsSampler(DNA: list, k: int, N: int):
    # Input: list of strings DNA, integer k, number of iterations N
    # Output: list of t k-mer motifs with best score
    t = len(DNA)
    m = []
    for i in range(t):
        j = np.random.randint(0, DNA[t] - k + 1)
        m.append(DNA[i][j:j + k])
    BestMotifs = m.copy()
    for i in range(N):
        j = np.random.randint(0, t)
        # TODO: write it further, motifs except j-th
        p = Profile(m)
        m = Motifs(p, DNA)
        if Score(m) < Score(BestMotifs):
            BestMotifs = m
        else:
            return BestMotifs
