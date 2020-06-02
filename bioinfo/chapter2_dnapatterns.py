import numpy as np


# CHAPTER 2. Which DNA patterns play the role of molecular clocks?
# Randomized algorithms

def Profile(M: list) -> list:
    # Input: list Motifs
    # Output: Profile matrix for these Motifs
    pass


def Motifs(P: list, DNA: list) -> list:
    # Input: profile matrix P, list of strings DNA
    # Output: profle-most probable k-mers motifs
    pass


def Score(M: list) -> int:
    # Input: list Motifs
    # Output: calculate score of Motifs
    pass


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
        # TODO: write it futher, motifs except j-th
        p = Profile(m)
        m = Motifs(p, DNA)
        if Score(m) < Score(BestMotifs):
            BestMotifs = m
        else:
            return BestMotifs
