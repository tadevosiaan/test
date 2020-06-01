import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
from tqdm import tqdm


# CHAPTER 1. Where in the genome doest DNA replication begin?

def FASTA_to_lists(filename) -> tuple:
    seq_names = []
    seqs = []
    with open('{}'.format(filename), mode='r', encoding='utf-8') as f:
        buf = f.readline().rstrip()
        while buf:
            seq_name, seq = buf[1:], ''
            buf = f.readline().rstrip()
            while not buf.startswith('>') and buf:
                seq = seq + buf
                buf = f.readline().rstrip()
            seq_names.append(seq_name)
            seqs.append(seq)
    return seq_names, seqs


# 1A Code challenge
def PatternCount(Text: str, Pattern: str) -> int:
    # Input: a string Text and a k-mer Pattern
    # Output: find number of occurrences (with overlapping) of Pattern in Text
    # Assuming: len(Pattern)<<len(Text)
    # Complexity: O(n * k), so it's better to use KMP in linear time
    count = 0
    for i in range(0, len(Text) - len(Pattern) + 1):
        t = Text[i: i + len(Pattern)]
        if hash(t) == hash(Pattern):
            if t == Pattern:
                count += 1
    return count


# 1B Code challenge
def FrequentWords(Text: str, k: int) -> tuple:
    # Input: a string Text and an integer k;
    # Output: All most frequent k-mers in Text.
    # Complexity: O(n^2 * k)
    FrequentPatterns = []
    CountArray = []
    for i in range(0, len(Text) - k + 1):
        Pattern = Text[i: i + k]
        CountArray.append(PatternCount(Text=Text, Pattern=Pattern))
    maxCount = max(CountArray)  # calculating argmax
    for i in range(0, len(Text) - k + 1):
        if CountArray[i] == maxCount:
            FrequentPatterns.append(Text[i: i + k])
    return list(FrequentPatterns), maxCount


# 1C Code challenge
def ReverseComplement(Pattern: str) -> str:
    # Input: a DNA string Pattern
    # Output: the reverse complement of Pattern
    # Complexity: O(n)
    return Pattern[::-1].translate(str.maketrans('ACGT', 'TGCA'))


# 1D Code challenge
def PatternMatching(Pattern: str, Genome: str) -> list:
    from test.bioinfo.bio import KMP
    # Input: string Pattern and Genome
    # Output: all starting positions in Genome where Patterns appears as a substring
    # Complexity: using KMP we get O(|Genome| + |Pattern|) time and similar space
    return KMP(text=Genome, pattern=Pattern)


# 1E Code challenge
def ClumpFindingProblem(Genome: str, k: int, L: int, t: int) -> set:
    # Input: string Genome, integers k,L,t
    # Output: all distinct k-mers forming (L,t)-clumps in Genome
    # Complexity: O(L^2 k |Genome|) if bad implementation of FrequentWords, like O(L^2 k)
    candidates = set()
    for i in range(len(Genome) - L):
        s = Genome[i: i + L]
        words = FrequentWords(Text=s, k=k)[0]  # index 0 means we get a list
        words = Counter(words)
        # print('i: {}, str: {}, words:{}'.format(i, s, words))
        words = list(words.items())
        for f in words:
            if f[1] >= t:
                candidates.add(f[0])
    return candidates


# 1F Code challenge
def skew(genome: str, k: int):
    # Input: str genome, int k
    # Output: plot the difference #G-#C for the first k nucleotides of genome
    difference = []
    if genome[0] == 'C':
        difference.append(-1)
    elif genome[0] == 'G':
        difference.append(1)
    else:
        difference.append(0)

    for i in tqdm(range(1, k)):
        if genome[i] == 'C':
            difference.append(difference[i - 1] - 1)
        elif genome[i] == 'G':
            difference.append(difference[i - 1] + 1)
        else:
            difference.append(difference[i - 1])
    argmin_ = np.argmin(difference)
    argmax_ = np.argmax(difference)
    min_ = difference[argmin_]
    max_ = difference[argmax_]
    plt.plot(difference, label='Min at k={}: {}\nMax at k={}: {}'.format(argmin_, min_,
                                                                         argmax_, max_))
    plt.legend()
    plt.show()


# 1G Code challenge
def HammingDistance(u: str, v: str) -> int:
    # Input: strings u,v
    # Output: their Hamming distance
    distance = 0
    assert len(u) == len(v)
    for i in range(len(u)):
        if u[i] != v[i]:
            distance += 1
    return distance
