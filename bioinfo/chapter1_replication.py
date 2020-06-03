import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
from tqdm import tqdm


# CHAPTER 1. Where in the genome doest DNA replication begin?
# Algorithmic warmup

def ReverseComplement(Pattern: str) -> str:
    # Input: a DNA string Pattern
    # Output: the reverse complement of Pattern
    # Complexity: O(n)
    return Pattern[::-1].translate(str.maketrans('ACGT', 'TGCA'))


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
    # Complexity: O(nk), so it's better to use KMP in linear time
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
    # Complexity: O(n^2 k)
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
    # Complexity: using KMP we get O(n + m) time and similar space
    return KMP(text=Genome, pattern=Pattern)


# 1E Code challenge
def ClumpFindingProblem(Genome: str, k: int, L: int, t: int) -> set:
    # Input: string Genome, integers k,L,t
    # Output: all distinct k-mers forming (L,t)-clumps in Genome
    # Complexity: O(L^2 k n) if bad implementation of FrequentWords, like O(L^2 k)
    candidates = set()
    for i in range(len(Genome) - L + 1):
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
def skew_gc(genome: str, k: int, save=None):
    # Input: str genome, int k
    # Output: plot the difference #G-#C for the first k nucleotides of genome, GC content diagram
    # Complexity: O(min(n,k))
    difference = []
    gc_content = []
    if genome[0] == 'C':
        difference.append(-1)
        gc_content.append(1)
    elif genome[0] == 'G':
        difference.append(1)
        gc_content.append(1)
    else:
        difference.append(0)
        gc_content.append(0)

    for i in tqdm(range(1, k)):
        delta = 1 if genome[i] == 'G' or genome[i] == 'C' else 0
        if genome[i] == 'C':
            difference.append(difference[i - 1] - 1)
        elif genome[i] == 'G':
            difference.append(difference[i - 1] + 1)
        else:
            difference.append(difference[i - 1])
        gc_content.append((gc_content[i - 1] * i + delta) / (i + 1))
    gc_content = np.array(gc_content)
    gc_content_mean = gc_content.mean()
    gc_content_std = gc_content.std()
    argmin_ = np.argmin(difference)
    argmax_ = np.argmax(difference)
    min_ = difference[argmin_]
    max_ = difference[argmax_]
    fig, (ax1, ax2) = plt.subplots(2)
    print('Skew diagram:\nMin at k={}: {}\nMax at k={}: {}'.format(argmin_, min_, argmax_, max_))
    print('GC content:\nMean: {}\nSTD: {}'.format(gc_content_mean, gc_content_std))
    ax1.plot(difference, label='Skew diagram')
    ax1.plot([argmin_, argmax_], [min_, max_], 'rx')
    ax1.legend(loc='lower left')
    ax2.set_ylim([0, 1])
    ax2.plot(gc_content, label='GC content', color='green')
    ax2.legend(loc='lower left')
    if save is not None:
        plt.savefig('{}.png'.format(save))
    plt.show()


def plot_EColi():
    a, b = FASTA_to_lists('EColi')
    skew_gc(b[0], len(b[0]))


# 1G Code challenge
def HammingDistance(u: str, v: str) -> int:
    # Input: strings u,v
    # Output: their Hamming distance
    # Complexity: O(n)
    distance = 0
    assert len(u) == len(v)
    for i in range(len(u)):
        if u[i] != v[i]:
            distance += 1
    return distance


# 1H Code challenge
def ApproximatePatternMatching(Text: str, Pattern: str, d: int, print_mode=False) -> list:
    # Input: strings Pattern, Text, integer d
    # Output: all starting positions in Text where Pattern appears as a substring with at most d mismatches
    # Complexity: O(nm)
    pos = []
    for i in range(len(Text) - len(Pattern) + 1):
        s = Text[i: i + len(Pattern)]
        if print_mode:
            print('iter: {}, s: {}, hamdist: {}'.format(i, s, HammingDistance(s, Pattern)))
        if HammingDistance(s, Pattern) <= d:
            pos.append(i)
    return pos


def ApproximatePatternCount(Text: str, Pattern: str, d: int, print_mode=False) -> int:
    # Input: strings Pattern, Text, integer d
    # Output: count of positions in Text where Pattern appears as a substring with at most d mismatches
    # Complexity: O(nm)
    return len(ApproximatePatternMatching(Text, Pattern, d, print_mode))


def ComputeFrequencyArray(Text: str, k: int) -> Counter:
    # Input: string Text, k integer
    # Output: frequency array for k-mers
    # Complexity: O(n)
    c = Counter()
    for i in range(len(Text) - k + 1):
        c[Text[i:i + k]] += 1
    return c


def ImmediateNeighbors(Pattern: str) -> set:
    # Input: Pattern string
    # Output: get 1-neighborhood if the Pattern in HammDist
    # TODO: check
    neighborhood = {Pattern}
    for i in range(len(Pattern)):
        s = Pattern[i]
        for c in ['A', 'T', 'G', 'C']:
            if c != s:
                p = list(Pattern)
                p[i] = c
                neighbor = ''.join(p)
                neighborhood.add(neighbor)
    return neighborhood


# 1N Code challenge
def Neighbors(Pattern: str, d: int) -> set:
    # Input: Pattern string, d int
    # Output: d-neighborhood of Pattern in HammDist
    if d == 0:
        return {Pattern}
    if len(Pattern) == 1:
        return {'A', 'T', 'G', 'C'}
    neighborhood = set()
    suffix_neighbors = Neighbors(Pattern[1:], d)
    for Text in suffix_neighbors:
        if HammingDistance(Pattern[1:], Text) < d:
            for x in ['A', 'T', 'G', 'C']:
                neighborhood.add(x + Text)
        else:
            neighborhood.add(Pattern[0] + Text)
    return neighborhood


def IterativeNeighbors(Pattern: str, d: int) -> set:
    # Input: Pattern string, d int
    # Output: iterative version of Neighbors()
    # TODO: check
    neighborhood = {Pattern}
    for j in range(1, d + 1):
        for p in neighborhood:
            imm_neighbors = ImmediateNeighbors(p)
            for elem in imm_neighbors:
                neighborhood.add(elem)
    return neighborhood


# 1I Code challenge
def FrequentWordsWithMismatches(Text: str, k: int, d: int) -> list:
    # Input:  string Text, integers k,d
    # Output: all most frequent k-mers with up to d mismatches in Text
    # Complexity: slow
    c = Counter()
    for i in range(len(Text) - k + 1):
        s = Text[i:i + k]
        neighbors = Neighbors(s, d)
        for n in neighbors:
            c[n] += 1
    c = list(c.items())
    c.sort(key=lambda x: -x[1])
    cmax = c[0][1]
    for i in range(1, len(c)):
        if c[i][1] < cmax:
            c = c[:i]
            break
    c.sort(key=lambda x: x[0])
    return c


def ApproximatePatternCount_viaCounter(Text: str, Pattern: str, d: int) -> list:
    # Input: strings Pattern, Text, integer d
    # Output: count of positions in Text where Pattern appears as a substring with at most d mismatches
    # Complexity: O(nm)
    # TODO: it does something different than the output's expectation
    c = Counter()
    for i in range(len(Text) - len(Pattern) + 1):
        s = Text[i: i + len(Pattern)]
        neighbors = Neighbors(s, d)
        for n in neighbors:
            c[n] += 1
    c = list(c.items())
    return len(c), sorted(c, key=lambda x: -x[1])


# 1J Code challenge
def FrequentWordsWithMismatchesAndReverseComplement(Text: str, k: int, d: int) -> list:
    # Input:  a DNA string Text, integers k,d
    # Output: all k-mers Pattern that maximize ApproximatePatternCount(Text,Pattern,d)
    # + ApproximatePatternCount(Text,rc(Pattern),d) over all possible k-mersx
    # Complexity: slow
    candidats = []
    for i in range(len(Text) - k + 1):
        s = Text[i: i + k]
        neighbors = Neighbors(s, d)
        candidats += neighbors
    score = lambda s: ApproximatePatternCount(Text, s, d) + ApproximatePatternCount(Text, ReverseComplement(s), d)
    candidats = list(set(candidats))
    candidats_score = list(map(score, candidats))
    patterns = list(zip(candidats, candidats_score))
    patterns.sort(key=lambda x: -x[1])
    score_max = patterns[0][1]
    for i in range(1, len(patterns)):
        if patterns[i][1] < score_max:
            patterns = patterns[:i]
            break
    patterns.sort(key=lambda x: x[0])
    return patterns


def gen_DNA(n: int, k: int) -> list:
    # Input: n, k integers
    # Output: n DNA strings each len = k
    from random import choice
    l = []
    for i in range(n):
        s = ''
        for j in range(k):
            s += choice(['A', 'G', 'T', 'C'])
        l.append(s)
    return l


def experiment():
    np.random.seed(1337)
    a = gen_DNA(10, 20)
    print(a)
    s = 'AAAGGG'
    from random import choice

    for i in range(len(a)):
        j = np.random.randint(0, len(a[0]) - len(s) + 1)
        k = np.random.randint(0, 6)
        l = np.random.randint(0, 6)
        t = list(a[i])
        t_str = list(s)
        t_str[k] = choice(['A', 'G', 'T', 'C'])
        t_str[l] = choice(['A', 'G', 'T', 'C'])
        t_str = ''.join(t_str)
        t[j: j + len(s)] = t_str
        print('iter: {}, s: {}, s_new: {}, pos: {}'.format(i, s, t_str, j))
        a[i] = ''.join(t)
    print(a)
    long_Str = ''.join(a)
    print(long_Str)
    print(FrequentWordsWithMismatches(long_Str, 6, 2))
