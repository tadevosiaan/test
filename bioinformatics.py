from collections import Counter
import constants
import numpy as np


def NucleotidesCount(s: str) -> tuple:
    # Input: a DNA of letters A, T, C, G
    # Output: counts for each letter in the form (AGCT)
    # at most O(n)?
    return s.count('A'), s.count('G'), s.count('C'), s.count('T')

def FASTA_to_lists()->tuple:
    seq_names = []
    seqs = []
    with open('input.txt', mode='r', encoding='utf-8') as f:
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

def DNA_to_mRNA(s: str) -> str:
    # Input: string of DNA
    # Output: transcribe of DNA into RNA
    # at most O(n)

    return s.replace('T', 'U')


def PatternCount(Text, Pattern: str) -> int:
    # Input: a string Text and a k-mer Pattern
    # Output: find number of occurrences (with overlapping) of Pattern in Text
    # Assuming: len(Pattern)<<len(Text)
    # at most O(n * k)

    count = 0
    for i in range(0, len(Text) - len(Pattern) + 1):
        t = Text[i: i + len(Pattern)]
        if hash(t) == hash(Pattern):
            if t == Pattern:
                count += 1
    return count


def FrequentWords(Text: str, k: int) -> tuple:
    # Input: a string Text and an integer k;
    # Output: All most frequent k-mers in Text.
    # at most O(n^2 * k)

    FrequentPatters = set()
    CountArray = []
    for i in range(0, len(Text) - k + 1):
        Pattern = Text[i: i + k]
        CountArray.append(PatternCount(Text=Text, Pattern=Pattern))
    maxCount = max(CountArray)
    for i in range(0, len(Text) - k + 1):
        if CountArray[i] == maxCount:
            FrequentPatters.add(Text[i: i + k])
    return list(FrequentPatters), maxCount


def ReverseComplement(Pattern: str) -> str:
    # Input: a DNA string Pattern
    # Output: the reverse complement of Pattern
    # at most O(n)

    return Pattern[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def PatternMatching(Pattern, Genome: str) -> list:
    # Input: string Pattern and Genome
    # Output: all starting positions in Genome where Patterns appears as a substring

    pass


def RabbitsProblem(n, k: int) -> int:
    # Input: n month, k pairs new children after sex, 1 month is reproductive age
    # Input: number of alive rabbit pairs after n months, assuming no rabbit dies

    F = [1, 1]
    for i in range(n - 2):
        t = k * F[0] + F[1]
        F[0], F[1] = F[1], t
    return F[1]


def GC_content(s: str) -> float:
    # Input: DNA string
    # Output: its GC content, ie the percentage of C,G in DNA

    return (s.count('G') + s.count('C')) / len(s)


def Hamming(s1, s2: str) -> int:
    # Input: DNA strings s1,s2 of same len
    # Output: hamming distance s1,s2

    return sum(a != b for a, b in zip(s1, s2))


def GC_content():
    # Input: FASTA-encoded DNAs as text in input file(>Rosalind_xxxx, DNA)
    # Output: id of DNA with highest GC_content

    with open('input.txt', mode='r', encoding='utf-8') as f:
        max_gc_name, max_gc_content = '', 0
        buf = f.readline().rstrip()
        while buf:
            seq_name, seq = buf[1:], ''
            buf = f.readline().rstrip()
            while not buf.startswith('>') and buf:
                seq = seq + buf
                buf = f.readline().rstrip()
            seq_gc_content = GC_content(seq)
            if seq_gc_content > max_gc_content:
                max_gc_name, max_gc_content = seq_name, seq_gc_content
        print('{}\n{}'.format(max_gc_name, max_gc_content * 100))


def MendelsFirstLaw(k, m, n: int) -> float:
    # Input: int k, m, n > 0 representing a population containing k+m+n organisms:
    # k individuals are homozygous dominant for a factor, m are heterozygous, and n are homozygous recessive.
    # Output: The probability that two randomly selected mating
    # organisms will produce an individual possessing a dominant allele

    pAB = [[(k - 1) * k, k * m, k * n],
           [m * k, (m - 1) * m, m * n],
           [n * k, n * m, (n - 1) * n]]
    norm = 1.0 / ((m + n + k) * (m + n + k - 1))

    p_yyAB = [[0, 0, 0],
              [0, 0.25, 0.5],
              [0, 0.5, 1]]
    p_yyAB = np.array(p_yyAB)
    pAB = np.array(pAB) * norm
    # (pAB)ij = prob(A = i, B = j), where i=0 <-> YY, i=1 <-> Yy, i=2 <-> yy
    # (p_yyAB)ij = prob(yy | A=i, B=j)
    return 1 - np.sum(p_yyAB * pAB)


def RNA_to_Protein(s: str) -> str:
    # Input: mRNA string s
    # Output: transcribe of s into protein chains, aka polypeptide

    proteins = \
        {'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V',
         'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V',
         'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V',
         'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V',
         'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A',
         'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
         'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
         'UCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
         'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D',
         'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
         'UAA': 'Z', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
         'UAG': 'Z', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
         'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G',
         'UGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
         'UGA': 'Z', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
         'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}
    # there is no protein Z avaliable, here 'Z' means 'STOP'
    polypeptide = ''
    assert len(s) % 3 == 0
    assert s[:3] == 'AUG'
    assert s[-3:] in ['UAA', 'UAG', 'UGA']
    # assume there is no stop codon at the middle

    for i in range((len(s) // 3) - 1):
        polypeptide += proteins[s[3 * i: 3 * i + 3]]
    return polypeptide


def prefix_function(s: str) -> list:
    # Input: string s
    # Output: its prefix function, ie prefix_function[i]=len of the longest suffix equal to prefix at string s[:i]
    # at most O(n)

    pi = [0]
    for i in range(1, len(s)):
        j = pi[i - 1]
        while j > 0 and s[i] != s[j]:
            j = pi[j - 1]
        if s[i] == s[j]:
            j += 1
        pi.append(j)
    return pi


def KMP(Text, Pattern: str) -> list:
    # Knuth–Morris–Pratt algorithm
    # Input: strings Text, Pattern, len = n, m
    # Output: search for occurrences of a Pattern in Text (return indexies)
    # at most O(n + m) time, O(n+m) (can be O(m), because PF is online algorithm)

    w = Pattern + '#' + Text
    pf = prefix_function(w)
    print(w)
    print(pf)
    positions = []
    for i in range(len(Pattern), len(w)):
        if pf[i] == len(Pattern):
            positions.append(i - 2 * len(Pattern))
    return positions


def poly_hash(s: str, p: int) -> int:
    # Input: string s, integer p
    # Output: compute polynomial hash of s with power p
    n = len(s)
    poly = np.array([p ** i for i in range(n)])
    s = np.array([ord(s[i]) for i in range(n)])
    return np.sum(poly * s)


def RK(t, s: str, p=31) -> list:
    # Rabin-Karp algorithm
    # Input: text t , pattern s, len = n, m; power value p for hash
    # Output: search for occurrences of a Pattern in Text (return indexies)
    # at most O(n+m) time

    positions = []
    n = len(t)
    k = len(s)
    print(n,k)
    hash_pattern = poly_hash(s, p)
    hash_text = [poly_hash(t[:k], p)]
    # computing all the hashes of len k in text t
    for i in range(0, n - k):
        hash_text.append((hash_text[i] + p**k * ord(t[i+k]) - ord(t[i]))/p)
    # checking hashes
    for i in range(n - k + 1):
        if hash_text[i] == hash_pattern:
            if t[i:i + k] == s:
                positions.append(i)
    return positions


def Consensus():
    # Input: FASTA-encoded DNAs as text in input file(>Rosalind_xxxx, DNA)
    # Output: consesus string computed by profile matrix

    seq_names, seqs = FASTA_to_lists()
    # A <-> 0, C <-> 1, G <-> 2, T <-> 3
    n = len(seqs[0]) #len of dna
    seqs = list(map(lambda x: x.translate(str.maketrans('ACGT', '0123')), seqs))
    seqs = list(map(lambda x: list(x), seqs))
    seqs = list(map(lambda x: list(map(int, x)), seqs))
    seqs = np.array(seqs)
    profile = np.zeros((4, n))
    for i in range(4):
        for j in range(n):
            profile[i, j] = np.sum(seqs[:,j] == i)
    cons = ''
    matrix_string = ''
    for i in range(n):
        cons += str(np.argmax(profile[:, i]))
    cons = cons.translate(str.maketrans('0123', 'ACGT'))
    for i in range(4):
        matrix_string += str(i).translate(str.maketrans('0123','ACGT')) + ': '
        for j in range(n - 1):
            matrix_string += str(int(profile[i, j]))+' '
        matrix_string += str(int(profile[i, n - 1]))+'\n'
    print(cons+'\n'+matrix_string)


def MortalFib(n, m: int)->int:
    # Input: n months,m months to be alive int, n<=100, m<=20
    # Output: number of pairs after n months

    x = np.zeros(m, dtype=np.uint64)
    x[0] = 1
    Q = np.zeros((m,m), dtype=np.uint64)
    Q[0, 1:] = np.ones(m-1)
    for i in range(1, m):
        Q[i, i - 1] = 1
    for i in range(n - 1):
        x = np.dot(Q, x)
    print(np.sum(x))

class OverlapGraph:
    #Overlapping graph for set of strings DNA

    def __init__(self, k:int, seq:dict):
        self.k = k
        # vertices
        self.v = seq
        # edges
        self.e = self.calculate_edges()

    def calculate_edges(self):
        n = len(self.v)
        keys = list(self.v.keys())
        z = np.zeros((n,n))
        for i in range(n):
            for j in range(n):
                z[i, j] = self.v[keys[i]][-self.k:] == self.v[keys[j]][:self.k]
        for i in range(n):
            z[i, i] = 0
        return z

    def print_RosalindNames_vertices(self):
        x, y = np.nonzero(self.e) #positions where True
        vkeys = list(self.v.keys())
        for i in range(len(x)):
            print(vkeys[x[i]], vkeys[y[i]])

    def print_adjacency_matrix(self):
        print(self.e)


def ExpectedOffspring(a) -> float:
    # Input: 6 integer numbers a_1 a_2 a_3 a_4 a_5 a_6 as list, representing the amount of pairs
    # AA-AA
    # AA-Aa
    # AA-aa
    # Aa-Aa
    # Aa-aa
    # aa-aa
    # Output: expected number of offspring displaying the dominant phenotype in the next generation,
    # under the assumption that every couple has exactly two offspring.
    EY = np.array([1, 1, 1, 0.75, 0.5, 0])
    EX = np.sum(2 * EY * np.array(a))
    return EX


def main():
    pass

main()
