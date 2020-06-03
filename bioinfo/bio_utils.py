from collections import Counter
import numpy as np
from tqdm import tqdm
import time
from random import choice
from bio_constants import Proteins, Mass, Mass_int


def NucleotidesCount(s: str) -> tuple:
    # Input: a DNA of letters A, T, C, G
    # Output: counts for each letter in the form (AGCT)
    # Complexity: O(n)?
    return s.count('A'), s.count('G'), s.count('C'), s.count('T')


def HammingDistance(u: str, v: str) -> int:
    # Input: strings u, v; u may be shorter than v
    # Output: hamming distance between u and v
    assert len(u) <= len(v)
    if len(u) == len(v):
        dist = sum(a != b for a, b in zip(s1, s2))
    else:
        dist = float('inf')
        for i in range(len(v) - len(u) + 1):
            s = v[i: i + len(u)]
            d = HammingDistance(s, u)
            if d < dist:
                dist = d
    return dist


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


def gen_DNA(n: int, k: int) -> list:
    # Input: n, k integers
    # Output: n DNA strings each len = k
    l = []
    for i in range(n):
        s = ''
        for j in range(k):
            s += choice(['A', 'G', 'T', 'C'])
        l.append(s)
    return l


def gen_kMers(k: int) -> list:
    # Input: integer k
    # Output: all k-mers over alphabet {A,T,G,C}
    if k == 1:
        return ['A', 'T', 'G', 'C']
    else:
        kmers = gen_kMers(k - 1)
        a = list(map(lambda x: 'A' + x, kmers))
        t = list(map(lambda x: 'T' + x, kmers))
        c = list(map(lambda x: 'C' + x, kmers))
        g = list(map(lambda x: 'G' + x, kmers))
        return a + t + c + g


def DNA_to_RNA(s: str) -> str:
    # Input: string of DNA
    # Output: transcribe of DNA into RNA
    # Complexity: O(n)
    return s.replace('T', 'U')


def similar_parts(s1: str, s2: str) -> str:
    # Input: strings s1, s2
    # Output: string s, which matches both s1, s2 at appropriate positions
    n = len(s1)
    s = ''
    for i in range(n):
        if s1[i] == s2[i]:
            s += s1[i]
        else:
            s += ' '
    return s


def print_matrix(a):
    # Input: matrix a
    # Output: print it
    for list_ in a:
        for l in list_:
            print(l, end=' ')
        print()


def GC_content_value(s: str) -> float:
    # Input: DNA string
    # Output: its GC content, ie the percentage of C,G in DNA
    return (s.count('G') + s.count('C')) / len(s)


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
            seq_gc_content = GC_content_value(seq)
            if seq_gc_content > max_gc_content:
                max_gc_name, max_gc_content = seq_name, seq_gc_content
        print('{}\n{}'.format(max_gc_name, max_gc_content * 100))


def prefix_function(s: str) -> list:
    # Input: string s
    # Output: its prefix function, ie prefix_function[i]=len of the longest suffix equal to prefix at string s[:i]
    # Complexity: O(n)
    pi = [0]
    for i in range(1, len(s)):
        j = pi[i - 1]
        while j > 0 and s[i] != s[j]:
            j = pi[j - 1]
        if s[i] == s[j]:
            j += 1
        pi.append(j)
    return pi


def KMP(text: str, pattern: str) -> list:
    # Knuth–Morris–Pratt algorithm
    # Input: strings text, pattern, len = n, m
    # Output: search for occurrences of a pattern in text (return indexies)
    # Complexity: O(n + m) time, O(n+m) (can be O(m), because PF is online algorithm)
    w = pattern + '#' + text
    pf = prefix_function(w)
    positions = []
    for i in range(len(pattern), len(w)):
        if pf[i] == len(pattern):
            positions.append(i - 2 * len(pattern))
    return positions


def poly_hash(s: str, p=51, m=2 ** 64) -> int:
    # Input: string s, integer p, modulo m
    # Output: compute polynomial hash of s with power p modulo m
    n = len(s)
    poly = np.array([p ** i for i in range(n)])
    s = np.array([ord(s[i]) for i in range(n)])
    return np.sum(poly * s) % m


def poly_hash_array(s: str, l: int, p=51, m=2 ** 64) -> list:
    # Input: string Text, parametres p,M for poly hash function, len of substting lenght in Text
    # Output: array of l-length hashes of string Text
    n = len(s)
    poly = np.array([p ** i for i in range(l)], dtype=object)
    s = np.array([ord(s[i]) for i in range(n)], dtype=object)
    hashes = []
    hashes.append(np.sum(poly * s[:l] % m))

    for i in range(0, n - l):
        h = ((hashes[i] % m + (p * poly[l - 1] * s[i + l]) % m - s[i]) / p) % m
        hashes.append(h)
    return hashes


def RK(text: str, pattern: str, p=51, m=2 ** 64) -> list:
    # Rabin-Karp algorithm
    # Input: text t , pattern s, len = n, m; power value p for hash
    # Output: search for occurrences of a Pattern in Text (return indexies)
    # Complexity: O(n+m)
    positions = []
    n = len(text)
    k = len(pattern)
    hash_pattern = poly_hash(pattern, p)
    hash_text = poly_hash_array(s=text, l=k, p=p, m=m)
    # checking hashes
    for i in range(n - k + 1):
        if hash_text[i] == hash_pattern:
            if text[i:i + k] == pattern:
                positions.append(i)
    return positions


def Consensus():
    # Input: FASTA-encoded DNAs as text in input file(>Rosalind_xxxx, DNA)
    # Output: consesus string computed by profile matrix
    seq_names, seqs = FASTA_to_lists()
    # A <-> 0, C <-> 1, G <-> 2, T <-> 3
    n = len(seqs[0])  # len of dna
    seqs = list(map(lambda x: x.translate(str.maketrans('ACGT', '0123')), seqs))
    seqs = list(map(lambda x: list(x), seqs))
    seqs = list(map(lambda x: list(map(int, x)), seqs))
    seqs = np.array(seqs)
    profile = np.zeros((4, n))
    for i in range(4):
        for j in range(n):
            profile[i, j] = np.sum(seqs[:, j] == i)
    cons = ''
    matrix_string = ''
    for i in range(n):
        cons += str(np.argmax(profile[:, i]))
    cons = cons.translate(str.maketrans('0123', 'ACGT'))
    for i in range(4):
        matrix_string += str(i).translate(str.maketrans('0123', 'ACGT')) + ': '
        for j in range(n - 1):
            matrix_string += str(int(profile[i, j])) + ' '
        matrix_string += str(int(profile[i, n - 1])) + '\n'
    print(cons + '\n' + matrix_string)


def common_substrings(u: str, v: str, l: int, p=51, m=2 ** 64) -> list:
    # Input: two strings u,v, integer l
    # Output: list of positions, where u,v have common l-substring
    if l == 0:
        return [None]  # special case to observe
    if not is_there_common_substring(u, v, l):
        return []
    u_hashes = sorted(enumerate(poly_hash_array(s=u, l=l, p=p, m=m)), key=lambda x: x[1])
    v_hashes = sorted(enumerate(poly_hash_array(s=v, l=l, p=p, m=m)), key=lambda x: x[1])
    i, j = 0, 0
    common = []
    while i < len(u_hashes) and j < len(v_hashes):
        if u_hashes[i][1] > v_hashes[j][1]:
            j += 1
        elif u_hashes[i][1] < v_hashes[j][1]:
            i += 1
        else:
            common.append((u_hashes[i][0], v_hashes[j][0]))
            if i + 1 < len(u_hashes) and j + 1 < len(v_hashes):
                if u_hashes[i + 1][1] == u_hashes[i][1] and v_hashes[j + 1][1] > v_hashes[j][1]:
                    i += 1
                elif u_hashes[i + 1][1] > u_hashes[i][1] and v_hashes[j + 1][1] == v_hashes[j][1]:
                    j += 1
                elif u_hashes[i + 1][1] == u_hashes[i][1] and v_hashes[j + 1][1] == v_hashes[j][1]:
                    common.append((u_hashes[i][0], v_hashes[j + 1][0]))
                    common.append((u_hashes[i + 1][0], v_hashes[j][0]))
                    i += 1
                    j += 1
                else:
                    i += 1
                    j += 1
            elif i == len(u_hashes) - 1 and j + 1 < len(v_hashes):
                j += 1
            elif j == len(v_hashes) - 1 and i + 1 < len(u_hashes):
                i += 1
            else:
                i += 1
                j += 1
    return common


def is_there_common_substring(u: str, v: str, l: int, p=51, m=2 ** 64) -> bool:
    # Input: two strings u,v, integer l
    # Output: True, if there exists common l-substring, False either
    if l == 0:
        return True
    if l >= min(len(u), len(v)) + 1:
        return False
    u_hashes = sorted(enumerate(poly_hash_array(s=u, l=l, p=p, m=m)), key=lambda x: x[1])
    v_hashes = sorted(enumerate(poly_hash_array(s=v, l=l, p=p, m=m)), key=lambda x: x[1])
    i, j = 0, 0
    while i < len(u_hashes) and j < len(v_hashes):
        if u_hashes[i][1] > v_hashes[j][1]:
            j += 1
        elif u_hashes[i][1] < v_hashes[j][1]:
            i += 1
        else:
            return True
    return False


def RNA_from_Protein(s: str) -> int:
    # Input: pepdite chain s
    # Output: the number of different RNA chains encoding that polypeptide
    def f(p: str) -> int:
        # Input: protein p as letter of english alphabet
        # Output: amount of different RNA string which encode that protein
        cnt = 0
        for k in Proteins.keys():
            if Proteins[k] == p:
                cnt += 1
        return cnt

    s = s + 'Z'  # added stop codon Z
    return reduce(lambda x, y: x * y, [f(v) for v in s])


def longest_common_substrings(u: str, v: str) -> list:
    # Input: u, v string
    # Output: their longest common substrings
    def f(u: str, v: str, l_1: int, l_2: int) -> int:
        # Input: u,v, l_1, l_2
        # Output: find the lenght of longest common substring
        if l_2 - l_1 == 1:
            return l_1
        k = (l_1 + l_2) // 2
        if is_there_common_substring(u, v, k):
            return f(u, v, k, l_2)
        else:
            return f(u, v, l_1, k)

    k = f(u, v, 0, min(len(u), len(v)) + 1)
    common_substrings_list = common_substrings(u, v, k)
    common = []
    for t in common_substrings_list:
        common.append(u[t[0]:t[0] + k])
    return common


def CalculateProteinMass(p: str, mass_mode='float') -> float:
    # Input: protein p
    # Output: its mass in Da
    m = 0
    if mass_mode == 'int':
        for p_ in p:
            m += Mass_int[p_]
    elif mass_mode == 'float':
        for p_ in p:
            m += Mass[p_]
    return m


def longest_polyndrome(s: str) -> list:
    # Input: string s
    # Output: longest substring of s that is polyndrome(as list)
    def exists_polyndrome(w: str, length: int) -> bool:
        # Input: string w, length
        # Output: check is there exists a polyndrome of length in w
        x = sorted(enumerate(poly_hash_array(s=w, l=length)), key=lambda x: x[1])
        y = sorted(enumerate(poly_hash_array(s=w[::-1], l=length)), key=lambda x: x[1])
        i, j = 0, 0
        while i < len(x) and j < len(y):
            if x[i][1] > y[j][1]:
                j += 1
            elif x[i][1] < y[j][1]:
                i += 1
            else:
                return True
        return False

    def bin_search(s: str, l_1: int, l_2: int) -> int:
        # Input: string s, l_1, l_2 borders
        # Output: length of the longest polyndrome inbetween l_1, l_2

        # print("current (l1, l2) = ({},{})".format(l_1, l_2))
        if l_2 - l_1 == 1:
            return l_1
        k = (l_1 + l_2) // 2
        if exists_polyndrome(w=s, length=k):
            return bin_search(s, k, l_2)
        else:
            return bin_search(s, l_1, k)

    k = bin_search(s, 0, len(s) + 1)
    common = common_substrings(u=s, v=s[::-1], l=k)  # common if looking at hashes
    polys = []
    seen = set()
    for obj in common:
        if s[obj[0]:obj[0] + k] == s[obj[0]:obj[0] + k][::-1] and not (obj[0], s[obj[0]:obj[0] + k]) in seen:
            polys.append((obj[0], s[obj[0]:obj[0] + k]))
            seen.add((obj[0], s[obj[0]:obj[0] + k]))
    return polys


if __name__ == '__main__':
    main()
