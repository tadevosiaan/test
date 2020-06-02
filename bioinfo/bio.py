from collections import Counter
import numpy as np
from tqdm import tqdm
from itertools import permutations
import time


def NucleotidesCount(s: str) -> tuple:
    # Input: a DNA of letters A, T, C, G
    # Output: counts for each letter in the form (AGCT)
    # Complexity: O(n)?
    return s.count('A'), s.count('G'), s.count('C'), s.count('T')


def FASTA_to_lists(filename) -> tuple:
    seq_names = []
    seqs = []
    with open('{}.txt'.format(filename), mode='r', encoding='utf-8') as f:
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


def DNA_to_RNA(s: str) -> str:
    # Input: string of DNA
    # Output: transcribe of DNA into RNA
    # Complexity: O(n)
    return s.replace('T', 'U')


def RabbitsProblem(n: int, k: int) -> int:
    # Input: n month, k pairs new children after sex, 1 month is reproductive age
    # Output: number of alive rabbit pairs after n months, assuming no rabbit dies
    F = [1, 1]
    for i in range(n - 2):
        t = k * F[0] + F[1]
        F[0], F[1] = F[1], t
    return F[1]


def GC_content_value(s: str) -> float:
    # Input: DNA string
    # Output: its GC content, ie the percentage of C,G in DNA
    return (s.count('G') + s.count('C')) / len(s)


def Hamming(s1: str, s2: str) -> int:
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
            seq_gc_content = GC_content_value(seq)
            if seq_gc_content > max_gc_content:
                max_gc_name, max_gc_content = seq_name, seq_gc_content
        print('{}\n{}'.format(max_gc_name, max_gc_content * 100))


def MendelsFirstLaw(k: int, m: int, n: int) -> float:
    # Input: int k, m, n > 0 representing a population containing k+m+n organisms:
    # k individuals are homozygous dominant for a factor, m are heterozygous, and n are homozygous recessive.
    # Output: The probability that two randomly selected mating
    # organisms will produce an individual possessing a dominant allele
    p_AB = [[(k - 1) * k, k * m, k * n],
            [m * k, (m - 1) * m, m * n],
            [n * k, n * m, (n - 1) * n]]
    norm = 1.0 / ((m + n + k) * (m + n + k - 1))
    p_yyAB = [[0, 0, 0],
              [0, 0.25, 0.5],
              [0, 0.5, 1]]
    p_yyAB = np.array(p_yyAB)
    p_AB = np.array(p_AB) * norm
    # (p_AB)ij = prob(A = i, B = j), where i=0 <-> YY, i=1 <-> Yy, i=2 <-> yy
    # (p_yyAB)ij = prob(yy | A=i, B=j)
    return 1 - np.sum(p_yyAB * p_AB)


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


def MortalFib(n: int, m: int) -> int:
    # Input: n months,m months to be alive int, n<=100, m<=20
    # Output: number of pairs after n months
    x = np.zeros(m, dtype=np.uint64)
    x[0] = 1
    Q = np.zeros((m, m), dtype=np.uint64)
    Q[0, 1:] = np.ones(m - 1)
    for i in range(1, m):
        Q[i, i - 1] = 1
    for i in range(n - 1):
        x = np.dot(Q, x)
    print(np.sum(x))


class OverlapGraph:
    # Overlapping graph for set of strings DNA
    def __init__(self, k: int, seq: dict):
        self.k = k
        # vertices
        self.v = seq
        # edges
        self.e = self.calculate_edges()

    def calculate_edges(self):
        n = len(self.v)
        keys = list(self.v.keys())
        z = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                z[i, j] = self.v[keys[i]][-self.k:] == self.v[keys[j]][:self.k]
        for i in range(n):
            z[i, i] = 0
        return z

    def print_RosalindNames_vertices(self):
        x, y = np.nonzero(self.e)  # positions where True
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


def InferGenotypeFromPedigree(par1, par2):
    # node is encoded by list (x1, x2, x3), where x1=Prob(AA), x2=Prob(Aa), x3=Prob(aa)
    def get_kernel(par1, par2):
        P = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                P[i, j] = par1[i] * par2[j]
        return P

    C_AA = np.array([[1, 0.5, 0],
                     [0.5, 0.25, 0],
                     [0, 0, 0]])
    C_Aa = np.array([[0, 0.5, 1],
                     [0.5, 0.5, 0.5],
                     [1, 0.5, 0]])
    C_aa = np.array([[0, 0, 0],
                     [0, 0.25, 0.5],
                     [0, 0.5, 1]])
    C = [C_AA, C_Aa, C_aa]

    def get_probabilities(par1, par2):
        P = get_kernel(par1, par2)
        z = np.zeros(3)
        for k in range(3):
            for i in range(3):
                for j in range(3):
                    z[k] += C[k][i, j] * P[i, j]
        return z

    print(get_probabilities(par1, par2))


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

    def f(p: str) -> int:
        # Input: protein p as letter of english alphabet
        # Output: amount of different RNA string which encode that protein
        cnt = 0
        for k in proteins.keys():
            if proteins[k] == p:
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
        if is_there_common_substring(u=u, v=v, l=k):
            return f(u=u, v=v, l_1=k, l_2=l_2)
        else:
            return f(u=u, v=v, l_1=l_1, l_2=k)

    k = f(u=u, v=v, l_1=0, l_2=min(len(u), len(v)) + 1)
    common_substrings_list = common_substrings(u=u, v=v, l=k)
    common = []
    for t in common_substrings_list:
        common.append(u[t[0]:t[0] + k])
    return common


def someshitcodeinsteadofsuffixtree():
    seq_names, seqs = FASTA_to_lists()
    with open('output.txt', mode='w', encoding='utf-8') as g:
        seqs.sort()
        common = []
        for u in tqdm(seqs):
            for v in seqs:
                if u != v:
                    common.append(longest_common_substrings(u, v))

        # common = множество общих строк для всех
        def t_m(a):
            b = []
            for a_ in a:
                b.append(len(a_))
            return max(b)

        common_len = list(map(t_m, common))  # длины строк
        min_val = min(common_len)  # кратчайшие
        A = []  # те списки строк, где t_m > min_val
        B = []  # те, что t_m == val
        for i in range(len(common)):
            if t_m(common[i]) > min_val:
                A.append(common[i])
            else:
                B.append(common[i])
        # print("A: {}\n".format(A), "B: {}\n".format(B))
        # print(common, common_len)
        s0 = set(B[0])
        # print(s0)
        for b in B:
            b0 = set(b)
            s0.intersection_update(b0)
        print(s0)

        # s0 -- множество потенциальных кандидатов
        def is_lies(a: str, A: list) -> bool:
            # A is list of strings
            for a_ in A:
                if len(KMP(text=a_, pattern=a)) > 0:
                    return True
            return False

        def is_lies_in_all(s: str, A: list) -> bool:
            # A is list of lists
            for a in A:
                if not is_lies(s, a):
                    return False
            return True

        flag = False
        for s_ in s0:
            print("Is {} lies in all {}? Anser is {}".format(s_, A, is_lies_in_all(s_, A)))
            if is_lies_in_all(s_, A) and not flag:
                g.write(str(s_))
                flag = True
        if not flag:
            g.write("nixua net")


def DNA_to_Protein(s: str) -> str:
    # Input: dna string s
    # Output: every distinct candidate protein string that can be translated from ORFs of s
    # TODO: FIX
    mRNA = DNA_to_RNA(s)

    def f(a: np.ndarray, x: int) -> int:
        # Input: array a of integers, integer x
        # Output: smallest number y of a, greater or equal to x

        candidats = a[a >= x]
        print("candidats for {} in {}: ".format(x, candidats))
        for c in candidats:
            if (c - x) % 3 == 0:
                return c
        return None

    print(mRNA, '\n')
    start_codon = np.array(KMP(text=mRNA, pattern='AUG'))
    stop_codon = np.array(sorted(KMP(text=mRNA, pattern='UAA') +
                                 KMP(text=mRNA, pattern='UGA') +
                                 KMP(text=mRNA, pattern='UAG')))
    stop_codon += 3
    print(start_codon)
    print(stop_codon)

    for start in start_codon:
        x = f(stop_codon, start)

        print("Part of mRNA: ", mRNA[start:x])
        print("Positions: ", start, x)
        print("Protein: ", RNA_to_Protein(mRNA[start:x]))
        print('\n')


def permutations_dna(n: int) -> list:
    p = permutations(range(1, n + 1))

    def factorial(n):
        m = 1
        for i in range(1, n + 1):
            m *= i
        return m

    print(factorial(n))
    for p_ in p:
        print(*p_)


def CalculateProteinMass(p: str, mass_mode='float') -> float:
    # Input: protein p
    # Output: its mass in Da
    mass = {'A': 71.03711,
            'C': 103.00919,
            'D': 115.02694,
            'E': 129.04259,
            'F': 147.06841,
            'G': 57.02146,
            'H': 137.05891,
            'I': 113.08406,
            'K': 128.09496,
            'L': 113.08406,
            'M': 131.04049,
            'N': 114.04293,
            'P': 97.05276,
            'Q': 128.05858,
            'R': 156.10111,
            'S': 87.03203,
            'T': 101.04768,
            'V': 99.06841,
            'W': 186.07931,
            'Y': 163.06333}
    m = 0
    if mass_mode == 'int':
        for p_ in p:
            m += int(mass[p_])
    elif mass_mode == 'float':
        for p_ in p:
            m += mass[p_]
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
            # print("{} -> ...{} {} {}...".format(s[obj[0]:obj[0]+k],
            #                                     s[obj[0]-30:obj[0]],
            #                                     s[obj[0]:obj[0]+k],
            #                                     s[obj[0]+k:obj[0]+k+30]))
    return polys


class Node:
    def __init__(self, value):
        # value: str
        self.value = value
        self.links = []

    def add_link(self, link):
        # link: Link
        self.links.append(link)

    def __str__(self) -> str:
        node = '(%s):\n' % self.value
        for link in self.links:
            node += '\t' + link.__str__() + '\n'
        return node

    def equals(self, node) -> bool:
        # node: Node
        flag = (self.value == node.value)
        if len(self.links) == len(node.links):
            for i in range(len(self.links)):
                flag = flag and (self.links[i] == node.links[i])
            return flag
        else:
            return False


class Link:
    def __init__(self, tail, mark: str, head):
        # tail, head : Node
        # mark: str
        self.tail = tail
        self.head = head
        self.mark = mark

    def __str__(self) -> str:
        return '(%s --%s--> %s)' % (self.tail.value, self.mark, self.head.value)

    def equals(self, link) -> bool:
        # link: Link
        return (self.tail == link.tail) and (self.mark == link.mark) and (self.head == link.head)


class Automata:
    def __init__(self, initial_node: Node, nodes: list, terminal_node: Node):
        self.initial_node = initial_node
        self.nodes = nodes
        self.terminal_node = terminal_node

    def get_next_node(self, current_node: Node, mark: str):
        for link in current_node.links:
            if link.mark == mark:
                return link.head
        return None

    def accepts(self, string: str) -> bool:
        node = self.initial_node
        for character in string:
            node = self.get_next_node(node, character)
        return self.terminal_node.equals(node)

    def __str__(self):
        automata = "Initial node: %s\nTerminal node: %s\n" % (self.initial_node.value, self.terminal_node.value)
        for node in self.nodes:
            automata += node.__str__()
        return automata


def main():
    with open('input.txt', mode='r', encoding='utf-8') as f, open('output.txt', mode='w', encoding='utf-8') as g:
        s0 = Node("q0")
        s1 = Node("q1")
        s2 = Node("q2")

        s0_0_s0 = Link(s0, '0', s0)
        s0_1_s1 = Link(s0, '1', s1)
        s1_0_s2 = Link(s1, '0', s2)
        s1_1_s0 = Link(s1, '1', s0)
        s2_0_s1 = Link(s2, '0', s1)
        s2_1_s2 = Link(s2, '1', s2)

        s0.add_link(s0_0_s0)
        s0.add_link(s0_1_s1)
        s1.add_link(s1_0_s2)
        s1.add_link(s1_1_s0)
        s2.add_link(s2_0_s1)
        s2.add_link(s2_1_s2)

        a = Automata(s0, [s0, s1, s2], s0)


if __name__ == '__main__':
    main()
