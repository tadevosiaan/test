from functools import reduce
import bio
import matplotlib.pyplot as plt
from collections import Counter
import time

# CHAPTER 4. How Do We Sequence Antibiotics?
# ______________________________
# 1 The discovery of antibiotics
# ______________________________


# there is no protein Z avaliable, here 'Z' means 'STOP'~'STP' codon
AminoAcids = \
    {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K', 'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F',
     'Asn': 'N', 'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W', 'Ala': 'A', 'Val': 'V', 'Glu': 'E',
     'Tyr': 'Y', 'Met': 'M', 'STP': 'Z'}
Mass = {'A': 71.03711,
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
Mass_int = {'A': 71,
            'C': 103,
            'D': 115,
            'E': 129,
            'F': 147,
            'G': 57,
            'H': 137,
            'I': 113,
            'K': 128,
            'L': 113,
            'M': 131,
            'N': 114,
            'P': 97,
            'Q': 128,
            'R': 156,
            'S': 87,
            'T': 101,
            'V': 99,
            'W': 186,
            'Y': 163}
Proteins = {'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V',
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

tyrocidine = 'VKLFPWFNQY'

"""
1928, Alexander Fleming, Penicillium - Staphylococcus bacteria
"""

# ____________________________________
# 2 How do bacteria make antibiotics?
# ____________________________________

""" 
Tyrocidine B1 -- one of many antibiotics produced by Baciliius brevis
Three-letter notation: Val Lys Leu Phe Pro Trp Phe Asn Gln Tyr
One-letter notation:    V   K   L   F   P   W   F   N   Q   Y


Central dogma of molecular biology:
DNA --- transcription ---> RNA --- translation ---> Protein
        RNA-polymerase              ribosome, NRP(non-ribosomal peptide) synthetase
"""


# 4A Code Challenge
def RNA_to_Protein(s: str, start=False, stop=False) -> str:
    # Input: RNA string s in form AUG...STOPCODON
    # Output: transcribe of s into peptide
    # Complexity: O(n)

    peptide = ''
    if start:
        assert s[:3] == 'AUG'
    if stop:
        assert s[-3:] in ['UAA', 'UAG', 'UGA']

    # assuming there is no stop codon at the middle

    for i in range((len(s) // 3)):
        peptide += Proteins[s[3 * i: 3 * i + 3]]
    return peptide


# 4B Code Challenge
def PeptideEncoding(text: str, peptide: str) -> list:
    # Input: : A DNA string Text and an amino acid string Peptide.
    # Output: All substrings of Text encoding Peptide (if any such substrings exist),
    # actually pairs (pos, substr in text/text_rc starting at pos), '#' splits positions at text and at text_rc
    # Complexty:  O(n)
    #
    # Def: We say that a DNA string Pattern encodes an amino acid string Peptide if the RNA
    # string transcribed from either Pattern or its reverse complement translates into
    # Peptide

    substrings = []
    substrings_rc = []
    text_rc = bio.ReverseComplement(text)

    print(text)  # 5' ---> 3'
    print(text_rc)  # 3' ---> 5'

    for shift in range(0, 3):
        positions_in_dna = []
        dna = text[shift:]
        rna = bio.DNA_to_RNA(dna)
        protein = RNA_to_Protein(rna)
        positions_in_protein = bio.KMP(text=protein, pattern=peptide)
        if len(positions_in_protein) > 0:
            positions_in_dna = [3 * p + shift for p in positions_in_protein]
        for pos in positions_in_dna:
            substrings.append((pos, text[pos:pos + 3 * len(peptide)]))

        positions_in_dna_rc = []
        dna_rc = text_rc[shift:]
        rna_rc = bio.DNA_to_RNA(dna_rc)
        protein_rc = RNA_to_Protein(rna_rc)
        positions_in_protein_rc = bio.KMP(text=protein_rc, pattern=peptide)
        if len(positions_in_protein_rc) > 0:
            positions_in_dna_rc = [3 * p + shift for p in positions_in_protein_rc]
        for pos in positions_in_dna_rc:
            substrings_rc.append((pos, text_rc[pos:pos + 3 * len(peptide)]))
    return substrings + ['#'] + substrings_rc


# ____________________________________
# 3 Sequencing antibiotics by shattering them into pieces
# ____________________________________

"""
Glycine = C2H3ON -> 57 mass, but Glucine~amino acid G

EXERCISE BREAK: How many subpeptides does a cyclic peptide of length n
have? 
For example: NQEL has 12 subpeptides: N, Q, E, L, NQ,
QE, EL, LN, NQE, QEL, ELN, and LNQ. 
ANSWER: n(n-1)
"""


# 4C Code challenge:
def GenerateTheoreticalCycloSpectrum(peptide: str, mass_mode='float') -> tuple:
    # Input: amino acid peptide (mass_mode = int or float, values of mass in Da)
    # Output: theoretical spectrum of peptide (CycloSpectrum)
    # Complexty: O(n^2), there is solution by prefix mass precount

    modified_peptide = peptide + peptide
    subpeptides = []
    mass = lambda p: bio.CalculateProteinMass(p, mass_mode)
    for l in range(1, len(peptide)):  # 1, 2, ..., n-1 ~~ O(|peptide|-1)
        for shift in range(0, len(peptide)):  # 0, 1, 2, ..., n-1 ~~O(|peptide|)
            subpeptides.append(modified_peptide[shift: shift + l])
    keys = ['_'] + subpeptides + [peptide]
    values = [0] + [mass(v) for v in subpeptides] + [
        mass(peptide)]
    z = list(zip(keys, values))
    z.sort(key=lambda x: x[1])
    return z, sorted(values)


def GenerateTheoreticalLinearSpectrum(peptide: str, mass_mode='float') -> tuple:
    # Input: amino acid peptide (mass_mode = int or float, values of mass in Da)
    # Output: theoretical spectrum of peptide (LinearSpectrum)
    # Complexty: O(n^2), there is solution by prefix mass precount

    subpeptides = []
    mass = lambda p: bio.CalculateProteinMass(p, mass_mode)
    for l in range(1, len(peptide)):  # 1, 2, ..., n-1 ~~ O(|peptide| - 1)
        for shift in range(0, len(peptide) - l + 1):  # 0, 1, 2, ..., n- l ~~ O(|peptide| - l)
            subpeptides.append(peptide[shift: shift + l])
    keys = ['_'] + subpeptides + [peptide]
    values = [0] + [mass(v) for v in subpeptides] + [mass(peptide)]
    z = list(zip(keys, values))
    z.sort(key=lambda x: x[1])
    return z, sorted(values)


def plot_spectrum(spectrum: list, compare_spectrum=None):
    # Input: theoretical spectrum
    # Output: print and plot spectrum

    plt.figure()
    plt.grid(linestyle=':', linewidth=1)
    plt.xlabel('n')
    plt.ylabel('Da')
    plt.plot(spectrum, color='blue')
    if compare_spectrum is not None:
        plt.plot(compare_spectrum, color='red')
    plt.show()


# ____________________________________
# 4 a brute force algorithm for cyclopeptide sequencing
# ____________________________________


# 4D Code challenge:
def CountingPeptideswithGivenMass(m: int) -> int:
    # Input: An integer m.
    # Output: The number of linear peptides having integer mass m.
    #
    # about ~~ k C^m for some constants C,k
    pass


# ____________________________________
# 5 cyclopeptide sequencing with branch-and-bound
# ____________________________________

"""
BRANCH-AND-BOUND algorithms: each such algorithm consists of a branching
step to increase the number of candidate solutions, followed by a bounding step to
remove hopeless candidates


EXERCISE BREAK: How many subpeptides does a linear peptide of length n
have?
ANSWER: 2+3+...+n = n(n-1)/2 - 1

BOUNDING STEP:
Given an experimental spectrum Spectrum of a cyclic peptide,
a linear peptide is consistent with Spectrum 
if every mass in its theoretical spectrum is contained in Spectrum.
If a mass appears more than once in the theoretical spectrum of the linear peptide, then
it must appear at least that many times in Spectrum in order for the linear peptide to be consistent with Spectrum.
"""


# 4E Code challenge
def CyclopeptideSequencing(spectrum: list) -> str:
    # Input: theoretical spectrum = list of integers
    # Output: reconstruct a cyclic peptide, such that CycloSpectrum(peptide) = spectrum(return None if no peptide)
    # Complexty: O(n^a), but in practice faster
    # Answer is accurate to a cyclic shift and/or reverse

    peptides = ['']

    def get_alphabet(spectrum: list) -> list:
        alphabet = []
        for k in Mass_int.keys():
            if Mass_int[k] in spectrum:
                alphabet.append(k)
        return alphabet

    def list_equivalence(l1: list, l2: list) -> bool:
        a = list((Counter(l1) - Counter(l2)).elements())
        b = list((Counter(l2) - Counter(l1)).elements())
        return len(a) == 0 and len(b) == 0

    def is_consistent(p: str) -> bool:
        sp_c = spectrum.copy()
        sp = GenerateTheoreticalLinearSpectrum(peptide=p, mass_mode='int')[1]
        list_diff = list((Counter(sp) - Counter(sp_c)).elements())
        return not len(list_diff) > 0

    a = get_alphabet(spectrum)
    mass = lambda p: bio.CalculateProteinMass(p, mass_mode='int')

    def expand(peptides: list) -> list:
        pep = []
        for p in peptides:
            for a_ in a:
                pep.append(p + a_)
        return pep

    iteration = 0
    while len(peptides) > 0:
        iteration += 1
        peptides = expand(peptides)
        print('iter: {}, peptides: {}'.format(iteration, len(peptides)))
        i = 0
        while i < len(peptides):
            p = peptides[i]
            i += 1
            if mass(p) == max(spectrum):
                if list_equivalence(GenerateTheoreticalCycloSpectrum(peptide=p, mass_mode='int')[1], spectrum):
                    print(p)
                    return p
                peptides.remove(p)
                i -= 1
            elif not is_consistent(p):
                peptides.remove(p)
                i -= 1
    return None


# ____________________________________
# 6 adapting sequencing for spectra with errors
# ____________________________________

# 4F Code challenge
def PeptideScoring(peptide_type: str, peptide: str, spectrum: list) -> float:
    # Input: peptide type: 'linear' or 'cyclic' amino acid cyclic or linear peptide, list of (int) spectrum
    # Output: compute the SCORE(Peptide, Spectrum)

    peptide_sp = []
    if peptide_type == 'linear':
        peptide_sp = GenerateTheoreticalLinearSpectrum(peptide, mass_mode='int')[1]
    elif peptide_type == 'cyclic':
        peptide_sp = GenerateTheoreticalCycloSpectrum(peptide, mass_mode='int')[1]
    score = list((Counter(peptide_sp) - Counter(spectrum)).elements())
    return len(peptide_sp) - len(score)


def CyclopeptideScoring(peptide: str, spectrum: list) -> float:
    return PeptideScoring(peptide_type='cyclic',
                          peptide=peptide,
                          spectrum=spectrum)


def LinearPeptideScoring(peptide: str, spectrum: list) -> float:
    return PeptideScoring(peptide_type='linear',
                          peptide=peptide,
                          spectrum=spectrum)


# 4G Code challenge
def LeaderBoardCyclopeptideSequencing(spectrum: int, N: int) -> str:
    # Input: A collection of integers Spectrum.
    # Output: A cyclic peptide Peptide maximizing SCORE(Peptide, Spectrum) over
    # all peptides Peptide with mass equal to PARENTMASS(Spectrum)

    Leaderboard = ['']
    LeaderPeptide = ''
    mass = lambda p: bio.CalculateProteinMass(p, mass_mode='int')

    def expand(lb: list) -> list:
        lb_ = []
        a = list(Mass_int.keys())
        for p in lb:
            for a_ in a:
                lb_.append(p + a_)
        return lb_

    def trim(lb: list, spectrum: list, N: int) -> list:
        linear_scores = []
        for j in range(len(lb)):
            linear_scores.append(LinearPeptideScoring(lb[j], spectrum))
        N -= 1
        z = list(zip(lb, linear_scores))
        z.sort(key=lambda x: -x[1])
        zkeys = []
        for j in range(N + 1, len(lb)):
            if z[N][1] > z[j][1]:
                z = z[:j]
                break
        for i in range(len(z)):
            zkeys.append(z[i][0])
        return zkeys

    iteration = 0
    while len(Leaderboard) > 0:
        iteration += 1
        Leaderboard = expand(Leaderboard)
        i = 0
        while i < len(Leaderboard):
            p = Leaderboard[i]  # peptide from Leaderboard
            i += 1
            if mass(p) == max(spectrum):
                if CyclopeptideScoring(p, spectrum) > CyclopeptideScoring(LeaderPeptide, spectrum):
                    LeaderPeptide = p
            elif mass(p) > max(spectrum):
                Leaderboard.remove(p)
                i -= 1
        l = len(Leaderboard)
        Leaderboard = trim(lb=Leaderboard, spectrum=spectrum, N=N)
        print('iter: {}, LB: {} .trim. {}, Score: {}'.format(iteration, l, len(Leaderboard),
                                                             CyclopeptideScoring(LeaderPeptide, spectrum)))
    print(LeaderPeptide)
    return LeaderPeptide


s = [0, 97, 99, 113, 114, 115, 128, 128, 147, 147, 163, 186, 227, 241, 242, 244, 244, 256, 260, 261, 262, 283, 291, 309,
     330, 333, 340, 347, 357, 385, 388, 389, 390, 390, 405, 430, 430, 435, 447, 485, 487, 503, 504, 518, 543, 544, 552,
     575, 577, 584, 599, 608, 631, 632, 650, 651, 653, 671, 672, 690, 691, 717, 738, 745, 747, 770, 778, 779, 804, 818,
     819, 827, 835, 837, 875, 892, 892, 917, 932, 932, 933, 934, 965, 982, 989, 1031, 1039, 1060, 1061, 1062, 1078,
     1080, 1081, 1095, 1136, 1159, 1175, 1175, 1194, 1194, 1208, 1209, 1223, 1225, 1322]
start = time.time()
g = LeaderBoardCyclopeptideSequencing(spectrum=s, N=50)
print('{} s'.format(time.time() - start))
# ____________________________________
# 7 from 20 to more than 100 amino acids
# ____________________________________

# ____________________________________
# 8 the spectral convolution saves the day
# ____________________________________

# ____________________________________
# 9 the truth about spectra
# ____________________________________
