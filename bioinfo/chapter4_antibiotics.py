import matplotlib.pyplot as plt
from collections import Counter
import time
from bio_constants import Proteins
from bio_utils import ReverseComplement, CalculateProteinMass, DNA_to_RNA, KMP


# CHAPTER 4. How Do We Sequence Antibiotics?
# Brute force algorithms


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
    # Complexity:  O(n)
    #
    # Def: We say that a DNA string Pattern encodes an amino acid string Peptide if the RNA
    # string transcribed from either Pattern or its reverse complement translates into
    # Peptide
    substrings = []
    substrings_rc = []
    text_rc = ReverseComplement(text)

    print(text)  # 5' ---> 3'
    print(text_rc)  # 3' ---> 5'

    for shift in range(0, 3):
        positions_in_dna = []
        dna = text[shift:]
        rna = DNA_to_RNA(dna)
        protein = RNA_to_Protein(rna)
        positions_in_protein = KMP(text=protein, pattern=peptide)
        if len(positions_in_protein) > 0:
            positions_in_dna = [3 * p + shift for p in positions_in_protein]
        for pos in positions_in_dna:
            substrings.append((pos, text[pos:pos + 3 * len(peptide)]))

        positions_in_dna_rc = []
        dna_rc = text_rc[shift:]
        rna_rc = DNA_to_RNA(dna_rc)
        protein_rc = RNA_to_Protein(rna_rc)
        positions_in_protein_rc = KMP(text=protein_rc, pattern=peptide)
        if len(positions_in_protein_rc) > 0:
            positions_in_dna_rc = [3 * p + shift for p in positions_in_protein_rc]
        for pos in positions_in_dna_rc:
            substrings_rc.append((pos, text_rc[pos:pos + 3 * len(peptide)]))
    return substrings + ['#'] + substrings_rc


# 4C Code challenge:
def GenerateTheoreticalCycloSpectrum(peptide: str, mass_mode='float') -> tuple:
    # Input: amino acid peptide (mass_mode = int or float, values of mass in Da)
    # Output: theoretical spectrum of peptide (CycloSpectrum)
    # Complexity: O(n^2), there is solution by prefix mass precount
    modified_peptide = peptide + peptide
    subpeptides = []
    mass = lambda p: CalculateProteinMass(p, mass_mode)
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
    # Complexity: O(n^2), there is solution by prefix mass precount
    subpeptides = []
    mass = lambda p: CalculateProteinMass(p, mass_mode)
    for l in range(1, len(peptide)):  # 1, 2, ..., n-1 ~~ O(|peptide| - 1)
        for shift in range(0, len(peptide) - l + 1):  # 0, 1, 2, ..., n- l ~~ O(|peptide| - l)
            subpeptides.append(peptide[shift: shift + l])
    keys = ['_'] + subpeptides + [peptide]
    values = [0] + [mass(v) for v in subpeptides] + [mass(peptide)]
    z = list(zip(keys, values))
    z.sort(key=lambda x: x[1])
    return z, sorted(values)


def plot_spectrum(spectrum: list, compare_spectrum=None):
    # Input: spectrum
    # Output: print and plot spectrum
    plt.figure()
    plt.grid(linestyle=':', linewidth=1)
    plt.xlabel('n')
    plt.ylabel('Da')
    plt.plot(spectrum, color='blue')
    if compare_spectrum is not None:
        plt.plot(compare_spectrum, color='red')
    plt.show()


# 4D Code challenge:
def CountingPeptideswithGivenMass(m: int) -> int:
    # Input: An integer m.
    # Output: The number of linear peptides having integer mass m.
    #
    # about ~~ k C^m for some constants C,k
    pass


# 4E Code challenge
def CyclopeptideSequencing(spectrum: list) -> str:
    # Input: theoretical spectrum = list of integers
    # Output: reconstruct a cyclic peptide, such that CycloSpectrum(peptide) = spectrum(return None if no peptide)
    # Complexity: O(n^a), but in practice faster
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
    mass = lambda p: CalculateProteinMass(p, mass_mode='int')

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
    print('None')
    return None


# 4F Code challenge
def PeptideScoring(peptide_type: str, peptide: str, spectrum: list) -> float:
    # Input: peptide type: 'linear' or 'cyclic' amino acid cyclic or linear peptide, list of (int) experimental spectrum
    # Output: compute the SCORE(Peptide, Spectrum)
    peptide_sp = []
    # compute theoretical spectrum
    if peptide_type == 'linear':
        peptide_sp = GenerateTheoreticalLinearSpectrum(peptide, mass_mode='int')[1]
    elif peptide_type == 'cyclic':
        peptide_sp = GenerateTheoreticalCycloSpectrum(peptide, mass_mode='int')[1]
    # calculate score
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
def LeaderBoardCyclopeptideSequencing(spectrum: list, N: int, a=None) -> str:
    # Input: A collection of integers (experimental) Spectrum, int N(for trim), alphabet a (optional)
    # Output: A cyclic peptide Peptide maximizing SCORE(Peptide, Spectrum) over
    # all peptides Peptide with mass equal to PARENTMASS(Spectrum)
    Leaderboard = ['']
    LeaderPeptide = ''
    mass = lambda p: CalculateProteinMass(p, mass_mode='int')

    def expand(lb: list, a=None) -> list:
        lb_ = []
        alphabet = []
        if a is None:
            alphabet = list(Mass_int.keys())
        else:
            # assuming a is alphabet of spectral convolution <=> list of tuples (acidoacidnum, cnt)
            for elem in a:
                if Mass_int['K'] == elem[0] or Mass_int['Q'] == elem[0]:
                    alphabet.append('K')
                    # alphabet.append('Q') #they are the same
                elif Mass_int['I'] == elem[0] or Mass_int['L'] == elem[0]:
                    # alphabet.append('I') #they are the same
                    alphabet.append('L')
                else:
                    for k in Mass_int.keys():
                        if Mass_int[k] == elem[0]:
                            alphabet.append(k)
            alphabet = list(set(alphabet))

        for p in lb:
            for a_ in alphabet:
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
        Leaderboard = expand(Leaderboard, a)
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
        Leaderboard = trim(Leaderboard, spectrum, N)
        print('iter: {}, LB: {} .trim. {}, score: {}'.format(iteration, l, len(Leaderboard),
                                                             CyclopeptideScoring(LeaderPeptide, spectrum)))
    print(LeaderPeptide)
    return LeaderPeptide


# 4H Code challenge
def SpectralConvolution(spectrum: list) -> Counter:
    # Input: list of numbers spectrum (assuming experimental spectrum)
    # Output: The list of elements in the convolution of Spectrum in decreasing order of their multiplicities
    # I will cut the array such that values will lie in range [57, 200]
    n = len(spectrum)
    cnt = Counter()
    for i in range(n):
        for j in range(i):
            t = spectrum[i] - spectrum[j]
            if 201 > t > 56:
                cnt[t] += 1
    return cnt


def get_SpectralAlphabet(spectrum: list, M: int) -> list:
    # Input: list of nums spectrum, int M
    # Output: alphabet of M (or more with'tie') amino acids such that they frequently appear in spectrum convolution
    spectrum = SpectralConvolution(spectrum)
    spc = sorted(list(spectrum.items()), key=lambda x: -x[1])
    M -= 1
    for j in range(M + 1, len(spc)):
        if spc[M][1] > spc[j][1]:
            spc = spc[:j]
            return spc
    return spc


# 4I Code challenge:
def ConvolutionCyclopeptideSequencing(spectrum: list, N: int, M: int) -> str:
    # Case of LeaderBoardCyclopeptideSequencing, where alphabet computed by spectral convolution
    alphabet = get_SpectralAlphabet(spectrum, M)
    return LeaderBoardCyclopeptideSequencing(spectrum, N, alphabet)
