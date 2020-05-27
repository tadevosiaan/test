from collections import Counter
import constants
import itertools
import string


def NucleotidesCount(s: str) -> tuple:
    # Input: a DNA of letters A, T, C, G
    # Output: counts for each letter in the form (AGCT)
    # at most O(n)?
    return s.count('A'), s.count('G'), s.count('C'), s.count('T')


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
    return Pattern[::-1].translate(string.maketrans('ACGT', 'TGCA'))


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

def Hamming(s1, s2: str)->int:
    # Input: DNA strings s1,s2 of same len
    # Output: hamming distance s1,s2
    return sum(a!=b for a,b in itertools.izip(s1,s2))

def GC_content_amoungDNA():
    # Input: FASTA-encoded DNAs (>Rosalind_xxxx, DNA)
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

def main():
    with open('input.txt', mode='r', encoding='utf-8') as f:
        f1 = f.readline().rstrip()
        f2 = f.readline().rstrip()
        print(f1, f2)
        print(Hamming(f1,f2))
main()