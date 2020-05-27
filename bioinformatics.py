import collections
import constants


def NucleotidesCount(s: str) -> tuple:
    # Input: a DNA of letters A, T, C, G
    # Output: counts for each letter in the form (AGCT)
    # at most O(n)?
    return s.count('A'), s.count('G'), s.count('C'), s.count('T')


def DNAintoRNA(s: str) -> str:
    # Input: string of DNA
    # Output: transcribe of DNA into RNA
    # at most O(n)
    return s.replace('T', 'U')


def PatternCount(Text: str, Pattern: str) -> int:
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
    from string import maketrans
    return Pattern[::-1].translate(maketrans('ACGT', 'TGCA'))


def PatternMatching(Pattern: str, Genome: str) -> list:
    # Input: string Pattern and Genome
    # Output: all starting positions in Genome where Patterns appears as a substring
    pass
