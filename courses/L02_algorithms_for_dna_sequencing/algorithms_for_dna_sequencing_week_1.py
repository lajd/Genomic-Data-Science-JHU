import os
import random
from typing import Optional
from collections import defaultdict, deque
from matplotlib import pyplot as plt

from data import DATA_DIR


class BaseQualities:
    @staticmethod
    def QToPhred33(Q: int):
        return chr(Q + 33)

    @staticmethod
    def phred33ToQ(qual: int):
        return ord(qual) + 33


def random_embedded_genome(p: str, l: int = 100):
    rand_pos = random.randint(0, l - len(p))
    seq = ''
    for i in range(l):
        if i == rand_pos:
            seq += p
        else:
            seq += random.choice('ATCG')
    return seq, rand_pos


class NaiveAlignment:
    @staticmethod
    def _naive_alignment(p: str, t: str):
        occurences = []
        for i in range(len(t) - len(p) + 1):
            match = True
            for j in range(len(p)):
                if t[i + j] == p[j]:
                    # Same char
                    continue
                else:
                    match = False
                    break
            if match:
                occurences.append(i)
        return occurences

    @staticmethod
    def _naive_alignment_with_mismatches(p: str, t: str, n: int, ignore_alignments: dict=None):
        occurences = []
        matched_alignments = defaultdict(set)
        for i in range(len(t) - len(p) + 1):
            match = True
            num_mismatches = 0
            alignment = ''
            for j in range(len(p)):
                if t[i + j] == p[j]:
                    alignment += p[j]
                    # Same char
                    continue
                else:
                    alignment += p[j]
                    num_mismatches += 1
                    if num_mismatches > n:
                        match = False
                        break
            if match:
                if ignore_alignments and alignment in ignore_alignments and i in ignore_alignments[alignment]:
                    # Skip this alignment
                    continue
                else:
                    occurences.append(i)
                    matched_alignments[alignment].add(i)
        return occurences, matched_alignments

    @staticmethod
    def naive_alignments(p: str, t: str, include_reverse_compliment: bool = True,
                         allow_mismatches: bool = False, num_mismatches: int = None):

        if allow_mismatches is False:
            all_occurrences = []
            reverse_p = reverseComplement(p)
            occurrences = NaiveAlignment._naive_alignment(p, t)
            all_occurrences.extend(occurrences)

            if include_reverse_compliment and reverse_p != p:
                all_occurrences.extend(NaiveAlignment._naive_alignment(reverse_p, t))

            return all_occurrences

        else:
            all_occurrences = []
            all_alignments = defaultdict(set)
            reverse_p = reverseComplement(p)
            occurrences, alignments = NaiveAlignment._naive_alignment_with_mismatches(p, t, ignore_alignments=all_alignments, n=num_mismatches)
            all_occurrences.extend(occurrences)

            for k, v in alignments.items():
                if k in all_alignments:
                    all_alignments[k].union(v)
                else:
                    all_alignments.setdefault(k, v)

            if include_reverse_compliment and reverse_p != p:
                occurrences, alignments = NaiveAlignment._naive_alignment_with_mismatches(reverse_p, t, ignore_alignments=all_alignments, n=num_mismatches)
                all_occurrences.extend(occurrences)
                for k, v in alignments.items():
                    if k in all_alignments:
                        all_alignments[k].union(v)
                    else:
                        all_alignments.setdefault(k, v)
            return all_occurrences, all_alignments


def reverseComplement(s: str):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


def readGenome(filename: str):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


def readFastq(filename: str):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


def gcContentByPosition(reads: list, n: int):
    gc_content = [0] * n
    total_content = [0] * n
    for read in reads:
        for i, base in enumerate(read):
            if base in {'C', 'G'}:
                gc_content[i] += 1
            total_content[i] += 1

    gc_ratio = []
    for i in range(len(total_content)):
        if total_content[i] != 0:
            gc_ratio.append(gc_content[i] / total_content[i] )
        else:
            gc_ratio.append(0)
    return gc_ratio


if __name__ == '__main__':

    genome = readGenome(os.path.join(DATA_DIR, 'lambda_virus.fa'))

    # Problem 1
    occurrences = NaiveAlignment.naive_alignments('AGGT', genome)
    print("Problem 1 occurrence counts: %s", len(occurrences))

    # Problem 2
    occurrences = NaiveAlignment.naive_alignments('TTAA', genome)
    print("Problem 2 occurrence counts: %s", len(occurrences))

    # Problem 3
    occurrences = NaiveAlignment.naive_alignments('ACTAAGT', genome)
    print("Problem 3 occurrence counts: %s", len(occurrences))

    # Problem 4
    occurrences = NaiveAlignment.naive_alignments('AGTCGA', genome)
    print("Problem 4 occurrence counts: %s", len(occurrences))

    # Problem 5
    occurrences, alignments = NaiveAlignment.naive_alignments('TTCAAGCC', genome, allow_mismatches=True, num_mismatches=2, include_reverse_compliment=False)
    print("Problem 5 num occurrence counts: %s", len(occurrences))

    # Problem 6
    sequence_reads, qualities = readFastq(os.path.join(DATA_DIR, 'ERR037900_1.first1000.fastq'))
    gc_ratio = gcContentByPosition(sequence_reads, 100)
    plt.plot(range(len(gc_ratio)), gc_ratio)
    plt.show()

    print("Problem 6 offset: %s", 66)  # Identified from graph
