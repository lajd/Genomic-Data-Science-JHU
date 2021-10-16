import random
from collections import defaultdict

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


class StringMatching:
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
            occurrences = StringMatching._naive_alignment(p, t)
            all_occurrences.extend(occurrences)

            if include_reverse_compliment and reverse_p != p:
                all_occurrences.extend(StringMatching._naive_alignment(reverse_p, t))

            return all_occurrences

        else:
            all_occurrences = []
            all_alignments = defaultdict(set)
            reverse_p = reverseComplement(p)
            occurrences, alignments = StringMatching._naive_alignment_with_mismatches(p, t, ignore_alignments=all_alignments, n=num_mismatches)
            all_occurrences.extend(occurrences)

            for k, v in alignments.items():
                if k in all_alignments:
                    all_alignments[k].union(v)
                else:
                    all_alignments.setdefault(k, v)

            if include_reverse_compliment and reverse_p != p:
                occurrences, alignments = StringMatching._naive_alignment_with_mismatches(reverse_p, t, ignore_alignments=all_alignments, n=num_mismatches)
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


genome = readGenome('/Users/jonathan/PycharmProjects/genomics_data_science/data/lambda_virus.fa')

# test basic
s = 'ATCGCGCGATGCGCAT'
test_reads = [random_embedded_genome(s, 100)[0] for _ in range(100)]
counts = 0
for read in test_reads:
    occurrences = StringMatching.naive_alignments(s, read)
    counts += len(occurrences)

assert counts == 100

# test same reverse compliment
single_stranded_counts = 0
reverse_compliment_counts = 0
s = 'TTAA'
test_reads = [random_embedded_genome(s, 100)[0] for _ in range(100)]
for read in test_reads:
    single_stranded_counts += len(StringMatching.naive_alignments(s, read, include_reverse_compliment=False))
    reverse_compliment_counts += len(StringMatching.naive_alignments(s, read, include_reverse_compliment=True))

assert reverse_compliment_counts == single_stranded_counts


# Problem 1
counts = 0
occurrences = StringMatching.naive_alignments('AGGT', genome)
counts += len(occurrences)

# Problem 2
counts = 0
occurrences = StringMatching.naive_alignments('TTAA', genome)
counts += len(occurrences)

# Problem 3
occurrences = StringMatching.naive_alignments('ACTAAGT', genome)

# Problem 4
occurrences = StringMatching.naive_alignments('AGTCGA', genome)

# Problem 5
## Test 1
occurrences, alignments = StringMatching.naive_alignments('ACTTTA', 'ACTTACTTGATAAAGT', allow_mismatches=True, num_mismatches=2, include_reverse_compliment=False)
assert occurrences == [0, 4]

# Test 2
p = 'CTGT'
ten_as = 'AAAAAAAAAA'
t = ten_as + 'CTGT' + ten_as + 'CTTT' + ten_as + 'CGGG' + ten_as
occurrences, alignments = StringMatching.naive_alignments(p, t, allow_mismatches=True, num_mismatches=2, include_reverse_compliment=False)
assert occurrences == [10, 24, 38]
# Test 3
occurrences, alignments = StringMatching.naive_alignments('GATTACA', readGenome('/Users/jonathan/PycharmProjects/genomics_data_science/data/phix.fa'), allow_mismatches=True, num_mismatches=2, include_reverse_compliment=False)
assert len(occurrences) == 79
assert min(occurrences) == 10

## Real
occurrences, alignments = StringMatching.naive_alignments('TTCAAGCC', genome, allow_mismatches=True, num_mismatches=2, include_reverse_compliment=False)
print(len(occurrences))


# Problem 6
sequence_reads, qualities = readFastq('/Users/jonathan/PycharmProjects/genomics_data_science/data/ERR037900_1.first1000.fastq')

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

gc_ratio = gcContentByPosition(sequence_reads, 100)

from matplotlib import pyplot as plt
plt.plot(range(len(gc_ratio)), gc_ratio)
plt.show()


# Offset 66


### Boyer-Moore basics


class BoyerMoore:

    @staticmethod
    def boyer_moore_alignment(p: str, t: str):
        rev_p = p[::-1]
        alignment_start_idx = 0

        while alignment_start_idx < (len(t) - len(p) + 1):
            # Check from right to left
            for j in range(alignment_start_idx + len(p))
                j = alignment_start_idx
                if alignment_start_idx +



            for i in range(len(p)):
                if p[i] == t[i]:
                    continue

            for i, base in enumerate(rev_p):




