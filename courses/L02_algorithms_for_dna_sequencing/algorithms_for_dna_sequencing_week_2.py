import os
from abc import ABC, abstractmethod
import random
from typing import Optional, List, Type
import bisect
from collections import defaultdict, deque
from matplotlib import pyplot as plt
from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_1 import readGenome

from data import DATA_DIR

from courses.L02_algorithms_for_dna_sequencing.utils.boyer_moore_preproc import BoyerMoorePreprocessing


class ExactMatchingStrategy(ABC):
    """
    The Strategy interface declares operations common to all supported versions
    of some algorithm.

    The Context uses this interface to call the algorithm defined by Concrete
    Strategies.
    """

    @abstractmethod
    def get_occurrences(self, **kwargs):
        pass

    @abstractmethod
    def prepare(self, **kwargs):
        pass



### Boyer-Moore basics
class BoyerMooreExact(ExactMatchingStrategy):
    def __init__(self, p: str, p_bm: Optional[BoyerMoorePreprocessing] = None, alphabet: str = 'ACGT', **kwargs):
        if not p_bm:
            p_bm = BoyerMoorePreprocessing(p, alphabet=alphabet)
        self.p = p
        self.p_bm = p_bm

    def prepare(self, **kwargs):
        pass

    def get_occurrences(self, t: str, **kwargs):
        return self.query(t)[0]

    def query(self, t: str):
        occurrences = []
        skipped_alignments = []
        num_char_comparisons = 0
        alignment_start_idx = 0
        num_alignments_tried = 0

        while alignment_start_idx <= (len(t) - len(self.p)):
            num_alignments_tried += 1
            shift = 1
            match = True
            # Check from right to left
            for j in range(len(self.p) - 1, -1, -1):
                num_char_comparisons += 1
                k = alignment_start_idx + j
                if self.p[j] != t[k]:
                    match = False
                    bad_char_shift = self.p_bm.bad_character_rule(j, t[k])
                    good_suffix_shift = self.p_bm.good_suffix_rule(j)
                    shift = max(bad_char_shift, good_suffix_shift, shift)
                    if shift > 1:
                        skipped_alignments.append(
                            {
                                'shift': shift - 1,
                                'bad_char_shift':  bad_char_shift - 1 if bad_char_shift > 0 else 0,
                                'good_suffix_shift': good_suffix_shift - 1 if good_suffix_shift > 0 else 0
                            }
                        )
                    break

            if match:
                occurrences.append(alignment_start_idx)
                skip_gs = self.p_bm.match_skip()
                shift = max(shift, skip_gs)
            alignment_start_idx += shift
        return occurrences, num_alignments_tried, num_char_comparisons


class KMERIndex(ExactMatchingStrategy):
    def __init__(self, t, k=8, **kwargs):
        ''' Create index from all substrings of size 'length' '''
        self.t = t
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i + k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    def prepare(self, **kwargs):
        pass

    def _get_hits(self, p: str):
        ''' Return index hits for first k-mer of P '''
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

    def get_occurrences(self, p: str, **kwargs):
        k = self.k
        offsets = []
        for i in self._get_hits(p):
            if p[k:] == self.t[i+k:i+len(p)]:
                offsets.append(i)
        return sorted(offsets)


class Context:
    def __init__(self, strategy: ExactMatchingStrategy) -> None:
        self._strategy = strategy

    @property
    def strategy(self) -> ExactMatchingStrategy:
        return self._strategy

    @strategy.setter
    def strategy(self, strategy: ExactMatchingStrategy) -> None:
        self._strategy = strategy

    def get_occurrences(self, **kwargs) -> List[int]:
        occurrences = self._strategy.get_occurrences(**kwargs)
        return occurrences

    def prepare(self, **kwargs) -> None:
        self._strategy.prepare(**kwargs)


class PigeonHoleApproximateMatching:

    @staticmethod
    def query(p: str, t: str, m: int, alphabet='ACGT', method='boyer_moore', **kwargs):
        partition_length = int(round(len(p) / (m + 1)))
        occurrences = set()

        matcher = None
        if method == 'kmer_index':
            matcher = KMERIndex(t=t, **kwargs)

        for i in range(m + 1):
            partition_start = i * partition_length
            partition_end = min(partition_start + partition_length, len(p))
            sub_p = p[partition_start:partition_end]
            if not sub_p:
                break

            # matcher = self.matching_strategy(p=sub_p, alphabet=alphabet)
            if method == 'boyer_moore':
                matcher = BoyerMooreExact(p=sub_p, **kwargs)

            occurrences_ = matcher.get_occurrences(p=sub_p, alphabet=alphabet, t=t)
            # bm = BoyerMooreExact(sub_p, alphabet=alphabet)
            # occurrences_, alignments_tried_, num_comparisons_ = bm.query(t)

            # For any exact matches found, perform a validation step. Look around the exact match,
            # and see if the text matches (allowing a certain number of mismatches). If the text
            # matches within a certain number of mismatches, we have found a match
            for match in occurrences_:
                # This match occurs outside of the range of this partition, once aligned with t
                if match < partition_start or (match - partition_start + len(p)) > len(t):
                    continue
                else:
                    mismatches = 0

                    # Test the part of p before the partition we already compared abovce
                    for j in range(0, partition_start):
                        if p[j] != t[match-partition_start + j]:
                            mismatches += 1
                            if mismatches > m:
                                break
                    # Compare section after the segment already tested
                    for j in range(partition_end, len(p)):
                        if p[j] != t[match-partition_start + j]:
                            mismatches += 1
                            if mismatches > m:
                                break

                    if mismatches <= m:
                        occurrences.add(match - partition_start)
        return sorted(occurrences)


class Index(object):
    def __init__(self, t, k):
        ''' Create index from all substrings of size 'length' '''
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i + k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    def query(self, p):
        ''' Return index hits for first k-mer of P '''
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
