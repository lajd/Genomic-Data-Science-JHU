from typing import Optional, List, Iterable

import itertools
from typing import Union

from courses.course_3_algorithms_for_dna_seq.algorithms_for_dna_sequencing_week_3 import Overlap, OverlapSuffixContainedInRead


class KMerIndex(OverlapSuffixContainedInRead):
    def __init__(self):
        super().__init__()

    def get_candidate_pairs(self, reads: List[str], k: int) -> Iterable:
        kmer_to_read_id_map, read_id_mapping = self._get_map_kmer_to_reads(reads, k)

        for read_id, read in enumerate(reads):
            suffix = read[-k:]
            # Find which other reads contain this suffix
            # Don't match the read against itself
            matching_read_ids = kmer_to_read_id_map[suffix] - {read_id}
            for matching_read_id in matching_read_ids:
                yield read_id_mapping[read_id], read_id_mapping[matching_read_id]


class ShortestCommonSuperstring(OverlapSuffixContainedInRead):
    def __init__(self):
        self.overlap = Overlap()
        super().__init__()

    def brute_force_scs(self, reads: list, k: int) -> list:
        """ Returns shortest common superstring of given
            strings, which must be the same length """
        shortest_sup_length = None
        shortest_superstrings = []

        for ssperm in itertools.permutations(reads):
            sup = ssperm[0]  # superstring starts as first string
            for i in range(len(reads)-1):
                # overlap adjacent strings A and B in the permutation
                olen = self.overlap.overlap(ssperm[i], ssperm[i+1], min_length=k)
                # add non-overlapping portion of B to superstring
                sup += ssperm[i+1][olen:]
            if shortest_sup_length is None or len(sup) < shortest_sup_length:
                shortest_superstrings = [sup]
                # shortest_sup = sup  # found shorter superstring
                shortest_sup_length = len(sup)
            elif len(sup) == shortest_sup_length:
                shortest_superstrings += [sup]

        return shortest_superstrings  # return shortest

    def find_maximal_overlapping_reads(self, reads: Union[dict, list], k: int, kmer_index_k: Optional[int] = None) -> tuple:
        """ Find the pair of reads with the greatest overlap
        :param reads:
        :param k:
        :param kmer_index_k:
        :return:
        """
        a = b = None
        best_overlap_len = 0

        if kmer_index_k is not None:
            iterable = KMerIndex().get_candidate_pairs(reads, k=kmer_index_k)
        else:
            iterable = itertools.permutations(reads, 2)

        for a_, b_ in iterable:
            # Require k=best_overlap_len, since otherwise we will discard this pair
            # in favour of the best overlap
            overlap_len = self.overlap.overlap(a_, b_, k)
            if overlap_len > best_overlap_len:
                a = a_
                b = b_
                best_overlap_len = overlap_len
        return a, b, best_overlap_len

    def greedy_scs(self, reads: list, k: int, kmer_index_k: Optional[int] = None) -> str:
        a, b, overlap_len = self.find_maximal_overlapping_reads(reads, k, kmer_index_k=kmer_index_k)
        while overlap_len > 0:
            reads.remove(a)
            reads.remove(b)
            # C is the concatenation of a and b (with overlap removed)
            c = a + b[overlap_len:]
            reads.append(c)
            a, b, overlap_len = self.find_maximal_overlapping_reads(reads, k, kmer_index_k=kmer_index_k)
        return "".join(reads)
