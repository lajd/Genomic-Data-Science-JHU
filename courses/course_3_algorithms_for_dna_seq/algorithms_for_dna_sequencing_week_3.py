from collections import defaultdict
from typing import Optional, List, Iterable


import numpy as np

from courses.course_3_algorithms_for_dna_seq.algorithms_for_dna_sequencing_week_2 import SubseqIndex


class EditDistance:
    def __init__(self):
        self._cache = {}
        self.calls = 0

    def ed_recursive(self, a: str, b: str) -> int:
        """
        :param a: String 1
        :param b: String 2

        alpha -> a[0:-1]
        beta -> b[0:-1]

        :return:
        """
        if (a, b) in self._cache:
            return self._cache[(a, b)]
        else:
            if len(a) == 0:
                return len(b)
            elif len(b) == 0:
                return len(a)
            else:
                x = a[-1]
                y = b[-1]
                delta = 1 if x != y else 0
                alpha = a[0:-1]
                beta = b[0:-1]
                distance = min(
                    self.ed_recursive(alpha, beta) + delta,
                    self.ed_recursive(alpha + x, beta) + 1,
                    self.ed_recursive(alpha, beta + y) + 1,
                )
                self._cache[(a, b)] = distance
                self.calls += 1
                return distance

    def _ed_dp(self, a: str, b: str, d: Optional[np.ndarray] = None) -> np.ndarray:
        """

        :param a:
        :param b:
        :param d: Initial distance matrix
        :return:
        """
        if d is None:
            d = np.zeros((len(a) + 1, len(b) + 1))
            d[:, 0] = np.arange(0, len(a) + 1)
            d[0, :] = np.arange(0, len(b) + 1)

        for i in range(1, len(a) + 1):
            for j in range(1, len(b) + 1):
                self.calls += 1
                d_right = d[i][j-1] + 1
                d_down = d[i - 1][j] + 1
                delta = 0 if a[i-1] == b[j-1] else 1
                d_diag = d[i-1][j-1] + delta
                d[i][j] = min(d_right, d_down, d_diag)
        return d

    def ed_dp(self, a: str, b: str, d: Optional[np.ndarray] = None) -> np.ndarray:
        """

        :param a:
        :param b:
        :param d: Initial distance matrix
        :return:
        """
        return self._ed_dp(a, b, d)[-1][-1]


class ApproximateMatching(EditDistance):
    # TODO: Support indexing
    def __init__(self):
        super().__init__()

    def closest_match(self, a: str, b: str) -> float:
        # Top row is zeros
        d = np.zeros((len(a) + 1, len(b) + 1))
        d[:, 0] = np.arange(0, len(a) + 1)
        d = self._ed_dp(a, b, d=d)
        return min(d[-1])


class Overlap:
    def __init__(self):
        self._overlap_cache = {}
        self._num_cache_hits = 0

    def overlap(self, a, b, min_length=3) -> int:
        """ Return length of longest suffix of 'a' matching
            a prefix of 'b' that is at least 'min_length'
            characters long.  If no such overlap exists,
            return 0. """

        # If the strings aren't the same length, shorten
        # the strings such that they are. This will
        # help with caching
        max_compared_length = min(len(a), len(b))
        a = a[-max_compared_length:]
        b = b[0:max_compared_length]
        if (a, b) in self._overlap_cache:
            self._num_cache_hits += 1
            return self._overlap_cache[(a, b)]
        start = 0  # start all the way at the left
        while True:
            start = a.find(b[:min_length], start)  # look for b's prefix in a
            if start == -1:  # no more occurrences to right
                overlap = 0
                self._overlap_cache[(a, b)] = overlap
                return overlap
            # found occurrence; check for full suffix/prefix match
            if b.startswith(a[start:]):
                overlap = len(a)-start
                self._overlap_cache[(a, b)] = overlap
                return overlap
            start += 1  # move just past previous match


class MatchExactPrefixSuffix(Overlap):
    """
    Obtain matches where the suffix of one read exactly matches
     the prefix of pthers
    """
    def __init__(self):
        super().__init__()

    @staticmethod
    def _get_exact_prefix_suffix_index(reads: List[str], k: int):
        # Extract length k prefix/suffix from each read
        prefix_mapping = defaultdict(set)
        suffix_mapping = defaultdict(set)

        read_id_mapping = {}

        for read_id, read in enumerate(reads):
            read_prefix = read[0:k]
            read_suffix = read[-k:]
            prefix_mapping[read_prefix].add(read_id)
            suffix_mapping[read_suffix].add(read_id)
            read_id_mapping[read_id] = read

        return prefix_mapping, suffix_mapping, read_id_mapping

    @staticmethod
    def _convert_adjacency_to_pairs(adj_dict: dict, node_id_mapping: Optional[dict] = None):
        pairs = []
        for node_from in adj_dict:
            for node_to in adj_dict[node_from]:
                if node_id_mapping:
                    pairs.append((node_id_mapping[node_from], node_id_mapping[node_to]))
                else:
                    pairs.append((node_from, node_to))
        return pairs

    def get_overlapping_pairs(self, reads: List[str], k: int, return_pairs: bool = False):
        """ Get the pairs which overlap exactly with a prefix/suffix length of k
        :param reads:
        :param k:
        :param return_pairs:
        :return:
        """
        suffix_prefix_read_pairs = defaultdict(set)
        prefix_mapping, suffix_mapping, read_id_mapping = self._get_exact_prefix_suffix_index(reads, k)

        # For each suffix, identify the matching prefixes
        for suffix, reads_from_set in suffix_mapping.items():
            for read_from in reads_from_set:
                suffix_prefix_read_pairs[read_from] = prefix_mapping[suffix]

        if return_pairs:
            pairs = self._convert_adjacency_to_pairs(suffix_prefix_read_pairs, node_id_mapping=read_id_mapping)
            return pairs
        else:
            return suffix_prefix_read_pairs


class OverlapSuffixContainedInRead(MatchExactPrefixSuffix):
    """
    Obtain matches where the suffix of one read occurs within
    another read (does not require exact suffix/prefix matching)
    """
    def __init__(self):
        super().__init__()
        self._kmer_cache = defaultdict(list)

    def _get_kmers(self, read: str, k: int) -> list:

        if (read, k) in self._kmer_cache:
            return self._kmer_cache[(read, k)]
        else:
            kmers = []
            i = 0
            while i + k < len(read) + 1:
                kmers.append(read[i: i + k])
                i += 1
            self._kmer_cache[(read, k)] = kmers
            return kmers

    def _get_map_kmer_to_reads(self, reads: List[str], k: int):
        kmer_to_read_id_map = defaultdict(set)
        read_id_mapping = {}

        for read_id, read in enumerate(reads):
            for kmer in self._get_kmers(read, k):
                kmer_to_read_id_map[kmer].add(read_id)
            read_id_mapping[read_id] = read

        return kmer_to_read_id_map, read_id_mapping

    def get_overlapping_pairs(self, reads: List[str], k: int, return_pairs: bool = False):
        """ Get the pairs which overlap such that the suffix of read a
        is contained (as a kmer) in read b
        :param reads:
        :param k:
        :param return_pairs:
        :return:
        """
        kmer_to_read_id_map, read_id_mapping = self._get_map_kmer_to_reads(reads, k)

        adj = defaultdict(dict)

        for read_id, read in enumerate(reads):
            suffix = read[-k:]
            # Find which other reads contain this suffix
            # Don't match the read against itself
            matching_read_ids = kmer_to_read_id_map[suffix] - {read_id}
            for matching_read_id in matching_read_ids:
                matching_read = read_id_mapping[matching_read_id]
                overlap = self.overlap(read, matching_read, min_length=k)
                if overlap >= k:
                    adj[read_id][matching_read_id] = overlap

        if return_pairs:
            return self._convert_adjacency_to_pairs(adj, read_id_mapping)
        else:
            return adj
