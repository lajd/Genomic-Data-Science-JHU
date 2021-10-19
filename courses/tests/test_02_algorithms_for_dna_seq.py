import os
from unittest import TestCase

import pytest
from Bio.SeqRecord import SeqRecord
from parameterized import parameterized

from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_1 import NaiveAlignment, random_embedded_genome, readGenome
from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_2 import BoyerMoore
from courses.L02_algorithms_for_dna_sequencing.utils.boyer_moore_preproc import BoyerMoore as BoyerMoorePreprocessing

from data import DATA_DIR


class TestAlgorithmsForDNASequencingWeek1(TestCase):

    def test_naive_string_matching(self):
        # test basic
        s = 'TTTTTTTTTTAAAAAAAAAA'
        test_reads = [random_embedded_genome(s, 100)[0] for _ in range(100)]
        counts = 0
        for read in test_reads:
            occurrences, naive_alignments = NaiveAlignment.naive_alignments(s, read)
            counts += len(occurrences)

        self.assertEqual(counts, 100)

    def test_reverse_compliment(self):
        single_stranded_counts = 0
        reverse_compliment_counts = 0
        s = 'TTAA'
        test_reads = [random_embedded_genome(s, 100)[0] for _ in range(100)]
        for read in test_reads:
            single_stranded_counts += len(NaiveAlignment.naive_alignments(s, read, include_reverse_compliment=False)[0])
            reverse_compliment_counts += len(NaiveAlignment.naive_alignments(s, read, include_reverse_compliment=True)[0])

        self.assertEqual(reverse_compliment_counts, single_stranded_counts)

    def test_naive_alignment_with_mismatches(self):
        occurrences, naive_alignments = NaiveAlignment.naive_alignments('ACTTTA', 'ACTTACTTGATAAAGT', allow_mismatches=True,
                                                                  num_mismatches=2, include_reverse_compliment=False)
        self.assertEqual(occurrences, [0, 4])

    def test_naive_alignment_with_mismatches_2(self):
        p = 'CTGT'
        ten_as = 'AAAAAAAAAA'
        t = ten_as + 'CTGT' + ten_as + 'CTTT' + ten_as + 'CGGG' + ten_as
        occurrences, naive_alignments = NaiveAlignment.naive_alignments(p, t, allow_mismatches=True, num_mismatches=2,
                                                                  include_reverse_compliment=False)
        self.assertEqual(occurrences, [10, 24, 38])

    def test_naive_alignment_with_mismatches_3(self):
        occurrences, naive_alignments = NaiveAlignment.naive_alignments('GATTACA', readGenome(os.path.join(DATA_DIR, 'phix.fa')),
                                                                  allow_mismatches=True, num_mismatches=2,
                                                                  include_reverse_compliment=False)
        self.assertEqual(len(occurrences), 79)
        self.assertEqual(min(occurrences), 10)

    def test_naive_alignments_exact_1(self):

        p = 'word'
        t = 'there would have been a time for such a word'
        occurrences, num_comparisons = NaiveAlignment.naive_alignments(p, t, include_reverse_compliment=False)
        self.assertEqual(occurrences, [40])
        self.assertEqual(num_comparisons, 46)
        num_alignments = len(t) - len(p) + 1
        self.assertEqual(num_alignments, 41)

    def test_naive_alignments_exact_2(self):

        p = 'needle'
        t = 'needle need noodle needle'
        occurrences, num_comparisons = NaiveAlignment.naive_alignments(p, t, include_reverse_compliment=False)
        self.assertEqual(occurrences, [0, 19])
        self.assertEqual(num_comparisons, 35)
        num_alignments = len(t) - len(p) + 1
        self.assertEqual(num_alignments, 20)

    @staticmethod
    def get_boyer_moore(p: str, alphabet: str = 'ACGT'):
        p_bm = BoyerMoorePreprocessing(p, alphabet=alphabet)
        bm = BoyerMoore(p, p_bm=p_bm)
        return bm

    def test_boyer_moore(self):
        t = 'GTTATAGCTGATCGCGGCGTAGCGGCGAA'
        bm = self.get_boyer_moore('GTAGCGGCG')
        occurrences, alignment_tried, num_char_comparisons = bm.boyer_moore_alignment(t)
        self.assertEqual(occurrences, [18])
        self.assertEqual(
            alignment_tried,
            4
        )

    def test_boyer_moore2(self):
        t = 'CCGGTGTTTGAC'
        bm = self.get_boyer_moore('GATTATT')
        occurrences, alignment_tried, num_char_comparisons = bm.boyer_moore_alignment(t)
        self.assertEqual(
            alignment_tried,
            2
        )

    def test_boyer_moore3(self):
        p = 'word'
        t = 'there would have been a time for such a word'
        lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
        bm = self.get_boyer_moore(p, alphabet=lowercase_alphabet)
        occurrences, alignment_tried, num_char_comparisons = bm.boyer_moore_alignment(t)
        self.assertListEqual(
            occurrences,
            [40]
        )

        self.assertEqual(
            alignment_tried,
            12
        )

        self.assertEqual(num_char_comparisons, 15)

    def test_boyer_moore4(self):
        lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
        p = 'needle'
        t = 'needle need noodle needle'
        bm = self.get_boyer_moore(p, alphabet=lowercase_alphabet)
        occurrences, alignment_tried, num_character_comparisons = bm.boyer_moore_alignment(t)

        self.assertListEqual(
            occurrences,
            [0, 19]
        )

        self.assertEqual(
            alignment_tried,
            5
        )

        self.assertEqual(num_character_comparisons, 18)
