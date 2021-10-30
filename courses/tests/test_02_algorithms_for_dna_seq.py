import os
from itertools import product
from unittest import TestCase

from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_1 import NaiveAlignment, random_embedded_genome, readGenome
from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_2 import PigeonHoleApproximateMatching, BoyerMooreExact
from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_3 import EditDistance, ApproximateMatching, Overlap, MatchExactPrefixSuffix, OverlapSuffixContainedInRead
from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_4 import ShortestCommonSuperstring

from courses.L02_algorithms_for_dna_sequencing.utils.boyer_moore_preproc import BoyerMoorePreprocessing

from data import DATA_DIR


class TestNaiveStringMatching(TestCase):

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


class TestExactBoyerMooreStringMatching(TestCase):
    @staticmethod
    def get_boyer_moore(p: str, alphabet: str = 'ACGT'):
        p_bm = BoyerMoorePreprocessing(p, alphabet=alphabet)
        bm = BoyerMooreExact(p, p_bm=p_bm)
        return bm

    def test_boyer_moore(self):
        t = 'GTTATAGCTGATCGCGGCGTAGCGGCGAA'
        bm = self.get_boyer_moore('GTAGCGGCG')
        occurrences, alignment_tried, num_char_comparisons = bm.query(t)
        self.assertEqual(occurrences, [18])
        self.assertEqual(
            alignment_tried,
            4
        )

    def test_boyer_moore2(self):
        t = 'CCGGTGTTTGAC'
        bm = self.get_boyer_moore('GATTATT')
        occurrences, alignment_tried, num_char_comparisons = bm.query(t)
        self.assertEqual(
            alignment_tried,
            2
        )

    def test_boyer_moore3(self):
        p = 'word'
        t = 'there would have been a time for such a word'
        lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
        bm = self.get_boyer_moore(p, alphabet=lowercase_alphabet)
        occurrences, alignment_tried, num_char_comparisons = bm.query(t)
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
        occurrences, alignment_tried, num_character_comparisons = bm.query(t)

        self.assertListEqual(
            occurrences,
            [0, 19]
        )

        self.assertEqual(
            alignment_tried,
            5
        )

        self.assertEqual(num_character_comparisons, 18)


class TestApproximateBoyerMoore(TestCase):
    def test_approximate_boyer_moore1(self):
        p = 'AACTTG'
        t = 'CACTTAATTTG'
        matches, num_index_hits = PigeonHoleApproximateMatching().query_bm(p=p, t=t, m=2)
        self.assertEqual(matches, [0, 5])

    def test_exact_kmer_index2(self):
        p = 'TCTA'
        t = 'GCTACGATCTAGAATCTA'
        matches, num_index_hits = PigeonHoleApproximateMatching().query_bm(t=t, m=0, p=p)
        self.assertEqual(matches, [7, 14])


    def test_approximate_kmer_index3(self):
        p = 'AACTTG'
        t = 'CACTTAATTTG'
        matches, num_index_hits = PigeonHoleApproximateMatching().query_bm(t=t, m=2, p=p)
        self.assertEqual(matches, [0, 5])


class TestApproximateSubSeqIndex(TestCase):
    def test_approximate_subseq_index1(self):
        p = 'AACTTG'
        t = 'CACTTAATTTG'
        matches, num_index_hits = PigeonHoleApproximateMatching().query_subseq_index(t=t, m=2, p=p, k=2)
        self.assertEqual(matches, [0, 5])

    def test_exact_subseq_index2(self):
        p = 'TCTA'
        t = 'GCTACGATCTAGAATCTA'
        matches, num_index_hits = PigeonHoleApproximateMatching().query_subseq_index(t=t, m=0, p=p, k=2)
        self.assertEqual(matches, [7, 14])

    def test_approximate_subseq_index3(self):
        p = 'AACTTG'
        t = 'CACTTAATTTG'
        matches, num_index_hits = PigeonHoleApproximateMatching().query_subseq_index(t=t, m=2, p=p, k=2)
        self.assertEqual(matches, [0, 5])

    def test_approximate_subseq_index4(self):
        t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
        p = 'to-morrow and to-morrow '
        matches, num_index_hits = PigeonHoleApproximateMatching().query_subseq_index(t=t, m=0, p=p, k=8, ival=3)
        self.assertEqual(matches, [0, 14])
        self.assertEqual(num_index_hits, 2)


class TestApproximateMatchingParity(TestCase):

    def test_equality(self):
        t = 'ACTTACTTGATAAAGT'
        p = 'ACTTTA'
        naive_matches, _ = NaiveAlignment.naive_alignments(p, t, allow_mismatches=True,
                                                           num_mismatches=2, include_reverse_compliment=False)
        bm_matches, _ = PigeonHoleApproximateMatching().query_bm(t=t, m=2, p=p)
        subseq_matches, _ = PigeonHoleApproximateMatching().query_subseq_index(t=t, m=2, p=p, k=2)

        assert len(naive_matches) == len(bm_matches) == len(subseq_matches)


class TestEditDistance(TestCase):

    def test_recursive_1(self):
        a = 'Hi there how are you doing today'
        b = 'Hi there how are you doing today Jim?'
        ed = EditDistance()
        d = ed.ed_recursive(a, b)
        self.assertEqual(d, 5)
        self.assertEqual(ed.calls, 1184)

    def test_dp_1(self):
        a = 'Hi there how are you doing today'
        b = 'Hi there how are you doing today Jim?'
        ed = EditDistance()
        d = ed.ed_dp(a, b)
        self.assertEqual(d, 5)
        self.assertEqual(ed.calls, 1184)

    def test_dp_approx_1(self):
        a = 'GCGTATGC'
        b = 'TATTGGCTATACGGTT'
        ed = ApproximateMatching()
        d = ed.closest_match(a, b)
        self.assertEqual(d, 2)
        self.assertEqual(ed.calls, 128)

    def test_dp_approx_2(self):
        a = 'Hi there how are you doing today'
        b = 'Hi there how are you doing today Jim?'
        ed = ApproximateMatching()
        d = ed.closest_match(a, b)
        self.assertEqual(d, 0)
        self.assertEqual(ed.calls, 1184)

    def test_dp_approx_3(self):
        a = 'When is the store opening today?'
        b = 'When as the shtore opening today.'
        ed = ApproximateMatching()
        d = ed.closest_match(a, b)
        self.assertEqual(d, 3)
        self.assertEqual(ed.calls, 1056)


class TestOverlap(TestCase):

    def test_overlap_1(self):
        a = 'Hello there tom'
        b = 'there tom is sad'
        o = Overlap().overlap(a, b, 5)
        self.assertEqual(o, 9)

    def test_overlap_2(self):
        a = 'what is'
        b = 'is the location'
        o = Overlap().overlap(a, b, 5)
        self.assertEqual(o, 0)


class TestIdentifyMatchExactPrefixSuffix(TestCase):
    def test_overlap_1(self):
        reads = ['ABCDEFG', 'EFGHIJ', 'HIJABC']
        o = MatchExactPrefixSuffix().get_overlapping_pairs(reads, 3, return_pairs=True)
        self.assertEqual(o, [('ABCDEFG', 'EFGHIJ'), ('EFGHIJ', 'HIJABC'), ('HIJABC', 'ABCDEFG')])

    def test_overlap_2(self):
        reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
        o = MatchExactPrefixSuffix().get_overlapping_pairs(reads, 4, return_pairs=True)

        self.assertEqual(o,
             [('CGTACG', 'TACGTA'), ('CGTACG', 'TACGAT'),
              ('TACGTA', 'CGTACG'), ('GTACGT', 'ACGTAC'),
              ('ACGTAC', 'GTACGT'), ('ACGTAC', 'GTACGA')]
         )


class TestIdentifyOverlapSuffixContainedInRead(TestCase):
    def test_overlap_1(self):
        reads = ['ABCDEFG', 'EFGHIJ', 'HIJABC']
        o = OverlapSuffixContainedInRead().get_overlapping_pairs(reads, 3, return_pairs=True)
        self.assertEqual(o, [('ABCDEFG', 'EFGHIJ'), ('EFGHIJ', 'HIJABC'), ('HIJABC', 'ABCDEFG')])

    def test_overlap_2(self):
        reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
        o = OverlapSuffixContainedInRead().get_overlapping_pairs(reads, 4, return_pairs=True)

        self.assertSetEqual(set(o),
             set([('CGTACG', 'TACGTA'),
              ('CGTACG', 'GTACGT'),
              ('CGTACG', 'GTACGA'),
              ('CGTACG', 'TACGAT'),
              ('TACGTA', 'ACGTAC'),
              ('TACGTA', 'CGTACG'),
              ('GTACGT', 'TACGTA'),
              ('GTACGT', 'ACGTAC'),
              ('ACGTAC', 'GTACGA'),
              ('ACGTAC', 'GTACGT'),
              ('ACGTAC', 'CGTACG'),
              ('GTACGA', 'TACGAT')])
         )

    def test_overlap_3(self):
        reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
        o = OverlapSuffixContainedInRead().get_overlapping_pairs(reads, 5, return_pairs=True)

        self.assertSetEqual(set(o),
             set([('CGTACG', 'GTACGT'),
                 ('CGTACG', 'GTACGA'),
                 ('TACGTA', 'ACGTAC'),
                 ('GTACGT', 'TACGTA'),
                 ('ACGTAC', 'CGTACG'),
                 ('GTACGA', 'TACGAT')]
                 )
         )

    def test_overlap_4(self):
        reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT', 'GGACGC', 'TACCGC']

        k = 5
        expected_results = []
        for a, b in product(reads, reads):
            if a != b:
                overlap = OverlapSuffixContainedInRead().overlap(a, b, k)
                if overlap >= k:
                    expected_results.append((a, b))

        o = OverlapSuffixContainedInRead().get_overlapping_pairs(reads, 5, return_pairs=True)

        self.assertSetEqual(set(o),
             set(expected_results)
         )


class TestBruteForceSCS(TestCase):
    def test_scs_1(self):
        strings = ['ABC', 'BCA', 'CAB']
        shortest_superstrings = ShortestCommonSuperstring().brute_force_scs(strings, k=1)
        self.assertEqual(len(shortest_superstrings), 3)
        self.assertSetEqual(set(shortest_superstrings), {'ABCAB', 'BCABC', 'CABCA'})

    def test_scs_2(self):
        strings = ['GAT', 'TAG', 'TCG', 'TGC', 'AAT', 'ATA']
        shortest_superstrings = ShortestCommonSuperstring().brute_force_scs(strings, k=1)
        self.assertSetEqual(set(shortest_superstrings), {'AATAGATCGTGC',
         'AATAGATGCTCG',
         'AATAGTCGATGC',
         'AATCGATAGTGC',
         'AATGCTCGATAG',
         'TCGAATAGATGC',
         'TCGATAGAATGC',
         'TCGATGCAATAG',
         'TGCAATAGATCG',
         'TGCAATCGATAG'}
         )


class TestGreedySCS(TestCase):
    """ For these simple problems, the greedy SCS should be one of the true SCSs"""
    def test_scs_1(self):
        strings = ['ABC', 'BCA', 'CAB']
        shortest_superstring = ShortestCommonSuperstring().greedy_scs(strings, k=1)
        self.assertIn(shortest_superstring, {'ABCAB', 'BCABC', 'CABCA'})

    def test_scs_2(self):
        strings = ['GAT', 'TAG', 'TCG', 'TGC', 'AAT', 'ATA']
        shortest_superstring = ShortestCommonSuperstring().greedy_scs(strings, k=1)
        self.assertIn(shortest_superstring, {'AATAGATCGTGC',
         'AATAGATGCTCG',
         'AATAGTCGATGC',
         'AATCGATAGTGC',
         'AATGCTCGATAG',
         'TCGAATAGATGC',
         'TCGATAGAATGC',
         'TCGATGCAATAG',
         'TGCAATAGATCG',
         'TGCAATCGATAG'}
         )

    def test_scs_kmer_index_1(self):
        strings = ['ABC', 'BCA', 'CAB']
        shortest_superstring = ShortestCommonSuperstring().greedy_scs(strings, k=1, kmer_index_k=1)
        self.assertIn(shortest_superstring, {'ABCAB', 'BCABC', 'CABCA'})

        shortest_superstring = ShortestCommonSuperstring().greedy_scs(strings, k=1, kmer_index_k=2)
        self.assertIn(shortest_superstring, {'ABCAB', 'BCABC', 'CABCA'})

    def test_scs_kmer_index_2(self):
        strings = ['GAT', 'TAG', 'TCG', 'TGC', 'AAT', 'ATA']
        shortest_superstring = ShortestCommonSuperstring().greedy_scs(strings, k=1, kmer_index_k=1)
        self.assertIn(
            shortest_superstring,
            {'AATAGATCGTGC',
             'AATAGATGCTCG',
             'AATAGTCGATGC',
             'AATCGATAGTGC',
             'AATGCTCGATAG',
             'TCGAATAGATGC',
             'TCGATAGAATGC',
             'TCGATGCAATAG',
             'TGCAATAGATCG',
             'TGCAATCGATAG'
             })
