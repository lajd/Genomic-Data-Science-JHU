from unittest import TestCase

from Bio.SeqRecord import SeqRecord

from courses.L01_python.python_genomics import ORFs, RepeatCounter


class TestPythonGenomics(TestCase):

    def test_orfs_1(self):
        # Test 1
        orfs = ORFs.get_orf_ranges(SeqRecord('ATGAAATAG'))
        self.assertEqual(len(orfs), 1)
        self.assertEqual(len(orfs[1]), 1)
        self.assertEqual(orfs[1][0]['orf_length'], 9)

    def test_orfs_2(self):
        # Test 1
        orfs = ORFs.get_orf_ranges(SeqRecord('ATGAACATCATGAAATAG'))
        self.assertEqual(len(orfs), 1)
        self.assertEqual(len(orfs[1]), 1)
        self.assertEqual(orfs[1][0]['orf_length'], 18)

    def test_orfs_3(self):
        # Test 1
        orfs = ORFs.get_orf_ranges(SeqRecord('ATGAACTAGATCATGAAATAG'))
        self.assertEqual(len(orfs), 1)
        self.assertEqual(len(orfs[1]), 2)
        self.assertEqual(orfs[1][0]['orf_length'], 9)
        self.assertEqual(orfs[1][1]['orf_length'], 9)

    def test_orfs_4(self):
        # Test 1
        orfs = ORFs.get_orf_ranges(SeqRecord('ATGAACTAGATCATGAAATAG'))
        self.assertEqual(len(orfs), 1)
        self.assertEqual(len(orfs[1]), 2)
        self.assertEqual(orfs[1][0]['orf_length'], 9)
        self.assertEqual(orfs[1][1]['orf_length'], 9)

    def test_orfs_5(self):
        orfs = ORFs.get_orf_ranges(SeqRecord('ATGAACTAGATGATCATGAAATAGCATGAAATAG'))
        assert len(orfs) == 2
        assert len(orfs[1]) == 2
        assert len(orfs[2]) == 1
        assert orfs[1][0]['orf_length'] == 9
        assert orfs[1][1]['orf_length'] == 15
        assert orfs[2][0]['orf_length'] == 9

    def test_repeat_counter(self):
        repeats_counter, repeats_meta = RepeatCounter.count_repeats(SeqRecord('ACACA'), 3)
        self.assertEqual(repeats_counter['ACA'], 2)
        self.assertEqual(repeats_meta['ACA'][0]['start_idx'], 0)
        self.assertEqual(repeats_meta['ACA'][0]['end_idx'], 3)

