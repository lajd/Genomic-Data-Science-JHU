import os
import random
from typing import Optional
from collections import defaultdict, deque
from matplotlib import pyplot as plt
from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_1 import readGenome

from data import DATA_DIR

from courses.L02_algorithms_for_dna_sequencing.utils.boyer_moore_preproc import BoyerMoore as BoyerMoorePreprocessing


### Boyer-Moore basics
class BoyerMoore:
    def __init__(self, p: str, p_bm: BoyerMoorePreprocessing):
        self.p = p
        self.p_bm = p_bm

    def boyer_moore_alignment(self, t: str):
        occurrences = []
        skipped_alignments = []
        num_char_comparisons = 0
        alignment_start_idx = 0

        while alignment_start_idx <= (len(t) - len(self.p)):
            match = True
            # Check from right to left
            matching_suffix = deque()
            for j in range(len(self.p) - 1, -1, -1):
                num_char_comparisons += 1
                k = alignment_start_idx + j
                if self.p[j] == t[k]:
                    matching_suffix.appendleft(self.p[j])
                    continue
                else:
                    match = False
                    bad_char_shift = self.p_bm.bad_character_rule(j, t[k])
                    good_suffix_shift = self.p_bm.good_suffix_rule(j)
                    shift = max(bad_char_shift, good_suffix_shift, 1)
                    if shift > 1:
                        skipped_alignments.append(
                            {
                                'shift': shift - 1,
                                'bad_char_shift':  bad_char_shift - 1 if bad_char_shift > 0 else 0,
                                'good_suffix_shift': good_suffix_shift - 1 if good_suffix_shift > 0 else 0
                            }
                        )
                    alignment_start_idx += shift
                    break
            if match:
                occurrences.append(alignment_start_idx)
                alignment_start_idx += len(self.p)
        return occurrences, skipped_alignments, num_char_comparisons

