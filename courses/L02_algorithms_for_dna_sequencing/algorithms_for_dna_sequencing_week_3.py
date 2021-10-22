import os
from abc import ABC, abstractmethod
import random
from typing import Optional, List, Type
import bisect
from collections import defaultdict, deque

from matplotlib import pyplot as plt
import numpy as np
from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_1 import readGenome

from data import DATA_DIR

from courses.L02_algorithms_for_dna_sequencing.utils.boyer_moore_preproc import BoyerMoorePreprocessing


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

    def ed_dp(self, a: str, b: str, approximate_matching: bool = False):
        """

        :param a:
        :param b:
        :param approximate_matching: Whether to use the approximate matching alrgorithm, which
            does not have a bias for any particular offset
        :return:
        """
        d = np.zeros((len(a) + 1, len(b) + 1))
        d[:, 0] = np.arange(0, len(a) + 1)
        if not approximate_matching:
            # Leave first row as all zeros
            d[0, :] = np.arange(0, len(b) + 1)

        for i in range(1, len(a) + 1):
            for j in range(1, len(b) + 1):
                self.calls += 1
                d_right = d[i][j-1] + 1
                d_down = d[i - 1][j] + 1
                delta = 0 if a[i-1] == b[j-1] else 1
                d_diag = d[i-1][j-1] + delta
                d[i][j] = min(d_right, d_down, d_diag)

        if approximate_matching:
            return min(d[-1, :])
            pass
        else:
            return d[-1, -1]


class ApproximateMatching(EditDistance):
    """
    Identify the match of P in T with the fewest number of edits

    The "match" is a substring of T that has the same length as P
    """
    def __init__(self):
        super().__init__()
        pass
