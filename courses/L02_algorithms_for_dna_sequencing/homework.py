import os

from Levenshtein.StringMatcher import StringMatcher

from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_1 import readGenome, DATA_DIR, NaiveAlignment
from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_2 import BoyerMooreExact, PigeonHoleApproximateMatching


if __name__ == '__main__':

    """ Homework """
    genome = readGenome(os.path.join(DATA_DIR, 'chr1.GRCh38.excerpt.fasta'))

    p1 = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
    p2 = 'GGCGCGGTGGCTCACGCCTGTAAT'

    print('Problem 1')
    num_alignments = len(genome) - len(p1) + 1
    print('Num alignments in exact matching %s', num_alignments)

    print('Problem 2')
    occurrences, num_comparisons = NaiveAlignment.naive_alignments(p1, genome, include_reverse_compliment=False)
    print('Num comparisons in exact matching %s', num_comparisons)

    print('Problem 3')
    bm = BoyerMooreExact(p1)
    occurrences, num_alignments_tried, num_comparisons = bm.query(genome)
    print('Occurrences %s', occurrences)
    print('Num tried alignments %s', num_alignments_tried)

    print('Problem 4')
    matcher = PigeonHoleApproximateMatching()
    occurrences, _ = matcher.query_bm(p=p2, t=genome, m=2)
    print('Num occurrences %s', len(occurrences))

    print('Problem 5')
    approximate_matcher = PigeonHoleApproximateMatching()
    subseq_occurrences, total_index_hits = approximate_matcher.query_subseq_index(p=p2, t=genome, m=2, k=8, ival=1)
    print('Total index hits count %s', total_index_hits)
    print('Total matches: %s', len(occurrences))


    print('Problem 6')
    approximate_matcher = PigeonHoleApproximateMatching()
    occurrences, total_index_hits = approximate_matcher.query_subseq_index(p=p2, t=genome, m=2, k=8, ival=3)
    print('Total index hits count %s', total_index_hits)
    print('Total matches: %s', len(occurrences))
