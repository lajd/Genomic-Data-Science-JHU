import os
from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_1 import readGenome, DATA_DIR, NaiveAlignment
from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_2 import BoyerMoore

if __name__ == '__main__':

    """ Homework """
    genome = readGenome(os.path.join(DATA_DIR, 'chr1.GRCh38.excerpt.fasta'))

    p1 = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'

    print('Problem 1')
    num_alignments = len(genome) - len(p1) + 1
    print('Num alignments in exact matching %s', num_alignments)

    print('Problem 2')
    num_alignments, num_comparisons = NaiveAlignment.naive_alignments(p1, genome, include_reverse_compliment=False)
    print('Num comparisons in exact matching %s', num_comparisons)

    print('Problem 2')
    bm = BoyerMoore(p1)
    num_alignments, num_alignments_tried, num_comparisons = bm.boyer_moore_alignment(genome)
    print('Num tried alignments %s', num_alignments_tried)

