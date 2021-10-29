import os

from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_1 import readGenome, DATA_DIR, readFastq
from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_3 import ApproximateMatching


if __name__ == '__main__':

    """ Homework """
    human_genome_sample = readGenome(os.path.join(DATA_DIR, 'chr1.GRCh38.excerpt.fasta'))
    phix_reads = readFastq(os.path.join(DATA_DIR, 'ERR266411_1.for_asm.fastq'))

    p1 = 'GCTGATCGATCGTACG'
    p2 = 'GATTTACCAGATTGAG'

    print('Problem 1')
    ed = ApproximateMatching()
    min_distance = ed.closest_match(p1, human_genome_sample)
    print('Closest alignment distance %s', min_distance)

    print('Problem 2')
    ed = ApproximateMatching()
    min_distance = ed.closest_match(p2, human_genome_sample)
    print('Closest alignment distance %s', min_distance)
