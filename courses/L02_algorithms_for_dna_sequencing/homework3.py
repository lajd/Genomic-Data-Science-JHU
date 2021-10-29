import os

from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_1 import readGenome, DATA_DIR, readFastq
from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_3 import ApproximateMatching, OverlapSuffixContainedInRead


if __name__ == '__main__':

    """ Homework """
    human_genome_sample = readGenome(os.path.join(DATA_DIR, 'chr1.GRCh38.excerpt.fasta'))

    # Get only base data
    phix_reads = readFastq(os.path.join(DATA_DIR, 'ERR266411_1.for_asm.fastq'))[0]

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

    print('Problem 3')
    read_overlap = OverlapSuffixContainedInRead()
    pairs = read_overlap.get_overlapping_pairs(phix_reads, k=30, return_pairs=True)
    unique_pairs = set(pairs)
    print('Number of edges (distinct pairs) in overlap graph is %s', len(unique_pairs))

    print('Problem 4')
    read_overlap = OverlapSuffixContainedInRead()
    pairs = read_overlap.get_overlapping_pairs(phix_reads, k=30, return_pairs=True)
    unique_from_pairs = set([i[0] for i in pairs])
    print('Number of reads with suffix in edge %s', len(unique_from_pairs))
