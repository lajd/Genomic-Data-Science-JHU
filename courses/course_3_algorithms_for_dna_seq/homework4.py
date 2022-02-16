import time
import os

from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_1 import readGenome, DATA_DIR, readFastq
from courses.L02_algorithms_for_dna_sequencing.algorithms_for_dna_sequencing_week_4 import ShortestCommonSuperstring

from data import DATA_DIR


if __name__ == '__main__':

    """ Homework """

    mystery_genome = readFastq(os.path.join(DATA_DIR, 'ads1_week4_reads.fq'))[0]
    s1 = ["CCT", "CTT", "TGC", "TGG", "GAT", "ATT"]

    print('Problem 1')
    shortest_superstrings = ShortestCommonSuperstring().brute_force_scs(s1, k=1)
    print('SCS has length %s', len(shortest_superstrings[0]))

    print('Problem 2')
    shortest_superstrings = ShortestCommonSuperstring().brute_force_scs(s1, k=1)
    print('Number of SCSs is %s', len(shortest_superstrings))

    print('Problem 3')
    # Iterate over possible values of k, from largest to smallest, until
    # the desired genome size is obtained
    k = 100
    t1 = time.time()
    scs = ShortestCommonSuperstring()
    while True:
        shortest_superstring = scs.greedy_scs(mystery_genome, k=k, kmer_index_k=k)
        if len(shortest_superstring) == 15894:
            print("Number of As is {}".format(shortest_superstring.count('A')))
            print("Number of Ts is {}".format(shortest_superstring.count('T')))
            print("superstring is: {}".format(shortest_superstring))
            print('k is {}'.format(k))
            break
        k -= 1
    t2 = time.time()
    print('Time to compute is %s s', round(t2-t1))
