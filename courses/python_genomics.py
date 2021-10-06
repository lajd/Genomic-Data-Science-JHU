import sys
from typing import Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from collections import OrderedDict, defaultdict, Counter, deque

DATA_PATH = '/Users/jonathan/PycharmProjects/genomics_data_science/data/dna2.fasta'


class FastaMetadata:
    @staticmethod
    def get_fasta_info(data_path: str = DATA_PATH) -> dict:
        metadata = {}

        with open(data_path, 'r') as handle:
            for i, record in enumerate(SeqIO.parse(handle, "fasta")):
                metadata[record.id] = {
                    'record_length': len(record)
                }
        return metadata


class ORFs:
    START_CODONS = {'ATG'}
    STOP_CODONS = {'TAA', 'TAG', 'TGA'}

    @staticmethod
    def get_orfs(sequence_record: SeqRecord) -> dict:
        map_orf_to_ranges = defaultdict(list)

        for reading_frame in (1, 2, 3):
            # Get the indices of all the start/stop codons
            start_codon_indices = deque()
            stop_codon_indices = deque()

            pointer = reading_frame - 1
            codon_index = 0
            while pointer < len(sequence_record) - 2:
                codon = sequence_record[pointer: pointer + 3].seq

                if codon in ORFs.START_CODONS:
                    start_codon_indices.append(codon_index)
                elif codon in ORFs.STOP_CODONS:
                    stop_codon_indices.append(codon_index)
                pointer += 3
                codon_index += 1

            # ORFs are assumed to be the longest sequence bounded by a start/stop codon
            stop_codon_idx = -1
            while start_codon_indices and stop_codon_indices:
                while start_codon_indices:
                    start_codon_idx = start_codon_indices.popleft()
                    if start_codon_idx > stop_codon_idx:
                        break

                while stop_codon_indices:
                    stop_codon_idx = stop_codon_indices.popleft()
                    if stop_codon_idx > start_codon_idx:
                        break
                if start_codon_idx < stop_codon_idx:
                    orf_start_idx = start_codon_idx * 3
                    orf_end_idx = stop_codon_idx * 3 + 3
                    map_orf_to_ranges[reading_frame].append(
                        {
                            # Indices are 0-based, meant for python slicing
                            'orf_start_idx': orf_start_idx,
                            'orf_end_idx': orf_end_idx,
                            'orf_length': orf_end_idx - orf_start_idx
                        }
                    )
        return map_orf_to_ranges

    @staticmethod
    def get_all_orfs_from_fasta(data_path: str = DATA_PATH, sequence_id: Optional[str] = None) -> dict:
        orfs_metadata = {}
        with open(data_path, 'r') as handle:
            for i, record in enumerate(SeqIO.parse(handle, "fasta")):
                if sequence_id and record.id != sequence_id:
                    continue
                orfs = ORFs.get_orfs(record)
                orfs_metadata[record.id] = orfs
        return orfs_metadata

    @staticmethod
    def get_all_orfs_with_reading_frame(reading_frame: int, data_path: str = DATA_PATH, sequence_id: Optional[str] = None):
        all_orfs = ORFs.get_all_orfs_from_fasta(data_path, sequence_id=sequence_id)
        all_orfs_of_reading_frame = []
        for record_orf in all_orfs.values():
            all_orfs_of_reading_frame.extend(record_orf[reading_frame])
        return all_orfs_of_reading_frame


class RepeatCounter:

    @staticmethod
    def count_repeats(sequence_record: SeqRecord, n: int) -> tuple:
        pointer = 0
        counter = Counter()
        repeat_metadata = defaultdict(list)
        while pointer < len(sequence_record) - (n - 1):
            substring = sequence_record[pointer: pointer + n].seq
            counter[substring] += 1
            repeat_metadata[substring].append({'start_idx': pointer, 'end_idx': pointer + n})
            pointer += 1
        return counter, repeat_metadata

    @staticmethod
    def count_all_repeats(n: int, data_path: str = DATA_PATH, sequence_id: Optional[str] = None):
        global_counter = Counter()
        with open(data_path, 'r') as handle:
            for i, record in enumerate(SeqIO.parse(handle, "fasta")):
                if sequence_id and record.id != sequence_id:
                    continue
                counter, repeat_meta = RepeatCounter.count_repeats(record, n)
                global_counter.update(counter)
        return global_counter


if __name__ == '__main__':

    # Test 1
    orfs = ORFs.get_orfs(SeqRecord('ATGAAATAG'))
    assert len(orfs) == 1
    assert len(orfs[1]) == 1
    assert orfs[1][0]['orf_length'] == 9

    orfs = ORFs.get_orfs(SeqRecord('ATGAACATCATGAAATAG'))
    assert len(orfs) == 1
    assert len(orfs[1]) == 1
    assert orfs[1][0]['orf_length'] == 18

    orfs = ORFs.get_orfs(SeqRecord('ATGAACTAGATCATGAAATAG'))
    assert len(orfs) == 1
    assert len(orfs[1]) == 2
    assert orfs[1][0]['orf_length'] == 9
    assert orfs[1][1]['orf_length'] == 9

    orfs = ORFs.get_orfs(SeqRecord('GATGAACTAGATCATGAAATAG'))
    assert len(orfs) == 1
    assert len(orfs[2]) == 2
    assert orfs[2][0]['orf_length'] == 9
    assert orfs[2][1]['orf_length'] == 9

    orfs = ORFs.get_orfs(SeqRecord('ATGAACTAGATGATCATGAAATAGCATGAAATAG'))
    assert len(orfs) == 2
    assert len(orfs[1]) == 2
    assert len(orfs[2]) == 1
    assert orfs[1][0]['orf_length'] == 9
    assert orfs[1][1]['orf_length'] == 15
    assert orfs[2][0]['orf_length'] == 9

    # Test 2
    repeats_counter, repeats_meta = RepeatCounter.count_repeats(SeqRecord('ACACA'), 3)
    assert repeats_counter['ACA'] == 2
    assert repeats_meta['ACA'][0]['start_idx'] == 0
    assert repeats_meta['ACA'][0]['end_idx'] == 3



