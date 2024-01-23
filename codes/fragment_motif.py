#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import time
import argparse
import multiprocessing
from collections import Counter
from pyfaidx import Fasta
import itertools


def generate_dna_combinations():
    """
    Generates all combinations of DNA sequences of length 6 using nucleotides A, T, G, and C.

    :return: List of all DNA sequence combinations.
    """
    nucleotides = ['A', 'T', 'G', 'C']
    return [''.join(p) for p in itertools.product(nucleotides, repeat=6)]


def get_reverse_complimentary(seq):
    """
    Generates the reverse complement of a DNA sequence.

    :param seq: String representing a DNA sequence.
    :return: Reverse complement of the input DNA sequence.
    """
    rc_dict = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return ''.join(rc_dict[i] for i in reversed(seq))


def process_batch(batch, reference_path, combinations_set):
    """
    Processes a batch of lines from the BED file, counting motifs.

    :param batch: List of lines from the BED file.
    :param reference_path: Reference sequence path.
    :param combinations_set: Set of all valid DNA sequence combinations.
    :return: Two Counters, one for BPM and one for EDM motifs.
    """
    ref = Fasta(reference_path)
    bpm_counter = Counter({comb: 0 for comb in combinations_set})
    edm_counter = Counter({comb: 0 for comb in combinations_set})
    CHROMS = set(["chr" + str(i) for i in range(1, 23)] + ["chrX"])

    for line in batch:
        # Process each line to count motifs
        chrom, start, end = line.split()[:3]
        start, end = int(start), int(end)

        if chrom not in CHROMS:
            continue

        try:
            front_bpm_motif = ref[chrom][start - 3: start + 3].seq.upper()
            back_bpm_motif = get_reverse_complimentary(ref[chrom][end - 3: end + 3].seq.upper())
            front_edm_motif = ref[chrom][start: start + 6].seq.upper()
            back_edm_motif = get_reverse_complimentary(ref[chrom][end - 6: end].seq.upper())
        except ValueError:
            # Skipping positions beyond chromosome ranges
            continue

        if "N" not in front_bpm_motif and "N" not in back_bpm_motif and "N" not in front_edm_motif and "N" not in back_edm_motif:
            bpm_counter[front_bpm_motif] += 1
            bpm_counter[back_bpm_motif] += 1
            edm_counter[front_edm_motif] += 1
            edm_counter[back_edm_motif] += 1

    return bpm_counter, edm_counter


def process_batch_wrapper(args):
    """
    Wrapper for process_batch to unpack arguments.

    :param args: Tuple containing batch, reference_path and combinations_set.
    :return: Result of process_batch.
    """
    return process_batch(*args)


def read_fixmate_bed_parallel(fixmate_bed, num_processes, reference_path, combinations_set):
    """
    Reads and processes the fixmate BED file in parallel.

    :param fixmate_bed: Path to the fixmate BED file.
    :param num_processes: Number of processes to use for parallel processing.
    :param reference_path: Reference sequence path.
    :param combinations_set: Set of all valid DNA sequence combinations.
    :return: Two Counters, one for BPM and one for EDM motifs.
    """
    with open(fixmate_bed, "r") as fh:
        lines = fh.readlines()

    # Split the lines into batches for parallel processing
    batch_size = len(lines) // num_processes
    line_batches = [lines[i:i + batch_size] for i in range(0, len(lines), batch_size)]

    final_bpm_counter = Counter()
    final_edm_counter = Counter()

    with multiprocessing.Pool(processes=num_processes) as pool:
        # 创建一个元组列表，每个元组包含一个批次和combinations_set
        tasks = [(batch, reference_path, combinations_set) for batch in line_batches]
        results = pool.imap(process_batch_wrapper, tasks)
        for bpm_counter, edm_counter in results:
            final_bpm_counter.update(bpm_counter)
            final_edm_counter.update(edm_counter)

    return final_bpm_counter, final_edm_counter


def write_output(ofname, counter):
    """
    Writes the motif count to an output file.

    :param ofname: The output file name.
    :param counter: The counter containing motif counts.
    """
    with open(ofname, "w") as fw:
        for motif in counter.keys():
            fw.write(f"{motif}\t{counter[motif]}\n")


def main():
    """
    Main function to process the BED file and reference sequence.
    Utilizes multiprocessing for parallel data processing.
    Note: Initial startup may be slow due to process initialization overhead.
    """
    start_time = time.time()

    parser = argparse.ArgumentParser(description='Process BED file and reference sequence')
    parser.add_argument('-i', '--input', required=True, help='Path to the input fixmate BED file')
    parser.add_argument('-r', '--reference', default="/public/database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta",
                        help='Path to the reference sequence file')
    parser.add_argument('-o', '--output', required=True, help='Output file prefix')
    parser.add_argument('-p', '--num_processes', type=int, default=8, help='Number of processes, default is 8')
    args = parser.parse_args()

    # generate combinations
    combinations = generate_dna_combinations()
    combinations_set = set(combinations)

    # Process BED file in parallel
    bpm_counter, edm_counter = read_fixmate_bed_parallel(args.input, args.num_processes, args.reference, combinations_set)

    # Write output files
    write_output(args.output + ".BPM.txt", bpm_counter)
    write_output(args.output + ".EDM.txt", edm_counter)


if __name__ == "__main__":
    main()
