#!/usr/bin/env python

""" Unmapped primer directed read overlap variant caller. """

from argparse import ArgumentParser
from Bio import SeqIO
from itertools import (izip, chain, repeat)
from pyfaidx import Fasta
import csv
import logging
import sys

DEFAULT_KMER_THRESHOLD = 0

def parse_args():
    """ Find variants from fastqs via a mapping-free approach."""
    parser = ArgumentParser(description="Find variants from fastqs via a \
        mapping-free approach.""")
    parser.add_argument('--primer_coords', type=str, required=True, \
        help='Primer coordinates in TSV format.')
    parser.add_argument('--primer_sequences', metavar='FILE', type=str, \
        help='Primer base sequences as determined by a primer generating \
        program.')
    parser.add_argument('--kmer_threshold', type=int, \
        default=DEFAULT_KMER_THRESHOLD, \
        help='K-mer length. If set to zero, k-mer test is not performed by \
        undr rover. Defaults to {}'.format(DEFAULT_KMER_THRESHOLD))
    parser.add_argument('--log', metavar='FILE', type=str, \
        help='Logs progress in specified file, defaults to stdout.')
    parser.add_argument('--reference', metavar='FILE', type=str, \
        help='Reference sequences in Fasta format.')
    parser.add_argument('fastqs', nargs='+', type=str, \
        help='Fastq files containing reads.')
    return parser.parse_args()

def get_block_coords(primer_coords_file):
    """ Retrieves the coordinates of the start and end of each block, and
    return as a list."""
    with open(primer_coords_file) as primer_coords:
        return list(csv.reader(primer_coords, delimiter='\t'))

def get_primer_sequences(primer_sequences_file):
    """ Retrieves the sequences of bases for each primer, and return as a
    list."""
    with open(primer_sequences_file) as primer_sequences:
        return list(csv.reader(primer_sequences, delimiter='\t'))

def make_base_sequence(name, bases, qualities):
    """ Take a list of DNA bases and a corresponding list of quality scores
    and return a list of Base objects where the base and score are
    paired together."""
    if len(bases) <= len(qualities):
        return [Base(b, q) for (b, q) in izip(bases, qualities)]
    else:
        logging.warning("In read {} fewer quality scores {} than bases {}" \
            .format(name, len(qualities), len(bases)))
        # we have fewer quality scores than bases
        # pad the end with 0 scores
        return [Base(b, q) for (b, q) in izip(bases, chain(qualities, \
            repeat(0)))]

class Base(object):
    """ A DNA base paired with its quality score."""
    def __init__(self, base, qual):
        self.base = base # a string
        self.qual = qual # an int
    def as_tuple(self):
        """ Return the base and quality score of the base as a tuple."""
        return (self.base, self.qual)
    def __eq__(self, other):
        return self.as_tuple() == other.as_tuple()
    def __str__(self):
        return str(self.as_tuple())
    def __repr__(self):
        return str(self)
    def __hash__(self):
        return hash(self.as_tuple)

def initialise_blocks(args):
    """ Create blocks, initially containing block coordinates, primer sequences
    and an empty dictionary to which reads pairs will be added. The blocks
    themselves are a dictionary at this stage."""
    blocks = {}
    primer_sequences = {}
    block_coords = get_block_coords(args.primer_coords)
    primer_info = get_primer_sequences(args.primer_sequences)
    reference = Fasta(args.reference)
    for primer in primer_info:
        primer_sequences[primer[0]] = primer[1]
    for block in block_coords:
        ref_sequence = reference[block[0]][int(block[1]) - 1:int(block[2])]
        # Actual block, for which the key is the first 20 bases of the forward
        # primer.
        blocks[primer_sequences[block[3]][:20]] = [block[1], block[2], {}, \
            ref_sequence, block[3], block[4], primer_sequences[block[3]], \
            primer_sequences[block[4]]]
        # Reverse primer (not an actual block), contains the key for the forward
        # primer.
        blocks[primer_sequences[block[4]][:20]] = [primer_sequences[block[4]], \
        primer_sequences[block[3]][:20]]
    return blocks

def complete_blocks(args, blocks):
    """ Organise reads into blocks."""
    for fastq_file in args.fastqs:
        for read in SeqIO.parse(fastq_file, 'fastq'):
            read_bases = str(read.seq)
            primer_key = read_bases[:20]
            if primer_key in blocks:
                if len(blocks[primer_key]) == 8:
                    # Forward primer matched.
                    if blocks[primer_key][6] == read_bases[:len(blocks\
                        [primer_key][6])]:
                        if read.id not in blocks[primer_key][2]:
                            blocks[primer_key][2][read.id] = [read]
                        else:
                            blocks[primer_key][2][read.id].append(read)
                elif len(blocks[primer_key]) == 2:
                    # Reverse primer matched.
                    if blocks[primer_key][0] == read_bases[:len(blocks\
                        [primer_key][0])]:
                        forward_key = blocks[primer_key][1]
                        if read.id not in blocks[forward_key][2]:
                            blocks[forward_key][2][read.id] = [read]
                        else:
                            blocks[forward_key][2][read.id].append(read)
    # For the next stage, we take only the actual blocks.
    return [b[:4] for b in blocks.values() if len(b) > 2]

def process_blocks(args, blocks):
    """ Variant calling stage. Process blocks one at a time and call variants
    for each block."""
    for block_info in blocks:
        count = 0
        for x in block_info[2].values():
            if len(x) == 2:
                count += 1
        print block_info[0], count

def main():
    """ Main function."""
    args = parse_args()
    if args.log is None:
        logfile = sys.stdout
    else:
        logfile = args.log
    logging.basicConfig(filename=logfile, level=logging.DEBUG, filemode='w', \
        format='%(asctime)s %(message)s', datefmt='%a, %d %b %Y %H:%M:%S')
    logging.info('Program started.')
    logging.info('Command line: {0}'.format(' '.join(sys.argv)))
    blocks = initialise_blocks(args)
    final_blocks = complete_blocks(args, blocks)
    process_blocks(args, final_blocks)

if __name__ == '__main__':
    main()
