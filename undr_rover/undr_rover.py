#!/usr/bin/env python

""" Unmapped primer directed read overlap variant caller. """

from argparse import ArgumentParser
from Bio import (pairwise2, SeqIO)
from itertools import takewhile
from operator import itemgetter
from pyfaidx import Fasta
import csv
import datetime
import logging
import os
import sys
import vcf

DEFAULT_ABSOLUTE_THRESHOLD = 2
DEFAULT_MAX_VARIANTS = 25
DEFAULT_PRIMER_BASES = 3
DEFAULT_PROPORTION_THRESHOLD = 0.05
OUTPUT_HEADER = '\t'.join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", \
    "FILTER", "INFO"])

def parse_args():
    """ Find variants from fastqs via a mapping-free approach."""
    parser = ArgumentParser(description="Find variants from fastqs via a \
        mapping-free approach.")
    parser.add_argument('--primer_coords', type=str, required=True, \
        help='Primer coordinates in TSV format.')
    parser.add_argument('--primer_sequences', metavar='FILE', type=str, \
        required=True, help='Primer base sequences as determined by a primer \
        generating program.')
    parser.add_argument('--primer_bases', type=int, \
        default=DEFAULT_PRIMER_BASES, \
        help='Number of bases from primer region to use in gapped alignment.' \
        'Helps with variant calling near the edges of a block.')
    parser.add_argument('--proportionthresh', metavar='N', type=float, \
        default=DEFAULT_PROPORTION_THRESHOLD, \
        help='Keep variants which appear in this proportion of the read pairs '
             'for a given target region, and bin otherwise. '
             'Defaults to {}.'.format(DEFAULT_PROPORTION_THRESHOLD))
    parser.add_argument('--absthresh', metavar='N', type=int, \
        default=DEFAULT_ABSOLUTE_THRESHOLD, \
        help='Only keep variants which appear in at least this many \
        read pairs. ' \
        'Defaults to {}.'.format(DEFAULT_ABSOLUTE_THRESHOLD))
    parser.add_argument('--qualthresh', metavar='N', type=int, \
        help='Minimum base quality score (phred).')
    parser.add_argument('--max_variants', metavar='N', type=int, \
        default=DEFAULT_MAX_VARIANTS, \
        help='Ignore reads with greater than this many variants observed.' \
        'Defaults to {}.'.format(DEFAULT_MAX_VARIANTS))
    parser.add_argument('--reference', metavar='FILE', type=str, \
        required=True, help='Reference sequences in Fasta format.')
    parser.add_argument('--id_info', metavar='FILE', type=str, \
    help='File containing rs ID information in VCF format.')
    parser.add_argument('--out', metavar='FILE', type=str, \
        required=True, help='Name of output file containing called variants.')
    parser.add_argument('--log', metavar='FILE', type=str, \
        help='Logs progress in specified file, defaults to stdout.')
    parser.add_argument('--coverdir', required=False, \
        help='Directory to write coverage files, defaults to current working \
        directory.')
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

def reverse_complement(sequence):
    """ Return the reverse complement of a DNA string."""
    rc_bases = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join([rc_bases[base] for base in sequence[::-1]])

def nts(none_string):
    """ Returns an empty string for None."""
    return none_string or ''

class Base(object):
    """ A DNA base paired with its quality score."""
    def __init__(self, base, qual):
        self.base = base  # a string
        self.qual = qual  # an int
    def __eq__(self, other):
        return self.as_tuple() == other.as_tuple()
    def __hash__(self):
        return hash(self.as_tuple)
    def __repr__(self):
        return str(self)
    def __str__(self):
        return str(self.as_tuple())
    def as_tuple(self):
        """ Return the base and quality score of the base as a tuple."""
        return (self.base, self.qual)

class SNV(object):
    """ Single nucleotide variant. Bases are represented as DNA strings."""
    def __init__(self, chrsm, pos, bases, qual):
        self.chrsm = chrsm
        self.pos = pos
        self.ref_base = bases[0]
        self.seq_base = bases[1]
        self.qual = qual
        self.filter_reason = None
        self.info = []
    def __eq__(self, other):
        return self.as_tuple() == other.as_tuple()
    def __hash__(self):
        return hash(self.as_tuple())
    def __repr__(self):
        return str(self)
    def __str__(self):
        return "S: {} {} {} {}".format(self.chrsm, self.pos, self.ref_base, \
            self.seq_base)
    def as_tuple(self):
        """ Return information about the SNV as a 4-tuple."""
        return (self.chrsm, self.pos, self.ref_base, self.seq_base)
    def position(self):
        """ SNV POS."""
        return self.pos
    def ref(self):
        """ REF base."""
        return self.ref_base
    def alt(self):
        """ ALT base."""
        return self.seq_base
    def fil(self):
        """ Return "PASS" if the SNV is not filtered, or the reason(s) for
        being discarded otherwise."""
        if self.filter_reason is None:
            return "PASS"
        return self.filter_reason[1:]

class Insertion(object):
    """ Insertion. Bases are represented as DNA strings."""
    def __init__(self, chrsm, pos, bases, qual):
        self.chrsm = chrsm
        self.pos = pos
        self.inserted_bases = bases[0]
        self.qual = qual
        self.filter_reason = None
        self.info = []
        self.context = bases[1]
    def __eq__(self, other):
        return self.as_tuple() == other.as_tuple()
    def __hash__(self):
        return hash(self.as_tuple())
    def __repr__(self):
        return str(self)
    def __str__(self):
        return "I: {} {} {}".format(self.chrsm, self.pos, self.inserted_bases)
    def as_tuple(self):
        """ Return information about the insertion as a 3-tuple."""
        return (self.chrsm, self.pos, self.inserted_bases)
    def position(self):
        """ Insertion POS."""
        return self.pos - 1
    def ref(self):
        """ REF base."""
        return self.context
    def alt(self):
        """ ALT (inserted) bases."""
        return self.context + self.inserted_bases
    def fil(self):
        """ Return "PASS" if the Insertion is not filtered, or the reason(s) for
        being discarded otherwise."""
        if self.filter_reason is None:
            return "PASS"
        return self.filter_reason[1:]

class Deletion(object):
    """ Deletion. Bases are represented as DNA strings."""
    def __init__(self, chrsm, pos, bases, qual):
        self.chrsm = chrsm
        self.pos = pos
        self.deleted_bases = bases[0]
        self.qual = qual
        self.filter_reason = None
        self.info = []
        self.context = bases[1]
    def __eq__(self, other):
        return self.as_tuple() == other.as_tuple()
    def __hash__(self):
        return hash(self.as_tuple())
    def __repr__(self):
        return str(self)
    def __str__(self):
        return "D: {} {} {}".format(self.chrsm, self.pos, self.deleted_bases)
    def as_tuple(self):
        """ Return infromation about the deletion as a 3-tuple."""
        return (self.chrsm, self.pos, self.deleted_bases)
    def position(self):
        """ Deletion POS."""
        return self.pos - 1
    def ref(self):
        """ REF (deleted) bases."""
        return self.context + self.deleted_bases
    def alt(self):
        """ ALT base."""
        return self.context
    def fil(self):
        """ Return "PASS" if the Deletion is not filtered, or the reason(s) for
        being discarded otherwise."""
        if self.filter_reason is None:
            return "PASS"
        return self.filter_reason[1:]

def read_snvs(args, chrsm, qual, pos, insert_seq, bases, direction):
    """ Find all the SNV's in a read. 0 if we find two SNV's in a row."""
    check = False
    pos -= args.primer_bases * direction
    result = []
    # Initiate the relevant indices given the direction of the read.
    i = 0 if direction == 1 else -1
    new = slice(1, None) if direction == 1 else slice(None, -1)
    # Check for identical insert sequence and read (no variants).
    if direction == 1 and insert_seq[args.primer_bases:-1 * args.primer_bases] \
    == bases[args.primer_bases:len(insert_seq) - args.primer_bases]:
        return result
    if direction == -1 and insert_seq[args.primer_bases:-1 * args.primer_bases]\
     == bases[-1 * len(insert_seq) + args.primer_bases:-1 * args.primer_bases]:
        return result
    while insert_seq and bases:
        if insert_seq[i] == bases[i]:
            insert_seq = insert_seq[new]
            bases = bases[new]
            qual = qual[new]
            check = False
            pos += direction
        else:
            # If we just saw an SNV, return 0 as there are now two in a row.
            if check is True:
                return 0
            result.append(SNV(chrsm, pos, [insert_seq[i], bases[i]], '.'))
            if not ((args.qualthresh is None) or (qual[i] >= \
            args.qualthresh)):
                result[-1].filter_reason = ''.join([nts(result[-1].\
                    filter_reason), ";qlt"])
            insert_seq = insert_seq[new]
            bases = bases[new]
            qual = qual[new]
            pos += direction
            check = True
    return result

def read_variants(args, chrsm, qual, pos, insert_seq, bases, direction):
    """ Find all the variants in a read (SNVs, Insertions, Deletions)."""
    pos -= args.primer_bases * direction
    result = []
    # Initialise the relevant indices given the direction of the read.
    i = 0 if direction == 1 else -1
    new = slice(1, None) if direction == 1 else slice(None, -1)
    # Pairwise2 uses a C implementation of the Needleman-Wunsch algorithm.
    aligned_insert, aligned_read = pairwise2.align.globalms(insert_seq, bases, \
        2, 0, -2, -1, penalize_end_gaps=(0, 0), one_alignment_only=1)[0][:2]
    context = '-'
    while aligned_insert and aligned_read:
        if aligned_insert[i] == aligned_read[i]:
            context = aligned_insert[i]
            aligned_insert = aligned_insert[new]
            aligned_read = aligned_read[new]
            qual = qual[new]
            pos += direction
        else:
            if not (aligned_insert[i] == '-' or aligned_read[i] == '-'):
                # Single Nucleotide Variation
                result.append(SNV(chrsm, pos, [aligned_insert[i], \
                    aligned_read[i]], '.'))
                if not ((args.qualthresh is None) or (qual[i] >= \
                args.qualthresh)):
                    result[-1].filter_reason = ''.join([nts(result[-1].\
                        filter_reason), ";qlt"])
                context = aligned_insert[i]
                aligned_insert = aligned_insert[new]
                aligned_read = aligned_read[new]
                qual = qual[new]
                pos += direction
            elif aligned_insert[i] == '-':
                # Insertion
                indel_length = sum(1 for _ in takewhile(lambda x: x == '-', \
                    aligned_insert[::direction]))
                skip = slice(indel_length, None) if direction == 1 else \
                slice(None, -1 * indel_length)
                if indel_length >= len(aligned_insert):
                    return result
                if direction == 1:
                    result.append(Insertion(chrsm, pos, \
                    [aligned_read[:indel_length], context], '.'))
                else:
                    result.append(Insertion(chrsm, pos + 1, [aligned_read[-1 * \
                        indel_length:], aligned_insert[-1 * \
                        indel_length - 1]], '.'))
                # insertion with QUAL data?
                if not ((args.qualthresh is None) or all([b >= \
                args.qualthresh for b in qual[skip]])):
                    result[-1].filter_reason = ''.join([nts(result[-1].\
                    filter_reason), ";qlt"])
                aligned_read = aligned_read[skip]
                aligned_insert = aligned_insert[skip]
                qual = qual[skip]
            elif aligned_read[i] == '-':
                # Deletion
                indel_length = sum(1 for _ in takewhile(lambda x: x == '-', \
                    aligned_read[::direction]))
                skip = slice(indel_length, None) if direction == 1 else \
                slice(None, -1 * indel_length)
                if indel_length >= len(aligned_insert):
                    return result
                if direction == 1:
                    result.append(Deletion(chrsm, pos, \
                    [aligned_insert[:indel_length], context], '.'))
                else:
                    result.append(Deletion(chrsm, pos - indel_length + 1, \
                        [aligned_insert[-1 * indel_length:], aligned_insert\
                        [-1 * indel_length - 1]], '.'))
                # deletion with QUAL data?
                if not ((args.qualthresh is None) or (qual[i] >= \
                args.qualthresh)):
                    result[-1].filter_reason = ''.join([nts(result[-1].\
                    filter_reason), ";qlt"])
                context = aligned_insert[indel_length - 1]
                aligned_insert = aligned_insert[skip]
                aligned_read = aligned_read[skip]
                pos += indel_length * direction
    return result

def initialise_blocks(args):
    """ Create blocks, initially containing block coordinates, primer sequences
    and an empty dictionary into which reads pairs will be added. The blocks
    themselves are a dictionary at this stage."""
    blocks = {}
    primer_sequences = {}
    block_coords = get_block_coords(args.primer_coords)
    primer_info = get_primer_sequences(args.primer_sequences)
    reference = Fasta(args.reference)
    for primer in primer_info:
        primer_sequences[primer[0]] = primer[1]
    for block in block_coords:
        # Take some extra bases to help with the variant calling near the edges
        # of the block.
        ref_sequence = reference[block[0]][int(block[1]) - 1 - \
        args.primer_bases:int(block[2]) + args.primer_bases]
        # Actual block, for which the key is the first 20 bases of the forward
        # primer.
        blocks[primer_sequences[block[3]][:20]] = [block[0], block[1], \
        block[2], {}, str(ref_sequence), block[3], block[4], \
        primer_sequences[block[3]], primer_sequences[block[4]]]
        # Reverse primer (not an actual block), contains the key for the forward
        # primer.
        blocks[primer_sequences[block[4]][:20]] = [primer_sequences[block[4]], \
        primer_sequences[block[3]][:20]]
    return blocks

def complete_blocks(args, blocks):
    """ Organise reads into blocks."""
    for fastq_file in args.fastqs:
        base = os.path.basename(fastq_file)
        sample = base.split('_')
        if len(sample) > 0:
            sample = '_'.join(sample[:3])
        else:
            exit('Cannot deduce sample name from fastq filename {}'.\
                format(fastq_file))
        for read in SeqIO.parse(fastq_file, 'fastq'):
            # Try to match each read with an expected primer.
            read_bases = str(read.seq)
            primer_key = read_bases[:20]
            if len(blocks.get(primer_key, [])) == 9:
                # Possible forward primer matched.
                fseq = blocks[primer_key][7]
                if fseq == read_bases[:len(fseq)]:
                    if read.id not in blocks[primer_key][3]:
                        blocks[primer_key][3][read.id] = [read, 0, \
                        len(fseq), 0, sample]
                    else:
                        blocks[primer_key][3][read.id][0] = read
                        blocks[primer_key][3][read.id][2] = len(fseq)
            elif len(blocks.get(primer_key, [])) == 2:
                # Possible reverse primer matched.
                rseq = blocks[primer_key][0]
                if rseq == read_bases[:len(rseq)]:
                    forward_key = blocks[primer_key][1]
                    if read.id not in blocks[forward_key][3]:
                        blocks[forward_key][3][read.id] = [0, read, \
                        0, len(rseq), sample]
                    else:
                        blocks[forward_key][3][read.id][1] = read
                        blocks[forward_key][3][read.id][3] = len(rseq)
    # For the next stage, we only need the actual blocks.
    return [b[:5] for b in blocks.values() if len(b) > 2]

def process_blocks(args, blocks, id_info, vcf_file):
    """ Variant calling stage. Process blocks one at a time and call variants
    for each block."""
    coverage_info = []
    for block_info in blocks:
        block_vars = {}
        num_pairs = 0
        chrsm, start, end, reads, insert_seq = block_info[:5]
        start = int(start)
        end = int(end)
        logging.info("Processing block chr: {}, start: {}, end: {}"\
            .format(chrsm, start, end))
        for read_pair in [r for r in reads.values() if 0 not in r]:
            num_pairs += 1
            read1, read2, fprimerlen, rprimerlen, sample = read_pair
            forward_bases = read1.seq[fprimerlen - args.primer_bases:]
            reverse_bases = read2.seq[rprimerlen - args.primer_bases:]

            forward_qual = read1.letter_annotations['phred_quality']\
            [fprimerlen - args.primer_bases:]
            reverse_qual = read2.letter_annotations['phred_quality']\
            [rprimerlen - args.primer_bases:]

            insert = insert_seq.upper()
            forward_seq = str(forward_bases.upper())
            reverse_seq = reverse_complement(str(reverse_bases)).upper()

            # For both reads, we initially assume that they only have single
            # nucleotide variants. If we detect successive SNV's in those
            # reads, we go back and do the gapped alignment.

            variants1 = read_snvs(args, chrsm, forward_qual, start, \
                insert, forward_seq, 1)
            variants2 = read_snvs(args, chrsm, reverse_qual, end, \
                insert, reverse_seq, -1)

            if variants1 == 0:
                variants1 = read_variants(args, chrsm, forward_qual, \
                    start, insert, forward_seq, 1)
            if variants2 == 0:
                variants2 = read_variants(args, chrsm, reverse_qual, end, \
                    insert, reverse_seq, -1)

            # Ignore reads which have an unusually high amount of variants.
            if len(variants1) > args.max_variants or len(variants2) > \
            args.max_variants:
                variants1, variants2 = [], []
                logging.info("Read {} discarded due to an unusually high \
amount of variants.".format(read1.id))

            # Consider variants each read in the pair share in common.
            for var in set(variants1).intersection(set(variants2)):
                # Only consider variants within the bounds of the block.
                if var.pos >= start and var.pos <= end:
                    block_vars[var] = block_vars.get(var, 0) + 1

        logging.info("Number of read pairs in block: {}".format(num_pairs))
        logging.info("Number of variants found in block: {}".\
            format(len(block_vars)))

        for var in block_vars:
            num_vars = block_vars[var]
            proportion = float(num_vars) / num_pairs
            var.info.extend(["Sample=" + str(sample), "NV=" + str(num_vars), \
                "NP=" + str(num_pairs), "PCT=" + str('{:.2%}'\
                    .format(proportion))])
            if num_vars < args.absthresh:
                var.filter_reason = ''.join([nts(var.filter_reason), ";at"])
            if proportion < args.proportionthresh:
                var.filter_reason = ''.join([nts(var.filter_reason), ";pt"])
            write_variant(vcf_file, var, id_info, args)

        coverage_info.append((chrsm, start, end, num_pairs))

    coverage_filename = sample + '.coverage'
    if args.coverdir is not None:
        coverage_filename = os.path.join(args.coverdir, coverage_filename)
    with open(coverage_filename, 'w') as coverage_file:
        write_coverage_data(coverage_file, coverage_info)

def write_variant(vcf_file, variant, id_info, args):
    """ Writes variant to vcf_file, while also finding the relevant rs number
    from dbsnp if applicable."""
    info = 0
    record_info = None
    # If the variant is deemed a "PASS", find the relevant rs number from dbsnp.
    if variant.fil() == "PASS" and args.id_info:
        for record in id_info.fetch(variant.chrsm, variant.position(), variant.\
        position() + max(len(variant.ref()), len(variant.alt())) + 1):
            if record.POS == variant.position() and record.REF == \
            variant.ref() and (variant.alt() in record.ALT):
                info = 1
                record_info = record
    if info == 1:
        vcf_file.write('\t'.join([variant.chrsm, str(variant.position()), \
str(record_info.ID), variant.ref(), variant.alt(), str(variant.qual), variant\
.fil(), ';'.join(variant.info)]) + '\n')
    else:
        vcf_file.write('\t'.join([variant.chrsm, str(variant.position()), \
'.', variant.ref(), variant.alt(), str(variant.qual), variant.fil(), ';'\
.join(variant.info)]) + '\n')

def write_coverage_data(coverage_file, coverage_info):
    """ Write coverage information to the coverage files."""
    coverage_file.write('chr\tblock_start\tblock_end\tnum_pairs\n')
    for chrsm, start, end, num_pairs in sorted(coverage_info, \
    key=itemgetter(3)):
        coverage_file.write('{}\t{}\t{}\t{}\n'.format(chrsm, start, \
            end, num_pairs))

def write_metadata(args, vcf_file):
    """ Write the opening lines of metadata to the vcf file."""
    vcf_file.write("##fileformat=VCFv4.2" + '\n')
    today = str(datetime.date.today())
    vcf_file.write("##fileDate=" + today[:4] + today[5:7] + today[8:] + '\n')
    vcf_file.write("##source=UNDR ROVER" + '\n')
    vcf_file.write("##INFO=<ID=Sample,Number=1,Type=String,Description=\
\"Sample Name\">" + '\n')
    vcf_file.write("##INFO=<ID=NV,Number=1,Type=Float,Description=\
\"Number of read pairs with variant\">" + '\n')
    vcf_file.write("##INFO=<ID=NP,Number=1,Type=Float,Description=\
\"Number of read pairs at POS\">" + '\n')
    vcf_file.write("##INFO=<ID=PCT,Number=1,Type=Float,Description=\
\"Percentage of read pairs at POS with variant\">" + '\n')
    if args.qualthresh:
        vcf_file.write("##FILTER=<ID=qlt,Description=\"Variant has phred \
quality score below " + str(args.qualthresh) + "\">" + '\n')
    vcf_file.write("##FILTER=<ID=at,Description=\"Variant does not appear \
in at least " + str(args.absthresh) + " read pairs\">" + '\n')
    vcf_file.write("##FILTER=<ID=pt,Description=\"Variant does not appear \
in at least " + str(args.proportionthresh*100) \
+ "% of read pairs for the given region\">" + '\n')

def main():
    """ Main function."""
    args = parse_args()
    if args.log is None:
        logfile = sys.stdout
    else:
        logfile = args.log
        logging.basicConfig(filename=logfile, level=logging.DEBUG, \
            filemode='w', format='%(asctime)s %(message)s', \
            datefmt='%a, %d %b %Y %H:%M:%S')
    logging.info('Program started.')
    logging.info('Command line: {}'.format(' '.join(sys.argv)))
    with open(args.out, 'w') as vcf_file:
        if args.id_info:
            vcf_reader = vcf.Reader(filename=args.id_info)
        else:
            vcf_reader = None
        write_metadata(args, vcf_file)
        vcf_file.write(OUTPUT_HEADER + '\n')
        blocks = initialise_blocks(args)
        final_blocks = complete_blocks(args, blocks)
        process_blocks(args, final_blocks, vcf_reader, vcf_file)

if __name__ == '__main__':
    main()
