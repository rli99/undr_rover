#!/usr/bin/env python

""" Unmapped primer directed read overlap variant caller. """

from argparse import ArgumentParser
from Bio import (pairwise2, SeqIO)
from itertools import (izip, chain, repeat)
from operator import itemgetter
from pyfaidx import Fasta
import csv
import datetime
import logging
import os
import sys
import vcf

DEFAULT_KMER_THRESHOLD = 0
DEFAULT_PROPORTION_THRESHOLD = 0.05
DEFAULT_ABSOLUTE_THRESHOLD = 2
OUTPUT_HEADER = '\t'.join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", \
    "FILTER", "INFO"])

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
    parser.add_argument('--reference', metavar='FILE', type=str, \
        help='Reference sequences in Fasta format.')
    parser.add_argument('--id_info', type=str, \
    help='File containing rs ID information.')
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

def make_base_sequence(name, bases, qualities):
    """ Take a list of DNA bases and a corresponding list of quality scores
    and return a list of Base objects where the base and score are
    paired together."""
    if len(bases) <= len(qualities):
        return [Base(b, q) for (b, q) in izip(bases, qualities)]
    else:
        logging.warning("In read {}, fewer quality scores ({}) than bases \
            ({}).".format(name, len(qualities), len(bases)))
        # we have fewer quality scores than bases
        # pad the end with 0 scores
        return [Base(b, q) for (b, q) in izip(bases, chain(qualities, \
            repeat(0)))]

def reverse_complement(sequence):
    """ Return the reverse complement of a DNA string."""
    complementary_bases = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}
    rc_bases = []
    for base in sequence:
        rc_bases.append(complementary_bases[str(base)])
    rc_seq = "".join([b for b in rc_bases])
    return rc_seq[::-1]

def nts(none_string):
    """ Turns None into an empty string."""
    if none_string is None:
        return ''
    return str(none_string)

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

class SNV(object):
    """ Single nucleotide variant. Bases are represented as DNA strings."""
    def __init__(self, chrsm, pos, ref_base, seq_base):
        self.chrsm = chrsm
        self.pos = pos
        self.ref_base = ref_base
        self.seq_base = seq_base
        self.qual = '.'
        self.filter_reason = None
        self.info = []
    def __str__(self):
        return "S: {} {} {} {}".format(self.chrsm, self.pos, self.ref_base, \
            self.seq_base)
    def __repr__(self):
        return str(self)
    def as_tuple(self):
        """ Return information about the SNV as a 4-tuple."""
        return (self.chrsm, self.pos, self.ref_base, self.seq_base)
    def __hash__(self):
        return hash(self.as_tuple())
    def __eq__(self, other):
        return self.as_tuple() == other.as_tuple()
    def ref(self):
        """ REF base."""
        return self.ref_base
    def alt(self):
        """ ALT base."""
        return self.seq_base
    def fil(self):
        """ Return "PASS" if the SNV is not filtered, or the reason for
        being discarded otherwise."""
        if self.filter_reason is None:
            return "PASS"
        else:
            return self.filter_reason[1:]
    def position(self):
        """ SNV POS."""
        return self.pos

class Insertion(object):
    """ Insertion. Bases are represented as DNA strings."""
    def __init__(self, chrsm, pos, inserted_bases, context):
        self.chrsm = chrsm
        self.pos = pos
        self.inserted_bases = inserted_bases
        self.qual = '.'
        self.filter_reason = None
        self.info = []
        self.context = context
    def __str__(self):
        return "I: {} {} {}".format(self.chrsm, self.pos, self.inserted_bases)
    def __repr__(self):
        return str(self)
    def as_tuple(self):
        """ Return information about the insertion as a 3-tuple."""
        return (self.chrsm, self.pos, self.inserted_bases)
    def __hash__(self):
        return hash(self.as_tuple())
    def __eq__(self, other):
        return self.as_tuple() == other.as_tuple()
    def ref(self):
        """ REF base."""
        return self.context
    def alt(self):
        """ ALT (inserted) bases."""
        return self.context + self.inserted_bases
    def fil(self):
        """ Return "PASS" if the Insertion is not filtered, or the reason for
        being discarded otherwise."""
        if self.filter_reason is None:
            return "PASS"
        else:
            return self.filter_reason[1:]
    def position(self):
        """ Insertion POS."""
        return self.pos - 1

class Deletion(object):
    """ Deletion. Bases are represented as DNA strings."""
    def __init__(self, chrsm, pos, deleted_bases, context):
        self.chrsm = chrsm
        self.pos = pos
        self.deleted_bases = deleted_bases
        self.qual = '.'
        self.filter_reason = None
        self.info = []
        self.context = context
    def __str__(self):
        return "D: {} {} {}".format(self.chrsm, self.pos, self.deleted_bases)
    def __repr__(self):
        return str(self)
    def as_tuple(self):
        """ Return infromation about the deletion as a 3-tuple."""
        return (self.chrsm, self.pos, self.deleted_bases)
    def __hash__(self):
        return hash(self.as_tuple())
    def __eq__(self, other):
        return self.as_tuple() == other.as_tuple()
    def ref(self):
        """ REF (deleted) bases."""
        return self.context + self.deleted_bases
    def alt(self):
        """ ALT base."""
        return self.context
    def fil(self):
        """ Return "PASS" if the Deletion is not filtered, or the reason for
        being discarded otherwise."""
        if self.filter_reason is None:
            return "PASS"
        else:
            return self.filter_reason[1:]
    def position(self):
        """ Deletion POS."""
        return self.pos - 1

def read_variants(args, chrsm, read, start, insert_seq):
    """Find all the variants in a single read (SNVs, Insertions, Deletions)."""
    result = []
    alignment = pairwise2.align.globalxx(''.join([b.base for b in read]), \
        str(insert_seq), penalize_end_gaps=(False, False), \
        one_alignment_only=True)[0]
    #if alignment[2] < len(insert_seq):
    #    print alignment
    return result

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
        blocks[primer_sequences[block[3]][:20]] = [block[0], block[1], \
        block[2], {}, ref_sequence, block[3], block[4], \
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
            read_bases = str(read.seq)
            primer_key = read_bases[:20]
            if primer_key in blocks:
                if len(blocks[primer_key]) == 9:
                    # Forward primer matched.
                    fseq = blocks[primer_key][7]
                    if fseq == read_bases[:len(fseq)]:
                        if read.id not in blocks[primer_key][3]:
                            blocks[primer_key][3][read.id] = [read, 0, \
                            len(fseq), 0, sample]
                        else:
                            blocks[primer_key][3][read.id][0] = read
                            blocks[primer_key][3][read.id][2] = len(fseq)
                elif len(blocks[primer_key]) == 2:
                    # Reverse primer matched.
                    rseq = blocks[primer_key][0]
                    if rseq == read_bases[:len(rseq)]:
                        forward_key = blocks[primer_key][1]
                        if read.id not in blocks[forward_key][3]:
                            blocks[forward_key][3][read.id] = [0, read, \
                            0, len(rseq), sample]
                        else:
                            blocks[forward_key][3][read.id][1] = read
                            blocks[forward_key][3][read.id][3] = len(rseq)
    # For the next stage, we take only the actual blocks.
    return [b[:5] for b in blocks.values() if len(b) > 2]

def process_blocks(args, blocks, id_info, vcf_file):
    """ Variant calling stage. Process blocks one at a time and call variants
    for each block."""
    coverage_info = []
    for block_info in blocks:
        block_vars = []
        num_pairs = 0
        chrsm, start, end, reads, insert_seq = block_info[:5]
        start = int(start)
        end = int(end)
        logging.info("Processing block chr: {}, start: {}, end: {}"\
            .format(chrsm, start, end))
        for read_pair in reads.values():
            if 0 not in read_pair:
                num_pairs += 1
                read1, read2, fprimerlen, rprimerlen, sample = read_pair
                forward_bases = read1.seq[fprimerlen:]
                reverse_bases = read2.seq[rprimerlen:]

                # make forward base sequence
                read1_bases = make_base_sequence(read1.id, forward_bases, \
                    read1.letter_annotations['phred_quality'])
                # make reverse base sequence
                read2_bases = make_base_sequence(read2.id, reverse_complement(\
                    reverse_bases), read2.letter_annotations['phred_quality'])

                #print ''.join([b.base for b in read1_bases])[:20], \
                #insert_seq[:20]
                #print ''.join([b.base for b in read2_bases])[:-20:-1], \
                #insert_seq[:-20:-1]

                # read variants for the forward read
                #variants1 = read_variants(args, chrsm, read1_bases, start, \
                    #insert_seq)
                # read variants for the reverse read
                #variants2 = read_variants(args, chrsm, read2_bases, start, \
                    #insert_seq)
                #set_variants1, set_variants2 = set(variants1), set(variants2)

                # find the variants each read in the pair share in common
                #same_variants = set_variants1.intersection(set_variants2)

        logging.info("Number of read pairs in block: {}".format(num_pairs))

        for var in block_vars:
            num_vars = block_vars[var]
            proportion = float(num_vars) / num_pairs
            var.info.append("Sample=" + str(sample))
            var.info.append("NV=" + str(num_vars))
            var.info.append("NP=" + str(num_pairs))
            var.info.append("PCT=" + str('{:.2%}'.format(proportion)))
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
        write_coverage_data(args, coverage_file, coverage_info)

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
str(record_info.ID), variant.ref(), variant.alt(), variant.qual, variant\
.fil(), ';'.join(variant.info)]) + '\n')
    else:
        vcf_file.write('\t'.join([variant.chrsm, str(variant.position()), \
'.', variant.ref(), variant.alt(), variant.qual, variant.fil(), ';'\
.join(variant.info)]) + '\n')

def write_coverage_data(args, coverage_file, coverage_info):
    """ Write coverage information to the coverage files."""
    coverage_file.write('chr\tblock_start\tblock_end\tnum_pairs\n')
    for chrsm, start, end, num_pairs in sorted(coverage_info, \
        key=itemgetter(3)):
        coverage_file.write('{}\t{}\t{}\t{}\n'.format(chrsm, start, \
            end, num_pairs))

def write_metadata(args, vcf_file):
    """ Write the opening lines of metadata to the vcf file."""
    vcf_file.write("##fileformat=VCFv4.2" + '\n')
    today = datetime.date.today()
    vcf_file.write("##fileDate=" + str(today)[:4] + str(today)[5:7] + \
        str(today)[8:] + '\n')
    vcf_file.write("##source=ROVER-PCR Variant Caller" + '\n')
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
    if args.absthresh:
        vcf_file.write("##FILTER=<ID=at,Description=\"Variant does not appear \
in at least " + str(args.absthresh) + " read pairs\">" + '\n')
    if args.proportionthresh:
        vcf_file.write("##FILTER=<ID=pt,Descroption=\"Variant does not appear \
in at least " + str(args.proportionthresh*100) \
+ "% of read pairs for the given region\">" + '\n')

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
    logging.info('Command line: {}'.format(' '.join(sys.argv)))
    with open(args.out, 'w') as vcf_file:
        if args.id_info:
            vcf_reader = vcf.Reader(filename=args.id_info)
        write_metadata(args, vcf_file)
        vcf_file.write(OUTPUT_HEADER + '\n')
        blocks = initialise_blocks(args)
        final_blocks = complete_blocks(args, blocks)
        process_blocks(args, final_blocks, vcf_reader, vcf_file)

if __name__ == '__main__':
    main()
