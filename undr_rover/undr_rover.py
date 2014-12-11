#!/usr/bin/env python

""" Unmapped primer directed read overlap variant caller. """

from argparse import ArgumentParser
from Bio import pairwise2
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from itertools import takewhile
from operator import itemgetter
from pyfaidx import Fasta
from version import undr_rover_version
import csv
import datetime
import logging
import os
import string
import sys
import vcf

DEFAULT_ABSOLUTE_THRESHOLD = 2
DEFAULT_KMER_LENGTH = 30
DEFAULT_MAX_VARIANTS = 20
DEFAULT_MINIMUM_READ_OVERLAP_BLOCK = 0.9
DEFAULT_PRIMER_BASES = 5
DEFAULT_PROPORTION_THRESHOLD = 0.05
DEFAULT_SNV_THRESHOLD = 1
OUTPUT_HEADER = '\t'.join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", \
    "FILTER", "INFO"])

def parse_args():
    """ Find variants from fastqs via a mapping-free approach."""
    parser = ArgumentParser(description="Find variants from fastqs via a \
        mapping-free approach.")
    parser.add_argument('--version', action='version', version='%(prog)s ' \
        + undr_rover_version)
    parser.add_argument('--primer_coords', type=str, required=True, \
        help='Primer coordinates in TSV format.')
    parser.add_argument('--primer_sequences', metavar='FILE', type=str, \
        required=True, help='Primer base sequences as determined by a primer \
        generating program.')
    parser.add_argument('--kmer_length', type=int, \
        default=DEFAULT_KMER_LENGTH, \
        help='Length of k-mer to check after the primer sequence.' \
        'Defaults to {}.'.format(DEFAULT_KMER_LENGTH))
    parser.add_argument('--kmer_threshold', type=int, \
        help='Number of single nucleotide variants deemed acceptable in kmer.')
    parser.add_argument('--primer_bases', type=int, \
        default=DEFAULT_PRIMER_BASES, \
        help='Number of bases from primer region to use in gapped alignment.' \
        'Helps with variant calling near the edges of a block.' \
        'Defaults to {}.'.format(DEFAULT_PRIMER_BASES))
    parser.add_argument('--proportionthresh', metavar='N', type=float, \
        default=DEFAULT_PROPORTION_THRESHOLD, \
        help='Keep variants which appear in this proportion of the read \
        pairs. For a given target region, and bin otherwise.' \
        'Defaults to {}.'.format(DEFAULT_PROPORTION_THRESHOLD))
    parser.add_argument('--absthresh', metavar='N', type=int, \
        default=DEFAULT_ABSOLUTE_THRESHOLD, \
        help='Only keep variants which appear in at least this many \
        read pairs. ' \
        'Defaults to {}.'.format(DEFAULT_ABSOLUTE_THRESHOLD))
    parser.add_argument('--qualthresh', metavar='N', type=int, \
        help='Minimum base quality score (phred).')
    parser.add_argument('--overlap', type=float, \
        default=DEFAULT_MINIMUM_READ_OVERLAP_BLOCK, \
        help='Minimum proportion of block which must be overlapped by a read.' \
        'Defaults to {}.'.format(DEFAULT_MINIMUM_READ_OVERLAP_BLOCK))
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
    parser.add_argument('--coverdir', \
        help='Directory to write coverage files, defaults to current working \
        directory.')
    parser.add_argument('--thorough', action='store_true', default=False, \
        help='Use gapped alignment more often.')
    parser.add_argument('--snvthresh', metavar='N', type=int, \
        default=DEFAULT_SNV_THRESHOLD, \
        help='Distance between two single nucleotide variants before going to \
        a gapped alignment.')
    parser.add_argument('fastqs', nargs='+', type=str, \
        help='Fastq files containing reads.')
    return parser.parse_args()

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

class Variant(object):
    """ Parent class representing all variants."""
    def __init__(self, chrsm, pos, _, qual):
        self.chrsm = chrsm
        self.pos = pos
        self.qual = qual
        self.filter_reason = None
        self.info = []
    def __eq__(self, other):
        return self.as_tuple() == other.as_tuple()
    def __hash__(self):
        return hash(self.as_tuple())
    def __repr__(self):
        return str(self)
    def fil(self):
        """ Return "PASS" if the Deletion is not filtered, or the reason(s) for
        being discarded otherwise."""
        if self.filter_reason is None:
            return "PASS"
        return self.filter_reason[1:]

class SNV(Variant):
    """ Single nucleotide variant. Bases are represented as DNA strings."""
    def __init__(self, chrsm, pos, bases, qual):
        super(SNV, self).__init__(chrsm, pos, bases, qual)
        self.ref_base = bases[0]
        self.seq_base = bases[1]
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

class Insertion(Variant):
    """ Insertion. Bases are represented as DNA strings."""
    def __init__(self, chrsm, pos, bases, qual):
        super(Insertion, self).__init__(chrsm, pos, bases, qual)
        self.inserted_bases = bases[0]
        self.context = bases[1]
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

class Deletion(Variant):
    """ Deletion. Bases are represented as DNA strings."""
    def __init__(self, chrsm, pos, bases, qual):
        super(Deletion, self).__init__(chrsm, pos, bases, qual)
        self.deleted_bases = bases[0]
        self.context = bases[1]
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

def ascii_to_phred(ascii):
    """ Quality score of a base is stored as a byte (ascii character) in
    "Qual plus 33 format". So we subtract off 33 from the ascii code to get
    the actual score."""
    return ord(ascii) - 33

def nts(none_string):
    """ Returns an empty string for None."""
    return none_string or ''

COMPLEMENT = string.maketrans('ATCGN', 'TAGCN')
def reverse_complement(sequence):
    """ Return the reverse complement of a DNA string."""
    return sequence.translate(COMPLEMENT)[::-1]

def read_snvs(args, chrsm, qual, pos, insert_seq, bases, direction):
    """ Find all the SNVs in a read. 0 if we find two SNVs within a certain
    distance from each other."""
    # If --thorough is set, any reads in which 2 SNVs are detected will undergo
    # a complete gapped alignment instead.
    min_distance = len(insert_seq) if args.thorough else args.snvthresh
    check = min_distance
    pos -= args.primer_bases * direction
    result = []
    # Initiate the relevant indices given the direction of the read.
    i = 0 if direction == 1 else -1
    new = slice(0, 1) if direction == 1 else slice(-1, None)
    insert = insert_seq[args.primer_bases:-1 * args.primer_bases]

    # Check for identical insert sequence and read (no variants).
    if direction == 1 and insert == bases[args.primer_bases:len(insert_seq) - \
    args.primer_bases]:
        return (result, 0)
    if direction == -1 and insert == bases[-1 * len(insert_seq) + \
    args.primer_bases:-1 * args.primer_bases]:
        return (result, 0)

    kmer_stop = pos + args.kmer_length * direction
    count = 0
    insert_seq, bases = list(insert_seq), list(bases)
    while insert_seq and bases:
        if insert_seq[i] == bases[i]:
            del insert_seq[new], bases[new], qual[new]
            check += 1
            pos += direction
        else:
            # Increment count if we see an SNV within the k-mer region.
            if pos < kmer_stop and direction == 1 or pos > kmer_stop and \
            direction == -1:
                count += 1
            if check <= min_distance:
                if not args.kmer_threshold:
                    return (0, 0)
                # If SNVs within a certain distance have been detected, and
                # we're past the k-mer region, return 0 as well as the count
                # of SNVs in the k-mer region.
                if pos >= kmer_stop and direction == 1 or pos <= kmer_stop and \
                direction == -1:
                    return (0, count)
            result.append(SNV(chrsm, pos, [insert_seq[i], bases[i]], '.'))
            if args.qualthresh and ascii_to_phred(qual[i]) < args.qualthresh:
                result[-1].filter_reason = ''.join([nts(result[-1].\
                    filter_reason), ";qlt"])
            del insert_seq[new], bases[new], qual[new]
            pos += direction
            check = 1
    # Return both the detected SNVs, and the count of SNVs within the k-mer
    # region, which will be used to determine which reads pass the k-mer test.
    return (result, count)

def read_variants(args, chrsm, qual, pos, insert_seq, bases, direction):
    """ Find all the variants in a read (SNVs, Insertions, Deletions)."""
    pos -= args.primer_bases * direction
    result = []
    # Pairwise2 uses a C implementation of the Needleman-Wunsch algorithm.
    aligned_insert, aligned_read = pairwise2.align.globalms(insert_seq, bases, \
        2, 0, -2, -1, penalize_end_gaps=(0, 0), one_alignment_only=1)[0][:2]
    context = '-'
    aligned_insert, aligned_read = list(aligned_insert), list(aligned_read)
    while aligned_insert and aligned_read:
        # Initialise the relevant indices given the direction of the read.
        i = 0 if direction == 1 else -1
        new = slice(0, 1) if direction == 1 else slice(-1, None)
        if aligned_insert[i] == aligned_read[i]:
            # Match
            context = aligned_insert[i]
            del aligned_insert[new], aligned_read[new], qual[new]
            pos += direction
        elif not (aligned_insert[i] == '-' or aligned_read[i] == '-'):
            # Single Nucleotide Variation
            result.append(SNV(chrsm, pos, [aligned_insert[i], \
                aligned_read[i]], '.'))
            if args.qualthresh and ascii_to_phred(qual[i]) < args.qualthresh:
                result[-1].filter_reason = ''.join([nts(result[-1].\
                    filter_reason), ";qlt"])
            context = aligned_insert[i]
            del aligned_insert[new], aligned_read[new], qual[new]
            pos += direction
        elif aligned_insert[i] == '-':
            # Insertion
            indel_length = sum(1 for _ in takewhile(lambda x: x == '-', \
                aligned_insert[::direction]))
            new = slice(0, indel_length) if direction == 1 else \
            slice(-1 * indel_length, None)
            if indel_length >= len(aligned_insert):
                return result
            if direction == 1:
                result.append(Insertion(chrsm, pos, \
                [''.join(aligned_read[:indel_length]), context], '.'))
            else:
                result.append(Insertion(chrsm, pos + 1, [''.join\
                    (aligned_read[-1 * indel_length:]), \
                    ''.join(aligned_insert[-1 * indel_length - 1])], '.'))
            # insertion with QUAL data?
            if args.qualthresh and any([ascii_to_phred(b) < args.qualthresh \
                for b in qual[new]]):
                result[-1].filter_reason = ''.join([nts(result[-1].\
                filter_reason), ";qlt"])
            del aligned_insert[new], aligned_read[new], qual[new]
        elif aligned_read[i] == '-':
            # Deletion
            indel_length = sum(1 for _ in takewhile(lambda x: x == '-', \
                aligned_read[::direction]))
            new = slice(0, indel_length) if direction == 1 else \
            slice(-1 * indel_length, None)
            if indel_length >= len(aligned_insert):
                return result
            if direction == 1:
                result.append(Deletion(chrsm, pos, \
                [''.join(aligned_insert[:indel_length]), context], '.'))
            else:
                result.append(Deletion(chrsm, pos - indel_length + 1, \
                    [''.join(aligned_insert[-1 * indel_length:]), \
                    ''.join(aligned_insert[-1 * indel_length - 1])], '.'))
            # deletion with QUAL data?
            if args.qualthresh and ascii_to_phred(qual[i]) < args.qualthresh:
                result[-1].filter_reason = ''.join([nts(result[-1].\
                filter_reason), ";qlt"])
            context = aligned_insert[indel_length - 1]
            del aligned_insert[new], aligned_read[new]
            pos += indel_length * direction
    return result

def initialise_blocks(args):
    """ Create blocks, initially containing block coordinates, primer sequences
    and an empty dictionary into which reads pairs will be added. The blocks
    themselves are a dictionary at this stage."""
    blocks = {}
    primer_sequences = {}
    block_coords = list(csv.reader(open(args.primer_coords), delimiter='\t'))
    primer_info = list(csv.reader(open(args.primer_sequences), delimiter='\t'))
    reference = Fasta(args.reference)
    for primer in primer_info:
        # Initiates a dictionary containing containing primer sequences.
        primer_sequences[primer[0]] = primer[1]
    for block in block_coords:
        # Get the insert sequence. Take some extra bases to help with the
        # variant calling near the edges of the block.
        ref_sequence = reference[block[0]][int(block[1]) - 1 - \
        args.primer_bases:int(block[2]) + args.primer_bases]
        # Actual block, for which the key is the first 20 bases of the forward
        # primer. Value contains [chr, start, end, {reads}, insert_seq, forward
        # primer name, reverse primer name, forward primer sequence, reverse
        # primer sequence]
        blocks[primer_sequences[block[3]][:20]] = [block[0], block[1], \
        block[2], {}, str(ref_sequence), block[3], block[4], \
        primer_sequences[block[3]], primer_sequences[block[4]]]
        # Reverse primer (not an actual block), contains the key for the forward
        # primer. Redirects to forward primer block.
        blocks[primer_sequences[block[4]][:20]] = [primer_sequences[block[4]], \
        primer_sequences[block[3]][:20]]
    return blocks

def complete_blocks(args, blocks, fastq_pair):
    """ Organise reads into blocks."""
    for block in blocks:
        if len(blocks[block]) > 2:
            blocks[block][3].clear()
    sample = os.path.basename(fastq_pair[0]).split('_')
    if len(sample) > 0:
        sample = '_'.join(sample[:3])
        logging.info("Processing sample {}".format(sample))
    else:
        exit('Cannot deduce sample name from fastq filename {}'.\
            format(fastq_pair[0]))
    for fastq_file in fastq_pair:
        with open(fastq_file, "rU") as fastq:
            for (title, seq, qual) in FastqGeneralIterator(fastq):
                # Each read is also stored as a dictionary.
                read = {'name': title.partition(' ')[0], 'seq': seq}
                read['qual'] = qual if args.qualthresh else []
                # Try to match each read (check the first 20 bases) with an
                # expected primer.
                read_bases = read['seq']
                primer_key = read_bases[:20]
                if len(blocks.get(primer_key, [])) == 9:
                    # Possible forward primer matched.
                    fseq = blocks[primer_key][7]
                    if fseq == read_bases[:len(fseq)]:
                        if read['name'] not in blocks[primer_key][3]:
                            blocks[primer_key][3][read['name']] = [read, 0, \
                            len(fseq), 0, sample]
                        else:
                            blocks[primer_key][3][read['name']][0] = read
                            blocks[primer_key][3][read['name']][2] = len(fseq)
                elif len(blocks.get(primer_key, [])) == 2:
                    # Possible reverse primer matched.
                    rseq = blocks[primer_key][0]
                    if rseq == read_bases[:len(rseq)]:
                        forward_key = blocks[primer_key][1]
                        if read['name'] not in blocks[forward_key][3]:
                            blocks[forward_key][3][read['name']] = [0, read, \
                            0, len(rseq), sample]
                        else:
                            blocks[forward_key][3][read['name']][1] = read
                            blocks[forward_key][3][read['name']][3] = len(rseq)
    # For the next stage, we only need the actual blocks.
    return [b[:5] for b in blocks.values() if len(b) > 2]

def process_blocks(args, blocks, id_info, vcf_file):
    """ Variant calling stage. Process blocks one at a time and call variants
    for each block."""
    coverage_info = []
    for block_info in sorted(blocks, key=itemgetter(0, int(1))):
        block_vars = {}
        num_pairs, total_pairs, kmer_fail = 0, 0, 0
        chrsm, start, end, reads, insert_seq = block_info[:5]
        start = int(start)
        end = int(end)
        logging.info("Processing block chr: {}, start: {}, end: {}"\
            .format(chrsm, start, end))
        for read_pair in [r for r in reads.values() if 0 not in r]:
            num_pairs += 1
            total_pairs += 1
            read1, read2, fprimerlen, rprimerlen, sample = read_pair
            variants1, variants2 = [], []

            forward_bases = read1['seq'][fprimerlen - args.primer_bases:]
            reverse_bases = read2['seq'][rprimerlen - args.primer_bases:]
            insert = insert_seq.upper()
            min_overlap = args.overlap * len(insert_seq)
            if len(forward_bases) - args.primer_bases > min_overlap and \
            len(reverse_bases) - args.primer_bases > min_overlap:
                # For both reads, we initially assume that they only have single
                # nucleotide variants. If we detect results which may suggest
                # otherwise, we go to a gapped alignment.
                variants1, kmer1 = read_snvs(args, chrsm, read1['qual'], \
                    start, insert, forward_bases, 1)
                variants2, kmer2 = read_snvs(args, chrsm, read2['qual'], end, \
                    insert, reverse_complement(reverse_bases), -1)
                # K-mer test. If either read fails, disregard the read pair.
                if args.kmer_threshold:
                    if kmer1 > args.kmer_threshold and kmer2 > \
                    args.kmer_threshold:
                        variants1, variants2 = [], []
                        num_pairs -= 1
                        kmer_fail += 1
            else:
                num_pairs -= 1

            # Gapped alignment for the reads in which we have detected the
            # possibility of indels.
            if variants1 == 0:
                variants1 = read_variants(args, chrsm, read1['qual'], start, \
                    insert, forward_bases, 1)
            if variants2 == 0:
                variants2 = read_variants(args, chrsm, read2['qual'], end, \
                    insert, reverse_complement(reverse_bases), -1)

            # Ignore reads which have an unusually high amount of variants.
            if len(variants1) > args.max_variants or len(variants2) > \
            args.max_variants:
                num_pairs -= 1
                variants1, variants2 = [], []

            # Consider variants each read in the pair share in common.
            for var in set(variants1).intersection(set(variants2)):
                # Only consider variants within the bounds of the block.
                if var.pos >= start and var.pos <= end:
                    block_vars[var] = block_vars.get(var, 0) + 1

        logging.info("Number of read pairs in block: {}".format(total_pairs))
        if kmer_fail > 0:
            logging.info("Number of read pairs which failed k-mer test: {}"\
                .format(kmer_fail))
        logging.info("Number of acceptable read pairs in block: {}"\
            .format(num_pairs))
        logging.info("Number of variants found in block: {}".\
            format(len(block_vars)))

        for var in block_vars:
            num_vars = block_vars[var]
            proportion = float(num_vars) / num_pairs
            var.info.extend([''.join(["Sample=", sample]), ''.join(["NV=", \
                str(num_vars)]), ''.join(["NP=", str(num_pairs)]), \
            ''.join(["PCT=", str('{:.2%}'.format(proportion))])])
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
    # If the variant is deemed a "PASS", find the relevant rs number from dbsnp.
    ref = '.'  # Default
    if variant.fil() == "PASS" and args.id_info:
        for record in id_info.fetch(variant.chrsm, variant.position(), variant.\
        position() + max(len(variant.ref()), len(variant.alt())) + 1):
            if record.POS == variant.position() and record.REF == \
            variant.ref() and (variant.alt() in record.ALT):
                ref = str(record.ID)
    vcf_file.write('\t'.join([variant.chrsm, str(variant.position()), \
ref, variant.ref(), variant.alt(), str(variant.qual), variant.fil(), \
';'.join(variant.info)]) + '\n')

def write_coverage_data(coverage_file, coverage_info):
    """ Write coverage information to the coverage files."""
    coverage_file.write('chr\tblock_start\tblock_end\tnum_pairs\n')
    for chrsm, start, end, num_pairs in sorted(coverage_info, \
    key=itemgetter(3)):
        coverage_file.write('{}\t{}\t{}\t{}\n'.format(chrsm, start, end, \
            num_pairs))

def write_metadata(args, vcf_file):
    """ Write the opening lines of metadata to the vcf file."""
    vcf_file.write("##fileformat=VCFv4.2\n")
    today = str(datetime.date.today())
    vcf_file.write("##fileDate={}{}{}\n".format(today[:4], today[5:7], \
        today[8:]))
    vcf_file.write("##source=UNDR ROVER\n")
    vcf_file.write("##INFO=<ID=Sample,Number=1,Type=String,Description=\
\"Sample Name\">\n")
    vcf_file.write("##INFO=<ID=NV,Number=1,Type=Float,Description=\"Number of \
read pairs with variant\">\n")
    vcf_file.write("##INFO=<ID=NP,Number=1,Type=Float,Description=\"Number of \
read pairs at POS\">\n")
    vcf_file.write("##INFO=<ID=PCT,Number=1,Type=Float,Description=\
\"Percentage of read pairs at POS with variant\">\n")
    if args.qualthresh:
        vcf_file.write("##FILTER=<ID=qlt,Description=\"Variant has phred \
quality score below {}\">\n".format(args.qualthresh))
    vcf_file.write("##FILTER=<ID=at,Description=\"Variant does not appear in \
at least {} read pairs\">\n".format(args.absthresh))
    vcf_file.write("##FILTER=<ID=pt,Description=\"Variant does not appear in \
at least {}% of read pairs for the given region\">\n"\
.format(args.proportionthresh * 100))

def main():
    """ Main function."""
    args = parse_args()
    if args.log:
        logfile = args.log
        logging.basicConfig(filename=logfile, level=logging.DEBUG, \
            filemode='w', format='%(asctime)s %(message)s', \
            datefmt='%a, %d %b %Y %H:%M:%S')
    else:
        logfile = sys.stdout
    logging.info('Program started.')
    logging.info('Command line: {}'.format(' '.join(sys.argv)))
    with open(args.out, 'w') as vcf_file:
        vcf_reader = vcf.Reader(filename=args.id_info) if args.id_info else None
        write_metadata(args, vcf_file)
        vcf_file.write(OUTPUT_HEADER + '\n')
        blocks = initialise_blocks(args)
        for fastq_pair in zip(*[iter(args.fastqs)]*2):
            final_blocks = complete_blocks(args, blocks, fastq_pair)
            process_blocks(args, final_blocks, vcf_reader, vcf_file)
    # from guppy import hpy
    # h = hpy()
    # print h.heap()

if __name__ == '__main__':
    main()
