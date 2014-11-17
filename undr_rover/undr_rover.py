#!/usr/bin/env python

""" Unmapped primer directed read overlap variant caller. """

from argparse import ArgumentParser
from Bio import (pairwise2, SeqIO)
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
        mapping-free approach.""")
    parser.add_argument('--primer_coords', type=str, required=True, \
        help='Primer coordinates in TSV format.')
    parser.add_argument('--primer_sequences', metavar='FILE', type=str, \
        help='Primer base sequences as determined by a primer generating \
        program.')
    parser.add_argument('--primer_bases', type=int, \
        default=DEFAULT_PRIMER_BASES, \
        help='Number of bases from primer region to use in gapped alignment.')
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
    parser.add_argument('--max_variants', metavar='N', type=int, \
        default=DEFAULT_MAX_VARIANTS, \
        help='Ignore reads with greater than this many variants observed.' \
        'Defaults to {}.'.format(DEFAULT_MAX_VARIANTS))
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

def reverse_complement(sequence):
    """ Return the reverse complement of a DNA string."""
    rc_bases = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
    return "".join([rc_bases[base] for base in sequence[::-1]])

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
    def __init__(self, chrsm, pos, ref_base, seq_base, qual):
        self.chrsm = chrsm
        self.pos = pos
        self.ref_base = ref_base
        self.seq_base = seq_base
        self.qual = qual
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
    def __init__(self, chrsm, pos, inserted_bases, context, qual):
        self.chrsm = chrsm
        self.pos = pos
        self.inserted_bases = inserted_bases
        self.qual = qual
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
    def __init__(self, chrsm, pos, deleted_bases, context, qual):
        self.chrsm = chrsm
        self.pos = pos
        self.deleted_bases = deleted_bases
        self.qual = qual
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

def read_fsnvs(args, chrsm, qual, pos, insert_seq, bases):
    """ Find all SNV's in a forward read. 0 if we find two SNV's in a row."""
    check = False
    pos -= args.primer_bases
    result = []
    # Identical insert sequence and read, so there are no variants.
    if insert_seq[args.primer_bases:-1 * args.primer_bases] == \
    bases[args.primer_bases:len(insert_seq) - args.primer_bases]:
        return result
    while insert_seq and bases:
        if insert_seq[0] == bases[0]:
            insert_seq = insert_seq[1:]
            bases = bases[1:]
            qual = qual[1:]
            check = False
            pos += 1
        else:
            # If we just saw an SNV, return 0 as there are now two in a row.
            if check == True:
                return 0
            else:
                check = True
            snv = SNV(chrsm, pos, insert_seq[0], bases[0], '.')
            if not ((args.qualthresh is None) or (qual[0] >= \
                args.qualthresh)):
                snv.filter_reason = ''.join([nts(snv.filter_reason), \
                ";qlt"])
            result.append(snv)
            insert_seq = insert_seq[1:]
            bases = bases[1:]
            qual = qual[1:]
            pos += 1
    return result

def read_rsnvs(args, chrsm, qual, pos, insert_seq, bases):
    """ Find all SNV's in a reverse read. 0 if we find two SNV's in a row."""
    check = False
    pos += args.primer_bases
    result = []
    # Identical insert sequence and read, so there are no variants.
    if insert_seq[args.primer_bases:-1 * args.primer_bases] == \
    bases[-1 * len(insert_seq) + args.primer_bases:-1 * args.primer_bases]:
        return result
    while insert_seq and bases:
        if insert_seq[-1] == bases[-1]:
            insert_seq = insert_seq[:-1]
            bases = bases[:-1]
            qual = qual[:-1]
            check = False
            pos -= 1
        else:
            # If we just saw an SNV, return 0 as there are now two in a row.
            if check == True:
                return 0
            else:
                check = True
            snv = SNV(chrsm, pos, insert_seq[-1], bases[-1], '.')
            if not ((args.qualthresh is None) or (qual[-1] >= \
                args.qualthresh)):
                snv.filter_reason = ''.join([nts(snv.filter_reason), \
                ";qlt"])
            result.append(snv)
            insert_seq = insert_seq[:-1]
            bases = bases[:-1]
            qual = qual[:-1]
            pos -= 1
    return result

def read_fvariants(args, chrsm, qual, pos, insert_seq, bases):
    """ Find all the variants in a forward read (SNVs, Insertions,
    Deletions)."""
    pos -= args.primer_bases
    result = []
    # Identical insert sequence and read, so there are no variants.
    if insert_seq[args.primer_bases:-1 * args.primer_bases] == \
    bases[args.primer_bases:len(insert_seq) - args.primer_bases]:
        return result
    # Pairwise2 uses a C implementation of the Needleman-Wunsch algorithm.
    alignment = pairwise2.align.globalms(insert_seq, bases, 2, 0, -2, -1, \
        penalize_end_gaps=(0, 0), one_alignment_only=1)[0]
    aligned_insert, aligned_read = alignment[0], alignment[1]
    context = '-'
    while aligned_insert and aligned_read:
        if aligned_insert[0] == aligned_read[0]:
            context = aligned_insert[0]
            aligned_insert = aligned_insert[1:]
            aligned_read = aligned_read[1:]
            qual = qual[1:]
            pos += 1
        else:
            if not (aligned_insert[0] == '-' or aligned_read[0] == '-'):
                # Single Nucleotide Variation
                snv = SNV(chrsm, pos, aligned_insert[0], aligned_read[0], '.')
                # snv = SNV(chrsm, pos, aligned_insert[0], aligned_read[0], \
                #     qual[0])
                if not ((args.qualthresh is None) or (qual[0] >= \
                    args.qualthresh)):
                    snv.filter_reason = ''.join([nts(snv.filter_reason), \
                    ";qlt"])
                result.append(snv)
                context = aligned_insert[0]
                aligned_insert = aligned_insert[1:]
                aligned_read = aligned_read[1:]
                qual = qual[1:]
                pos += 1
            elif aligned_insert[0] == '-':
                insert = True # Insertion
                ins_length = 0
                for base in aligned_insert:
                    if base == '-' and insert:
                        ins_length += 1
                    else:
                        insert = False
                if ins_length >= len(aligned_insert):
                    return result
                insertion = Insertion(chrsm, pos, aligned_read\
                    [:ins_length], context, '.')
                # insertion with QUAL data?
                if not ((args.qualthresh is None) or all([b >= \
                    args.qualthresh for b in qual[ins_length:]])):
                    insertion.filter_reason = ''.join([nts(insertion.\
                    filter_reason), ";qlt"])
                result.append(insertion)
                aligned_read = aligned_read[ins_length:]
                aligned_insert = aligned_insert[ins_length:]
                qual = qual[ins_length:]
            elif aligned_read[0] == '-':
                delete = True # Deletion
                del_length = 0
                for base in aligned_read:
                    if base == '-' and delete:
                        del_length += 1
                    else:
                        delete = False
                if del_length >= len(aligned_insert):
                    return result
                deletion = Deletion(chrsm, pos, aligned_insert\
                    [:del_length], context, '.')
                # deletion with QUAL data?
                if not ((args.qualthresh is None) or (qual[0] >= \
                    args.qualthresh)):
                    deletion.filter_reason = ''.join([nts(deletion.\
                    filter_reason), ";qlt"])
                result.append(deletion)
                context = aligned_insert[del_length - 1]
                aligned_insert = aligned_insert[del_length:]
                aligned_read = aligned_read[del_length:]
                pos += del_length
    return result

def read_rvariants(args, chrsm, qual, pos, insert_seq, bases):
    """ Find all the variants in a reverse read (SNVs, Insertions,
        Deletions)."""
    pos += args.primer_bases
    result = []
    # Identical insert sequence and read, so there are no variants.
    if insert_seq[args.primer_bases:-1 * args.primer_bases] == \
    bases[-1 * len(insert_seq) + args.primer_bases:-1 * args.primer_bases]:
        return result
    # Pairwise2 uses a C implementation of the Needleman-Wunsch algorithm.
    alignment = pairwise2.align.globalms(insert_seq, bases, 2, 0, -2, -1, \
        penalize_end_gaps=(0, 0), one_alignment_only=1)[0]
    aligned_insert, aligned_read = alignment[0], alignment[1]
    while aligned_insert and aligned_read:
        if aligned_insert[-1] == aligned_read[-1]:
            aligned_insert = aligned_insert[:-1]
            aligned_read = aligned_read[:-1]
            qual = qual[:-1]
            pos -= 1
        else:
            if not (aligned_insert[-1] == '-' or aligned_read[-1] == '-'):
                # Single Nucleotide Variation
                snv = SNV(chrsm, pos, aligned_insert[-1], aligned_read[-1], '.')
                # snv = SNV(chrsm, pos, aligned_insert[-1], aligned_read[-1], \
                #     qual[0])
                if not ((args.qualthresh is None) or (qual[-1] >= \
                    args.qualthresh)):
                    snv.filter_reason = ''.join([nts(snv.filter_reason), \
                    ";qlt"])
                result.append(snv)
                aligned_insert = aligned_insert[:-1]
                aligned_read = aligned_read[:-1]
                qual = qual[:-1]
                pos -= 1
            elif aligned_insert[-1] == '-':
                insert = True # Insertion
                ins_length = 0
                for base in aligned_insert[::-1]:
                    if base == '-' and insert:
                        ins_length += 1
                    else:
                        insert = False
                if ins_length >= len(aligned_insert):
                    return result
                insertion = Insertion(chrsm, pos + 1, aligned_read[-1 * \
                    ins_length:], aligned_insert[-1 * ins_length - 1], '.')
                # insertion with QUAL data?
                if not ((args.qualthresh is None) or all([b >= args.qualthresh \
                    for b in qual[:-1 * ins_length]])):
                    insertion.filter_reason = ''.join([nts(insertion.\
                    filter_reason), ";qlt"])
                result.append(insertion)
                aligned_read = aligned_read[:-1 * ins_length]
                aligned_insert = aligned_insert[:-1 * ins_length]
                qual = qual[:-1 * ins_length]
            elif aligned_read[-1] == '-':
                delete = True # Deletion
                del_length = 0
                for base in aligned_read[::-1]:
                    if base == '-' and delete:
                        del_length += 1
                    else:
                        delete = False
                if del_length >= len(aligned_insert):
                    return result
                deletion = Deletion(chrsm, pos - del_length + 1, \
                    aligned_insert[-1 * del_length:], aligned_insert\
                    [-1 * del_length - 1], '.')
                # deletion with QUAL data?
                if not ((args.qualthresh is None) or (qual[-1] >= \
                    args.qualthresh)):
                    deletion.filter_reason = ''.join([nts(deletion.\
                    filter_reason), ";qlt"])
                result.append(deletion)
                aligned_insert = aligned_insert[:-1 * del_length]
                aligned_read = aligned_read[:-1 * del_length]
                pos -= del_length
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
        for read_pair in reads.values():
            if 0 not in read_pair:
                num_pairs += 1
                read1, read2, fprimerlen, rprimerlen, sample = read_pair
                forward_bases = read1.seq[fprimerlen - args.primer_bases:]
                reverse_bases = read2.seq[rprimerlen - args.primer_bases:]

                forward_qual = read1.letter_annotations["phred_quality"]\
                [fprimerlen - args.primer_bases:]
                reverse_qual = read2.letter_annotations["phred_quality"]\
                [rprimerlen - args.primer_bases:]

                insert = insert_seq.upper()
                forward_seq = str(forward_bases.upper())
                reverse_seq = reverse_complement(str(reverse_bases)).upper()

                # For both reads, we initially assume that they only have single
                # nucleotide variants. If we detect successive SNV's in those
                # reads, we go back and do the gapped alignment.

                variants1 = read_fsnvs(args, chrsm, forward_qual, start, \
                    insert, forward_seq)
                if variants1 == 0:
                    variants1 = read_fvariants(args, chrsm, forward_qual, \
                        start, insert, forward_seq)

                variants2 = read_rsnvs(args, chrsm, reverse_qual, end, \
                    insert, reverse_seq)
                if variants2 == 0:
                    variants2 = read_rvariants(args, chrsm, reverse_qual, end, \
                        insert, reverse_seq)

                # Ignore reads which have an unusually high amount of variants.
                if len(variants1) > args.max_variants or len(variants2) > \
                args.max_variants:
                    variants1, variants2 = [], []
                    logging.info("Read {} discarded due to an unusually high \
amount of variants.".format(read1.id))

                # Find the variants each read in the pair share in common.
                set_variants1, set_variants2 = set(variants1), set(variants2)
                same_variants = set_variants1.intersection(set_variants2)

                for var in same_variants:
                    # Only consider variants within the bounds of the block.
                    if var.pos >= start and var.pos <= end:
                        if var in block_vars:
                            block_vars[var] += 1
                        else:
                            block_vars[var] = 1

        logging.info("Number of read pairs in block: {}".format(num_pairs))
        logging.info("Number of variants found in block: {}".\
            format(len(block_vars)))

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
    today = datetime.date.today()
    vcf_file.write("##fileDate=" + str(today)[:4] + str(today)[5:7] + \
        str(today)[8:] + '\n')
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
