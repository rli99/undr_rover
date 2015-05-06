--------------------------------------------------------------------------------
UNDR ROVER - Unmapped primer directed read overlap variant caller
--------------------------------------------------------------------------------

Version: 1.0.0

Authors: Bernard J Pope (2,3), Tú Nguyen-Dumont (1), Fleur Hammet (1), 
         Daniel J Park (1) and Roger Li (2)

         (1) Genetic Epidemiology Laboratory, Department of Pathology,
             The University of Melbourne.
         (2) Victorian Life Sciences Computation Initiative (VLSCI).
         (3) Department of Computing and Information Systems,
             The University of Melbourne.
         
Web:     https://github.com/rli99/undr_rover

License: <License>

Citation:

    Please cite Undr Rover as follows:

    <Citation>

Requirements: Python 2.7, and the PyVCF, Pyfaidx and Biopython libraries
(https://pypi.python.org/pypi/biopython)
(https://pypi.python.org/pypi/pyfaidx)
(https://pypi.python.org/pypi/PyVCF)

The above requirements can be installed by running:
    pip install -r REQUIREMENTS.txt

--------------------------------------------------------------------------------
General Description
--------------------------------------------------------------------------------

Undr Rover enables the user to call variants from decompressed FASTQ files 
directly without having to go through a mapping step. 

Reads are organised into blocks, which are initialised from information about
known primer sequences. For each block, we record the primer's DNA sequence, 
coordinates and reference insert sequence. Once the reads have been assigned
to blocks (by matching the bases in their primer region to a known primer), 
we call variants for each block one at a time. 

The approach Undr Rover takes when processing reads to call variants is to 
initially assume that all reads contain only single nucleotide variants. If 
we detect two single nucleotide variants in a row during this variant calling, 
Undr Rover immediatly ceases this step and instead does a full gapped alignment 
of the read against the reference insert sequence to detect possible insertions 
and deletions. 

Only variants detected in both reads of a read-pair are considered relevant. 
User-defined thresholds determine the minimum number and proportion of 
read-pairs a variant must be observed to be considered a 'PASS'. Coverage is 
also reported for each block to allow for regions which require additional 
screening to be easily identified.

--------------------------------------------------------------------------------
Command Line Usage
--------------------------------------------------------------------------------

usage: undr_rover [-h] [--version] --primer_coords COORDS --primer_sequences SEQ
                  [--primer_prefix_size N] [--kmer_length N] [--kmer_threshold] 
                  [--primer_bases N] [--proportionthresh N] [--absthresh N] 
                  [--qualthresh N] [--overlap OVERLAP] [--max_variants N] 
                  --reference FILE [--id_info FILE] [--fast] [--snvthresh N]
                  --out FILE [--log FILE] [--coverdir COVERDIR]
                  fastqs [fastqs ...]

positional arguments:
    fastqs                      FASTQ files grouped in pairs. If you have 
                                FASTQ.gz files, they will need to decompressed
                                prior to using with Undr Rover.

optional arguments:
    -h, --help                  Show this help message and exit.
    --version                   Show the program's version number and exit.
    --primer_coords COORDS      Primer coordinates in TSV format.
    --primer_sequences SEQ      Primer base sequences as determined by a primer
                                generating program.
    --primer_prefix_size N 		Size of primer prefix in primer dictionary.
    --kmer_length N             Length of k-mer to use in k-mer test.
    --kmer_threshold N          Number of SNVs in k-mer region deemed 
                                acceptable.
    --primer_bases N            Number of bases from primer region to use in
                                gapped alignment.
    --proportionthresh N        Keep variants which appear in this proportion of
                                the read-pairs for a given target region, and 
                                bin otherwise. Defaults to 0.05.
    --absthresh N               Only keep variants which appear in at least this
                                many read-pairs. Defaults to 2.
    --qualthresh N              Minimum base quality score (phred).
    --overlap OVERLAP           Minimum fraction overlap of read to block
                                region. Defaults to 0.9
    --max_variants N            Maximum amount of variants per read before the
                                read is discarded. Defaults to 25.
    --reference FILE            Reference sequences in FASTA format.
    --id_info FILE              File containing rs ID information in vcf format.
    --fast                  	If this parameter is set, gapped alignment will
                                be used for fewer reads, leading to potentially 
                                less accurate results but much faster run time. 
                                Defaults to False. 
    --snvthresh N               Only applicable if --fast is set. Gapped
                                alignment will be used for reads which have 2 
                                SNVs within N bases of each other. Defaults
                                to 1. 
    --out FILE                  Output file containing called variants.
    --log FILE                  Logs progress in specified file, defaults to 
                                stdout.
    --coverdir COVERDIR         Directory to write coverage files, defaults to
                                current working directory.

Explanation of the arguments
    
    -h

        Print a help message and exit.

    --version

        Print the version number of rover and exit.

    --primer_coords COORDS

        Required.

        A list of primers with their expected coordinates. TSV format with the
        following data:

            chromosome  start   end     forward_primer  reverse_primer

    --primer_sequences SEQ

        Required.

        A list of primers with their expected base sequences. TSV format with
        the following data:

            primer_name sequence

    --primer_prefix_size N

        Optional. Defaults to 20. 

    	Determines the size of the prefix used for the dictionary which stores 
    	primers. Will most likely be set to a value equal to the minimum length 
    	of a primer. 

    --kmer_length N

        Optional. Defaults to 30.

        Determines the length of k-mer to use in k-mer test. In this 
        test, if a certain number of SNVs is found in k-mers from both reads
        then the read pair is discarded.

    --kmer_threshold N

        Optional. Defaults to 2.

        The number of SNVs considered acceptable in the k-mer region of a read.
        If there are more than N SNVs in both k-mers of a read pair, then the 
        read will be discarded. 

    --primer_bases N

        Optional. Defaults to 3.

        When there are repeated bases or repeated sequences of bases near the 
        start of a block, detecting variants becomes slightly harder as some 
        variants won't be apparent. By default, Undr Rover will use 3 additional 
        bases from the end of the primer region during variant calling. This will 
        help to detect those variants. 

        Increasing the value of this parameter will mean more bases from the end 
        of the primer region are used, leading to slightly improved accuracy at 
        the cost of slower run time.  

    --proportionthresh N

        Optional. Defaults to 0.05.

        Only keep variants which appear in this proportion of the read-pairs for
        a given target region, and bin otherwise. A variant must appear in
        both reads of a pair to be counted. The proportion is calculated as
        follows:

            N = number of pairs containing this variant in both reads
            T = number of read-pairs overlapping the target region

            proportion = N/T

        Note: variants must pass BOTH the proportionthresh and absthresh 
        thresholds to be kept. If they fail either test then they are binned.

        That is to say:

            if N/T  >= proportionthresh and N >= absthresh:
                keep the variant
            else:
                bin the variant

    --absthresh N

        Optional. Defaults to 2 read-pairs.

        Only keep variants which appear in at least this many read-pairs
        for a given target region.

        See comments above about proportionthresh.

    --qualthresh N

        Optional. Minimum phred quality score for bases appearing in SNVs and
        insertions. If this argument is set Rover will only consider SNVs and
        insertions where all DNA bases in those variants have a quality score
        greater than or equal to the argument. For example, setting this
        argument to 35 will cause Rover to discard any SNVs or insertions
        containing any bases with a score below 35. If the argument is not
        set then Rover will not consider quality scores in its decision
        to keep or discard a variant.

    --overlap OVERLAP

       Optional. Defaults to 0.9.

       Minimum fraction overlap of read to block region.
       0.5 means at least half of a block must be overlapped by a read
       for the read to be considered for that region. 1.0 would mean the entire
       block must be overlapped by the read.

    --max_variants N

        Optional. Defaults to 25. 

        The maximum number of variants which can be observed for a single read
        before it is considered meaningless and discarded for the purposes of
        variant calling. 

    --reference FILE

        Required.

        Reference sequences in FASTA format. Undr Rover extracts insert 
        sequences for each block from the reference and compares this with the
        read sequences to call variants.

    --id_info FILE

        Optional. Rover will find the rs numbers in DBSNP and show them in the 
        output VCF file if it can find the relevant number in DBSNP, which is a 
        vcf format file containing rs numbers for known variants in .vcf.gz 
        format with an accompanying .tbi file (created by tabix).

    --fast

        Optional. Defaults to False.

        If set, gapped alignment will only be used when SNVs are detected within 
        snvthresh bases of each other. Leads to a significant decrease in the run 
        time of Undr Rover. 

    --snvthresh N

        Optional. Defaults to 1.

        If --fast is set, gapped alignment will be performed on reads
        which contain SNVs within N bases of each other. The default setting
        of 1 allows for the fastest run-time and should provide a sufficient 
        level of accuracy in the majority of cases.
     
    --coverdir COVERDIR

        Optional. Defaults to current working directory.

        Directory to write the coverage files. If the directory does not
        exist Rover will not try to create it. You must create the directory
        yourself.

        The format of the coverage files is a TSV with:

        chr     block_start     block_end       num_pairs

    fastqs [fastqs ...]

        One or more pairs of FASTQ files containing reads for which variant
        calling will be attempted. FASTQ.gz files will need to decompressed
        prior to use with Undr Rover.

--------------------------------------------------------------------------------
Example usage (should all be on one line)
--------------------------------------------------------------------------------

undr_rover.py   --primer_coords roverfile.txt --primer_sequences primers.txt 
                --id_info dbsnp.vcf.gz --log log_file --fast
                --reference fasta/hg19bis.fa --out variant_calls.vcf 
                --coverdir coverage_files sample1*.fastq sample2*.fastq

This assumes that there are two FASTQ files each for sample1* and sample2*, 
which contain the first and second reads respectively for each read-pair. 
Roverfile.txt and primers.txt are both TSV format files containing data as 
described earlier. The --fast flag is also set to decrease the run time
of Undr Rover. 

The variants called will be written to the file variant_calls.vcf. Coverage 
files containing the number of read-pairs which mapped to each region will be 
output in coverage_files/sample1.coverage and coverage_files/sample2.coverage.
A log file describing the actions taken by the program will be stored in 
log_file.
