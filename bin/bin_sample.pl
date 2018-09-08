#! /usr/bin/env perl

# This script is the executable for the each sample
#
# Scott Yourstone
# May 8, 2012

use strict;
use warnings;

use Getopt::Long;
#use Storable qw(store retrieve freeze thaw dclone); # uncomment when needed
use Carp;
use File::Basename;

use MTToolbox::Sample::SE;
use MTToolbox::Sample::PE::Overlap;
use MTToolbox::Sample::PE::NoOverlap;

use version; our $VERSION = qv('4.1.2');

sub _to_bool;

my $usage = "$0
--id --seq_type --fwd_file --rev_file --merged_file --barcode
--fwd_linker --rev_linker --fwd_primer --rev_primer --output_dir

[--fwd_mer_len] [--rev_mer_len] [--fwd_max_shifts] [--rev_max_shifts]
[--min_con_depth] [--diginorm_max] [--con_algo]

[--flash_m] [--flash_Max] [--flash_x] [--flash_p] [--flash_r] [--flash_f]
[--flash_s]

[--qc_con_trim_to_base] [--qc_con_min_len] [--qc_con_min_avg_qual]
[--qc_con_allow_ambig] [--qc_con_min_c_score] [--qc_SRC_trim_to_base]
[--qc_SRC_min_len] [--qc_SRC_min_avg_qual] [--qc_SRC_allow_gaps]
[--qc_SRC_allow_ambig]

[--merge_params_file] [--keep_tmp_files] [--debug] [--help]\n";

# NOTE: default variable are set up in the code

# Varaibles
my $id;
my $seq_type; # 1 == SE, 2 == PE with overlap, 3 == PE without overlap
my $barcode;
my $fwd_file;
my $rev_file;
my $merged_file;
my $fwd_linker;
my $rev_linker;
my $fwd_primer;
my $rev_primer;
my $fwd_mer_len;
my $rev_mer_len;
my $fwd_max_shifts;
my $rev_max_shifts;
my $min_con_depth;
my $diginorm_max;
my $con_algo;
my $output_dir;
my $flash_m;
my $flash_M;
my $flash_x;
my $flash_p;
my $flash_r;
my $flash_f;
my $flash_s;
my $qc_con_trim_to_base;
my $qc_con_min_len;
my $qc_con_min_avg_qual;
my $qc_con_allow_ambig;
my $qc_con_min_c_score;
my $qc_SRC_trim_to_base;
my $qc_SRC_min_len;
my $qc_SRC_min_avg_qual;
my $qc_SRC_allow_gaps;
my $qc_SRC_allow_ambig;
my $merge_params_file;
my $keep_tmp_files;
my $debug = 0;
my $help = 0;

GetOptions (
    "id:s" => \$id,
    "seq_type:i" => \$seq_type,
    "barcode:s" => \$barcode,
    "fwd_file:s" => \$fwd_file,
    "rev_file:s" => \$rev_file,
    "merged_file:s" => \$merged_file,
    "fwd_linker:s" => \$fwd_linker,
    "rev_linker:s" => \$rev_linker,
    "fwd_primer:s" => \$fwd_primer,
    "rev_primer:s" => \$rev_primer,
    "fwd_mer_len:i" => \$fwd_mer_len,
    "rev_mer_len:i" => \$rev_mer_len,
    "fwd_max_shifts:i" => \$fwd_max_shifts,
    "rev_max_shifts:i" => \$rev_max_shifts,
    "min_con_depth:i" => \$min_con_depth,
    "diginorm_max:s" => \$diginorm_max,
    "con_algo:s" => \$con_algo,
    "output_dir:s" => \$output_dir,
    "flash_m:i" => \$flash_m,
    "flash_Max:i" => \$flash_M,
    "flash_x:s" => \$flash_x,
    "flash_p:i" => \$flash_p,
    "flash_r:i" => \$flash_r,
    "flash_f:i" => \$flash_f,
    "flash_s:i" => \$flash_s,
    "qc_con_trim_to_base:s" => \$qc_con_trim_to_base,
    "qc_con_min_len:s" => \$qc_con_min_len,
    "qc_con_min_avg_qual:i" => \$qc_con_min_avg_qual,
    "qc_con_allow_ambig:s" => \$qc_con_allow_ambig,
    "qc_con_min_c_score:i" => \$qc_con_min_c_score,
    "qc_SRC_trim_to_base:s" => \$qc_SRC_trim_to_base,
    "qc_SRC_min_len:s" => \$qc_SRC_min_len,
    "qc_SRC_min_avg_qual:i" => \$qc_SRC_min_avg_qual,
    "qc_SRC_allow_gaps:s" => \$qc_SRC_allow_gaps,
    "qc_SRC_allow_ambig:s" => \$qc_SRC_allow_ambig,
    "merge_params_file:s" => \$merge_params_file,
    "keep_tmp_files:s" => \$keep_tmp_files,
    "debug" => \$debug,
    "help" => \$help,
);

if ( $help ) { die $usage; }


### Do some paramter checks
if ( ! defined $id              or
     ! defined $barcode         or
     ! defined $fwd_linker      or
     ! defined $rev_linker      or
     ! defined $fwd_primer      or
     ! defined $rev_primer      or
     ! defined $seq_type        or
     ! defined $output_dir
    ) {
    die $usage;
}

if (! -d $output_dir ) {
    mkdir $output_dir or die "Cannot create dir: $output_dir\nERROR: $!\n";
}

if ( $seq_type < 1 or $seq_type > 3 ) {
    die "seq_type must be either 1, 2, or 3.  1 == SE, 2 == PE with overlap, 3 == PE without overlap\n";
}

if ( $seq_type == 1 ) {
    if ( ! defined $fwd_file ) {
        craok("For seq_type == 1 a fwd_file must be provided\n");
    }
    if ( defined $merged_file ) {
        carp("For seq_type == 1 merged_file is ignored\n");
    }
    if ( defined $rev_file ) {
        carp("For seq_type == 1 rev_file is ignored\n");
    }
}

if ( $seq_type == 2 ) {
    if ( defined $merged_file and (defined $fwd_file or defined $rev_file) ) {
        carp("Using only merged_file for $id.  Other input files are ignored.\n");
    }
    elsif ( ! defined $merged_file and ! (defined $fwd_file and defined $rev_file) ){
        croak("For seq_type == 2 either a merged_file or (fwd_file, rev_file) pair must be provided\n");
    }
    elsif ( defined $fwd_file and ! defined $rev_file ) {
        croak("For seq_type == 2 a rev_file must be provided with the fwd_file");
    }
    elsif ( defined $rev_file and ! defined $fwd_file ) {
        croak("For seq_type == 2 a fwd_file must be provided with the rev_file");
    }
}

if ( $seq_type == 3 ) {
    if ( ! defined $fwd_file or ! defined $rev_file ) {
        croak("For seq_type == 3 a fwd_file and rev_file must be provided\n");
    }
    if ( defined $merged_file ) {
        carp("For seq_type == 3 merged_file is ignored\n");
    }
}

# Set default fwd_mer_len or rev_mer_len if not set by user
# You don't have to worry about this if ran from MTToolbox.pm
if ( ! defined $fwd_mer_len ) {
    $fwd_mer_len = 9;
}

if ( ! defined $rev_mer_len ) {
    $rev_mer_len = 4;
}

# Set the default fwd_max_shifts nad rev_max_shifts if not set by user
if ( ! defined $fwd_max_shifts ) {
    $fwd_max_shifts = 5;
}

if ( ! defined $rev_max_shifts ) {
    $rev_max_shifts = 5;
}

# Set defualt for min_con_depth
if ( ! defined $min_con_depth ) {
    $min_con_depth = 2;
}

# Set default for diginorm_max
if ( ! defined $diginorm_max ) {
    $diginorm_max = "NA";
}

if ( ! defined $con_algo ) {
    $con_algo = 'NoMSA';
}

# Set default for keep_tmp_files
if ( ! defined $keep_tmp_files ) {
    $keep_tmp_files = 0;
}
$keep_tmp_files = _to_bool($keep_tmp_files); # translate from string to bool
if ( ! ($keep_tmp_files == 0 or $keep_tmp_files == 1) ) {
    croak("--keep_tmp_files ($keep_tmp_files) must be 0 (do NOT keep) or 1 (keep)");
}



# Build the correct sample based on the input variable $seq_type
warn "### Building Sample Objects ###\n\n";
my $sample;
if ( $seq_type == 1 ) {
    $sample = MTToolbox::Sample::SE->new({
                                id => $id,
                                name => fileparse($fwd_file,  '\..*'),
                                barcode => $barcode,
                                fwd_linker => $fwd_linker,
                                fwd_primer => $fwd_primer,
                                rev_primer => $rev_primer,
                                fwd_mer_len => $fwd_mer_len,
                                fwd_max_shifts => $fwd_max_shifts,
                                min_con_depth => $min_con_depth,
                                diginorm_max => $diginorm_max,
                                output_dir => $output_dir,
                                reads_file => $fwd_file,
                            });
}
elsif ( $seq_type == 2 and defined $merged_file ) {
    $sample = MTToolbox::Sample::PE::Overlap->new({
                                        id => $id,
                                        name => fileparse($merged_file,  '\..*'),
                                        barcode => $barcode,
                                        fwd_linker => $fwd_linker,
                                        rev_linker => $rev_linker,
                                        fwd_primer => $fwd_primer,
                                        rev_primer => $rev_primer,
                                        fwd_mer_len => $fwd_mer_len,
                                        rev_mer_len => $rev_mer_len,
                                        fwd_max_shifts => $fwd_max_shifts,
                                        rev_max_shifts => $rev_max_shifts,
                                        min_con_depth => $min_con_depth,
                                        diginorm_max => $diginorm_max,
                                        output_dir => $output_dir,
                                        reads_file => $merged_file,
                                        });
}
elsif ( $seq_type == 2 and ! defined $merged_file ) {
    $sample = MTToolbox::Sample::PE::Overlap->new({
                                        id => $id,
                                        name => fileparse($fwd_file,  '\..*'),
                                        barcode => $barcode,
                                        fwd_linker => $fwd_linker,
                                        rev_linker => $rev_linker,
                                        fwd_primer => $fwd_primer,
                                        rev_primer => $rev_primer,
                                        fwd_mer_len => $fwd_mer_len,
                                        rev_mer_len => $rev_mer_len,
                                        fwd_max_shifts => $fwd_max_shifts,
                                        rev_max_shifts => $rev_max_shifts,
                                        min_con_depth => $min_con_depth,
                                        diginorm_max => $diginorm_max,
                                        output_dir => $output_dir,
                                        fwd_file => $fwd_file,
                                        rev_file => $rev_file,
                                        flash_m => $flash_m,
                                        flash_M => $flash_M,
                                        flash_x => $flash_x,
                                        flash_p => $flash_p,
                                        flash_r => $flash_r,
                                        flash_f => $flash_f,
                                        flash_s => $flash_s,
                                        });
    
    # FLASH params are generally included in the xml config
    # instead of a seperate merge_params_file
    $sample->merge_with_flash();
}
elsif ( $seq_type == 3 ) {
    $sample = MTToolbox::Sample::PE::NoOverlap->new({
                                            id => $id,
                                            name => fileparse($fwd_file,  '\..*'),
                                            barcode => $barcode,
                                            fwd_reads_file => $fwd_file,
                                            rev_reads_file => $rev_file,
                                            fwd_linker => $fwd_linker,
                                            rev_linker => $rev_linker,
                                            fwd_primer => $fwd_primer,
                                            rev_primer => $rev_primer,
                                            fwd_mer_len => $fwd_mer_len,
                                            rev_mer_len => $rev_mer_len,
                                            fwd_max_shifts => $fwd_max_shifts,
                                            rev_max_shifts => $rev_max_shifts,
                                            min_con_depth => $min_con_depth,
                                            diginorm_max => $diginorm_max,
                                            output_dir => $output_dir,
                                        });
}
else {
    die "seq_type must be either 1, 2, or 3.  1 == SE, 2 == PE with overlap, 3 == PE without overlap\n"
}

#print "outputdir: ", $sample->getOutputDir(), "\n";

### BINNING ###
warn "### BINNING ###\n\n";
# Bin the sample reads by Molecule Tag (MT)
$sample->categorize_by_MTs();

# Print categorizable reads
warn "### Print Categorizable reads ###\n\n";
$sample->print_categorizable_reads();

# Build MSAs
warn "### Build MSAs ###\n\n";
if ( $con_algo =~ m/NoMSA/i ) { warn "SKIP Build MSAs\n\n"; }
else { $sample->build_MSAs($con_algo); }

# Build Consensi
warn "### Build Consensi ###\n\n";
$sample->build_consensi($con_algo);

# Print Consensusi
warn "### Print Consensi ###\n\n";
$sample->print_consensi();

# Print single read categories
warn "### Print Single Read Categories ###\n\n";
$sample->print_single_read_categories();

# Print the MT distribution
warn "### Print MT Depth Output ###\n\n";
$sample->print_MT_depth_output();

# Print the unmatched seqeunces
warn "### Print Mismatch Seqs ###\n\n";
$sample->print_mismatch_seqs();

# Print the summary information
warn "### Print Sample Summaries ###\n\n";
$sample->print_sample_summary();

# Store this object for later use after the clustalw alignment step
#my $objectFile = "$output_dir/" . $sample->getBarcode() . ".obj";
#store($sample, $objectFile);


### Quality Control ###
warn "### Quality Control ###\n\n";
my $qc_params_href = {con_trim_to_base => $qc_con_trim_to_base,
                      con_min_len => $qc_con_min_len,
                      con_min_avg_qual => $qc_con_min_avg_qual,
                      con_allow_ambig => $qc_con_allow_ambig,
                      con_min_c_score => $qc_con_min_c_score,
                      SRC_trim_to_base => $qc_SRC_trim_to_base,
                      SRC_min_len => $qc_SRC_min_len,
                      SRC_min_avg_qual => $qc_SRC_min_avg_qual,
                      SRC_allow_gaps => $qc_SRC_allow_gaps,
                      SRC_allow_ambig => $qc_SRC_allow_ambig};
$sample->qc($qc_params_href);

### Delete Temp Files If Needed ###
if ( ! $keep_tmp_files ) {
    $sample->rm_tmp_files();
}

### FINISHED ###
warn "### FINISHED ###\n\n";


### Helper Subroutines ###
sub _to_bool {
    my ($parallelize) = @_;
    
    if ( $parallelize eq 1 or $parallelize eq 0 ) {
        return $parallelize;
    }
    
    my %good_yes_values = map { $_ => 1 } qw(Y YES Yes y yes);
    if ( defined $good_yes_values{$parallelize} ) {
        return 1;
    }
    
    # else -- meaning no parallel
    return 0;
}


__END__

# POD

=head1 NAME

bin_sample.pl - Bins sequences by molecule tag for one sample


=head1 VERSION

This documentation refers to bin_sample.pl version 4.1.2


=head1 USAGE

# Serial usage
perl bin_sample.pl

--id --seq_type --fwd_file --rev_file --merged_file --barcode
--fwd_linker --rev_linker --fwd_primer --rev_primer --output_dir

[--fwd_mer_len] [--rev_mer_len] [--fwd_max_shifts] [--rev_max_shifts]
[--min_con_depth] [--diginorm_max] [--con_algo]

[--flash_m] [--flash_Max] [--flash_x] [--flash_p] [--flash_r] [--flash_f]
[--flash_s]

[--qc_con_trim_to_base] [--qc_con_min_len] [--qc_con_min_avg_qual]
[--qc_con_allow_ambig] [--qc_con_min_c_score] [--qc_SRC_trim_to_base]
[--qc_SRC_min_len] [--qc_SRC_min_avg_qual] [--qc_SRC_allow_gaps]
[--qc_SRC_allow_ambig]

[--merge_params_file] [--keep_tmp_files] [--debug] [--help]

# Parallel usage
bsub -q week -n 1 -o bin_sample.log -e bin_sample.err
perl bin_sample.pl
--id --seq_type --fwd_file --rev_file --merged_file --barcode
--fwd_linker --rev_linker --fwd_primer --rev_primer --output_dir

[--fwd_mer_len] [--rev_mer_len] [--fwd_max_shifts] [--rev_max_shifts]
[--min_con_depth] [--diginorm_max] [--con_algo]

[--flash_m] [--flash_Max] [--flash_x] [--flash_p] [--flash_r] [--flash_f]
[--flash_s]

[--qc_con_trim_to_base] [--qc_con_min_len] [--qc_con_min_avg_qual]
[--qc_con_allow_ambig] [--qc_con_min_c_score] [--qc_SRC_trim_to_base]
[--qc_SRC_min_len] [--qc_SRC_min_avg_qual] [--qc_SRC_allow_gaps]
[--qc_SRC_allow_ambig]

[--merge_params_file] [--keep_tmp_files] [--debug] [--help]\n";


=head1 REQUIRED ARGUMENTS

=head2 --id

    A unique ID for this sample.  This ID will be incorperated into the headers
    of the output sequences.

=head2 --seq_type

    This parameter is a coded integer representing the type of sequences in the
    sequence files.
    
    1 => Single-end, non-overlapping reads
    2 => Paired-end, overlapping reads
    3 => Paired-end, non-overlapping reads

=head2 --fwd_file

    A path to the file contianing the forward seqeucnes.  If the reads are single
    end the --fwd_file parameter is the path to those reads.
    
=head2 --rev_file

    A path to the file containing the reverse sequences.  This parameter should
    only be provided if --seq_type is either 2 or 3.
    
=head2 --barcode

    The barcode sequence of this sample.  This sequence is used to name output
    directories.
    
=head2 --fwd_linker

    The forward linker sequence.
    
=head2 --rev_linker

    The reverse linker sequence.
    
=head2 --fwd_primer

    A regular expression to match the forward primer.  A valid input might look
    something like this (include quotes in command):
    
    "GTGCCAGC[AC]GCCGCGGTAA"
    
    A useful web tool for testing regular expression is: www.regextester.com/
    
=head2 --rev_primer

    A regular expression to match the reverse primer.  A valid input might look
    something like this (include quotes in command):
    
    "GGACTAC[ACT][ACG]GGGT[AT]TCTAAT"
    
    A useful web tool for testing regular expression is: www.regextester.com/
    

=head1 OPTIONAL ARGUMENTS

=head2 --fwd_mer_len

    The length of the random mer on the forward read.  This must be an integer
    > 0.  DEFAULT = 9.
    
=head2 --rev_mer_len

    The length of the random mer on the forward read.  This must be an integer
    > 0.  DEFAULT = 4.
    
=head2 --fwd_max_shifts

    The maximum number of frameshifts included in the fwd tag.  This must be an
    integer >= 0.  DEFAULT = 5.
    
=head2 --rev_max_shifts

    The maximum number of frameshifts included in the rev tag.  This must be an
    integer >= 0.  DEFAULT = 5.
    
=head2 --min_con_depth

    When building the consensus sequences min_con_depth is the minimum number
    of reads in a molecule tag category in order to build a consensus sequence.
    For example, if min_con_depth is set to 2 then a molecule tag category needs
    at least 2 sequences in order to build a consensus.  To have more correct
    sequences increase min_con_depth.  However, when min_con_depth is increased
    the total number of resulting consensus sequences is reduced.
    integer >= 2.  DEFAULT = 2.

=head2 --diginorm_max

    Some molecule tag categories can have hundreds of sequences.  With a large
    number of sequences the MSA takes a long time.  In practice such a larege
    number of sequences is not necessary to make a very accurate consensus.
    The diginorm_max parameter is the maximum number of reads to use for
    building a consensus sequence.
    
=head2 --con_algo

    The algorithm to use when building a consensus seuqnece.  It is generally
    recommended that this parameter be set to MSA meaning use the multiple
    sequence alignment generated by clustalw.

=head2 --output_dir

    A path to the directory where output files and directories can be printed.
    If the given output_dir does not exists it is automatically made.

=head2 --merge_params_file

    A file with the flash parameters.  DEPRECIATED.  It is now recommended to
    put the flash parameters in the config.xml file.

=head2 --flash_m

    The min overlap flash parameter.

=head2 --flash_Max

    The max overlap flash parameter.  I had to change this from flash_M to flash_Max
    because Getopt::Long is case insensitive.

=head2 --flash_x

    The ratio flash parameter.

=head2 --flash_p

    The phred flash parameter.

=head2 --flash_r

    The read length flash parameter.

=head2 --flash_f

    The fragment length flash parameter.

=head2 --flash_s

    The standard deviation flash parameter.
    
=head2 --qc_con_min_len
    
    The minimum length for high quality consensus sequences
    
=head2 --qc_con_min_avg_qual

    The minimum average quality value for high quality consensus sequences
    
=head2 --qc_con_allow_ambig

    Y or N indicating whether to allow ambiguous bases in high quality consensus
    sequences

=head2 --qc_min_c_score
    
    The minimum c-score for high quality consensus sequences.  C-score is
    calculated as follows:  sum the quotient of each of the consensus sequence
    bases and total depth for that consensus sequence for each base.
    i.e sum(bases[1..end] / depth)
    
=head2 --qc_SRC_min_len
    
    The minimum length for high quality single read categories
    
=head2 --qc_SRC_min_avg_qual

    The minimum average quality value for high quality single read categories
    
=head2 --qc_SRC_allow_gaps

    Y or N indicating whether to allow gaps in high quality single read
    categories
    
=head2 --qc_SRC_allow_ambig

    Y or N indicating whether to allow ambiguous bases in high quality single
    read categories
    
=head2 --keep_tmp_files

    Set to N to remove tmp files.  Set to Y to keep tmp files.  DEFAULT: N
    This variable is not stored in the sample object.  It is simply used in
    bin_sample.pl to optionally remove temporary files.

=head2 --debug

    An optional parameter to print debugging statments as the script runs.
    DEPRECIATED.
    
=head2 --help

    An optional parameter to print a usage statement.

=head1 DESCRIPTION

Limitting the number of sequencing errors is crutial for precise sequence
comparisons.  Sequencing errors are often introduced during the PCR steps
and/or the sequencing step.  In an attempt to limit the number of errors
we have developed a molecular biology technique in which individual DNA
molecules are tagged with a unique barcode-like sequences.  We call this
sequence a molecule tag (MT).  After each DNA molecule has been tagged it
can be PCRed and sequenced.  After sequencing we recieve raw reads that
contain a molecule tag at the beginning.  MTToolbox.pl sorts those sequences
based on their molecule tags and builds a consensus sequence to represent
the original DNA molecule that was tagged.

For large sequence sets across several samples, these operations can be quite
time intensive so MTToolbox.pl includes an option to run in parallel using LSF.  

=head2 INPUTS

See REQUIRED ARGUMENTS


=head2 OUTPUTS
    
Inside the specified output_dir a directory named using the barcode sequence
is created.  This directory contains all output files for this sample
generated by bin_sample.pl.

The main output file is called consensusSeqs.fastq which contains all the
consensus sequences for all the molecule tag bins.  Another fastq file
called mismatched.fastq contains all the fastq reads that did not match
the specified pattern.  A summary file (summary.txt) contains summary
information about the molecule tags and other summary numbers for the entire
sample.  A file called seqs_per_mt_dist.txt contains the distribution of the
number of reads found in each molecule tag bin.  A file called
categorizable_reads.fastq is a fastq file of all the reads that match the
expected pattern and can therefore be classified in catigories by their
molecule tag.  A file called single_read_categories.fastq is a fastq formated
file that contains all the reads that are categorizable but only have one
read that fits in that category.

There are two intermediate directories that are created:  fasta_files and
aln_files.  The fasta_files directory contains the fasta sequences split by
their molecule tag sequences.  These file are the inputs to clustalw to
build a multiple sequence alignment.  The aln_files directory contains the
outputs of running clustalw.  These are the multiple sequence alignment files
which are parsed to build the consensus sequence for each molecule tag bin.
    

=head1 CONFIGURATION AND ENVIRONMENT
    
Tested using perl5.8.9 and perl5.12.3.
    
    
=head1 DEPENDANCIES

    Getopt::Long
    Carp
    File::Basename
    MTToolbox::Sample::SE
    MTToolbox::Sample::PE::Overlap
    MTToolbox::Sample::PE::NoOverlap
    version
    

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to author


=head1 AUTHOR

Scott Yourstone  C<< <scott.yourstone81@gmail.com> >>


=head1 LICENCE AND COPYRIGHT

Copyright (c) 2013, Scott Yourstone
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.


=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.


=cut
