#! /usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Carp;
use version; our $VERSION = qv('4.1.2');

use MTToolbox 4.1.2;

# USAGE
my $usage = "$0 --params_file file [--help]\n";

# Variables to be set by GetOpts and their defaults
my $params_file = undef;
my $debug = 0;
my $help = 0;


my $options_okay = GetOptions (
    "params_file:s" => \$params_file,  # string
    "help" => \$help,                  # flag
);

# Fail if unknown options are encountered
if ( ! $options_okay ) {
    die $usage;
}

# Print a help message
if ( $help ) {
    die $usage;
}

# Fail if no params file is given
if ( ! defined $params_file ) {
    die $usage;
}

my $toolbox = MTToolbox->new({params_file => $params_file});
$toolbox->run();


__END__

# POD

=head1 NAME

MTDriver.pl - Runs the guts of a MTToolbox analysis


=head1 VERSION

This documentation refers to MTDriver.pl version 4.1.2


=head1 USAGE

MTDriver.pl --params_file my_params.xml


=head1 REQUIRED ARGUMENTS

=head2 --params_file
    
    A XML formated paramters file.  Descriptions and examples of such a file can
    be found on the MTToolbox webpage
    (https://sites.google.com/site/moleculetagtoolbox/).
    
    Also, the params file can be generated by using the MTToolbox GUI which is
    run by typing the command "MTToolbox.pl".
    
    
=head1 OPTIONAL ARGUMENTS

NA


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


=head1 OUTPUTS

The following is a list and summary of some of the major output files:

all_summary.txt - Summary numbers for the entire run and each sample
individually.

all_summary.ps - PostScript version of the above text file.

config.xml - The config file with all the run variables.

all_consensusSeqs_LQ.fastq - Consensus sequences that failed the QC test.

all_consensusSeqs_diagnostics.txt - Reasons why each of the low quality reads
failed the QC tests.  If large numbers of sequences failed the QC tests this
file can be used to diagnose what went wrong.

all_single_read_categories_HQ.fastq - MT categories with only one raw read.  It
is imposible to build a consensus sequences from just one read, however it these
may still be useful in downstream analysis.

all_single_read_categories_LQ.fastq - MT categories with only one raw read that
were classified as low quality

all_single_read_categories_diagnostics.txt - Reasons why each of the low quality
reads failed the QC tests.

all_categorizable_reads.fastq - All raw reads which match the expected pattern.

all_seqs_per_mt_dist.txt - Distribution of raw reads per MT category text file

all_seqs_per_mt_hist.txt - The above text file in histogram format

all_seqs_per_mt_hist.ps - The PostScript version of the above file

con_ambig_base_comp_dist.png - Graph of the kinds of ambiguous bases in the
consensus seqs

con_ambig_base_pos_dist.png - Graph of where the ambiguous base position in the
consensus seqs.

con_ambig_base_per_seq.png - Graph showing the distribution of ambiguous bases
in each read.
    

=head1 CONFIGURATION AND ENVIRONMENT
    
Tested using perl5.8.9 and perl5.12.3.
    
    
=head1 DEPENDANCIES

    Getopt::Long
    Carp
    version
    MTToolbox 4.1.2
    

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
