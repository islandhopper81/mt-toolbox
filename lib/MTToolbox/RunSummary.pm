package MTToolbox::RunSummary;

# This object is a data structure that holds summary information for
# an entire run

use strict;
use warnings;

use MTToolbox::SampleSummary;
use MTToolbox::SampleSummary::SE;
use MTToolbox::SampleSummary::PE::NoOverlap;
use MTToolbox::SampleSummary::PE::Overlap;

use Chart::Gnuplot;
use Class::Std::Utils;
use Carp qw(carp croak);
use version; our $VERSION = qv('4.1.2');

{
    
    # Attributes #
    my %sample_count_of;
    my %total_seqs_count_of;
    my %merged_seqs_count_of;
    my %merged_percent_of;
    my %match_count_of;
    my %mismatch_count_of;
    my %mt_count_of;
    my %SRC_count_of;
    my %samples_href_of;
    
    # Setters #
    sub set_sample_count;
    sub set_total_seqs_count;
    sub set_merged_seqs_count;
    sub set_match_count;
    sub set_mismatch_count;
    sub set_mt_count;
    sub set_SRC_count;
    sub set_samples_href;
    
    # Getters #
    sub get_sample_count;
    sub get_total_seqs_count;
    sub get_merged_seqs_count;
    sub get_not_merged_seqs_count;
    sub get_merged_percent;
    sub get_match_count;
    sub get_mismatch_count;
    sub get_mt_count;
    sub get_SRC_count;
    sub get_sample_href;
    
    # Others #
    sub to_string;
    sub load_sample_summaries;
    sub print_summary_files;
    sub _print_summary_graph;
    sub _get_sample_name;
    
    
    ###############
    # Constructor #
    ###############
    sub new {
        my ($class) = @_;

        # Croak if calling 'new' on an already blessed reference
        croak 'Constructor called on existing object instead of class'
            if ref $class;
        
        # Bless a scalar to instantiate an object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        # Initialize the objects attributes        
        
        return $new_obj;
    }
    
    
    ###########
    # Setters #
    ###########
    sub set_sample_count {
        my ($self, $c) = @_;
        $sample_count_of{ident $self} = $c;
        return 1;
    }
    
    sub set_total_seqs_count {
        my ($self, $c) = @_;
        $total_seqs_count_of{ident $self} = $c;
        return 1;
    }
    
    sub set_merged_seqs_count {
        my ($self, $c) = @_;
        $merged_seqs_count_of{ident $self} = $c;
        return 1;
    }
    
    sub set_match_count {
        my ($self, $c) = @_;
        $match_count_of{ident $self} = $c;
        return 1;
    }
    
    sub set_mismatch_count {
        my ($self, $c) = @_;
        $mismatch_count_of{ident $self} = $c;
        return 1;
    }
    
    sub set_mt_count {
        my ($self, $c) = @_;
        $mt_count_of{ident $self} = $c;
        return 1;
    }
    
    sub set_SRC_count {
        my ($self, $c) = @_;
        # SRC - single read categories
        $SRC_count_of{ident $self} = $c;
        return 1;
    }
    
    sub set_samples_href {
        my ($self, $href) = @_;
        $samples_href_of{ident $self} = $href;
        return 1;
    }
    
    ###########
    # Getters #
    ###########
    sub get_sample_count {
        my ($self) = @_;
        
        if ( my $c = $sample_count_of{ident $self} ) {
            return $c;
        }
        
        return "--";
    }
    
    sub get_total_seqs_count {
        my ($self) = @_;
        
        my $count = $total_seqs_count_of{ident $self};
        if ( defined $count ) {
            return $count;
        }
        
        return "--";
    }
    
    sub get_merged_seqs_count {
        my ($self) = @_;
        
        my $count = $merged_seqs_count_of{ident $self};
        if ( defined $count ) {
            return $count;
        }
        
        return "--";
    }
    
    sub get_not_merged_seqs_count {
        my ($self) = @_;
        
        my $total_seqs = $self->get_total_seqs_count();
        my $merged_seqs = $self->get_merged_seqs_count();
        if ( $total_seqs eq "--" or
             $merged_seqs eq "--" ) {
            return "--";
        }
        
        return $total_seqs - $merged_seqs;
    }
    
    sub get_merged_percent {
        my ($self) = @_;
        
        my $merged_seqs = $self->get_merged_seqs_count();
        if ( $merged_seqs =~ m/^(\d+)$/ ) {
            return int($1 / $self->get_total_seqs_count() * 100) . "%";
        }
        
        return "--";
    }
    
    sub get_match_count {
        my ($self) = @_;
        
        my $count = $match_count_of{ident $self};
        if ( defined $count ) {
            return $count;
        }
        
        return "--";
    }
    
    sub get_mismatch_count {
        my ($self) = @_;
        
        my $count = $mismatch_count_of{ident $self};
        if ( defined $count ) {
            return $count;
        }
        
        return "--";
    }
    
    sub get_mt_count {
        my ($self) = @_;
        
        my $count = $mt_count_of{ident $self};
        if ( defined $count ) {
            return $count;
        }
        
        return "--";
    }
    
    sub get_SRC_count {
        my ($self) = @_;
        
        my $count = $SRC_count_of{ident $self};
        if ( defined $count ) {
            return $count;
        }
        
        return "--";
    }
    
    sub get_samples_href {
        my ($self) = @_;
        
        return $samples_href_of{ident $self};
    }
    
    
    
    ##########
    # Others #
    ##########
    sub to_string {
        my ($self) = @_;
        
        my $total_seqs_count = $self->get_total_seqs_count();
        my $merged_seqs_count = $self->get_merged_seqs_count();
        my $not_merged_seqs_count = $total_seqs_count - $merged_seqs_count;
        
        my $str = "Sample count: " . $self->get_sample_count() . "\n";
        $str .= "Seq count: " . $total_seqs_count . "\n";
        $str .= "Merged count: " . $merged_seqs_count . "\n";
        $str .= "Not merged count: " . $not_merged_seqs_count . "\n";
        $str .= "Percent merged: " . $self->get_merged_percent() . "\n";
        $str .= "Match count: " . $self->get_match_count() . "\n";
        $str .= "Mismatch count: " . $self->get_mismatch_count() . "\n";
        $str .= "Molecule Tag (MT) count: " . $self->get_mt_count() . "\n";
        $str .= "Single read category count: " . $self->get_SRC_count() . "\n";
        
        $str .= "\n\n";  # this is important to seperate sections for gnuplot
        
        # add the indvidual sample summary information
        $str .= "# Name\tID\tTotal_Seq_Count\tMerged_Count\tNot_Merged_Count\t";
        $str .= "Percent_Merged\tMatch_Count\tMismatch_Count\tMT_Count\tSRC_Count\n";
        
        my $href = $self->get_samples_href();
        
        foreach my $id ( sort keys %{$href} ) {
            
            my $summ_obj = $href->{$id};
            
            $str .= $summ_obj->get_sample_name() . "\t";
            $str .= $summ_obj->get_sample_id() . "\t";
            
            # get the merged counts if necessary.  Only the
            # MTToolbox::SampleSummary::PE::Overlap objects have these methods so
            # I have to check the blessing of the object
            my $merged_count = "NA";
            my $not_merged_count = "NA";
            my $merged_perc = "NA";
            if ( ref $summ_obj eq "MTToolbox::SampleSummary::PE::Overlap" ) {
                $merged_count = $summ_obj->get_merged_count();
                $not_merged_count = $summ_obj->get_not_merged_count();
                $merged_perc = $summ_obj->get_perc_merged();
            }
            
            $str .= $summ_obj->get_seq_count() . "\t";
            $str .= $merged_count . "\t";
            $str .= $not_merged_count . "\t";
            $str .= $merged_perc . "\t";
            $str .= $summ_obj->get_match_count() . "\t";
            $str .= $summ_obj->get_mismatch_count() . "\t";
            $str .= $summ_obj->get_MT_count() . "\t";
            $str .= $summ_obj->get_SRC_count() . "\n";
        }

        return $str;                                                        
    }
    
    sub load_sample_summaries {
        my ($self, $output_root) = @_;
        
        # stores sample summary objects.
        # KEY: sample name  VALUE: MTToolbox::SampleSummary
        my %samples_hash = ();
        
        # global variables to calculate and store
        my $sample_count = 0;
        my $total_seqs_count = 0;
        my $merged_seqs_count = 0;
        my $percent_merged = 0;
        my $match_count = 0;
        my $mismatch_count = 0;
        my $mt_count = 0;
        my $SRC_count = 0;
        
        # load the samples summaries by their file
        my @files = `find $output_root -maxdepth 2 -name \"summary.txt\"`;
        FILE: foreach my $file_path ( @files ) {
            my $summ_obj;
            #my $sample_name;
            my $sample_id;
            
            # build the right type of MTToolbox::SampleSummary object
            if ( `grep "PE w/ Overlap" $file_path` ) {
                $summ_obj = MTToolbox::SampleSummary::PE::Overlap->new();
            }
            elsif ( `grep "SE" $file_path` ) {
                $summ_obj = MTToolbox::SampleSummary::SE->new();
            }
            elsif ( `grep "PE w/o Overlap" $file_path` ) {
                $summ_obj = MTToolbox::SampleSummary::PE::NoOverlap->new();
            }
            else {
                my $message = "Missing sample type in summary file: $file_path\n";
                $message .= "Skipping summary file: $file_path";
                carp($message);
                next FILE;
            }
            
            # load the data from the file
            $summ_obj->load_summary_file($file_path);
            
            # get the sample id
            #$sample_name = _get_sample_name($file_path);
            $sample_id = $summ_obj->get_sample_id();
            
            # store the sample
            #$samples_hash{$sample_name} = $summ_obj;
            $samples_hash{$sample_id} = $summ_obj;
            
            # calculate global variables
            $sample_count++;
            $total_seqs_count += $summ_obj->get_seq_count();
            $match_count += $summ_obj->get_match_count();
            $mismatch_count += $summ_obj->get_mismatch_count();
            $mt_count += $summ_obj->get_MT_count();
            $SRC_count += $summ_obj->get_SRC_count();
            
            if ( `grep "PE w/ Overlap" $file_path` ) {
                $merged_seqs_count += $summ_obj->get_merged_count();
            }
        }
        
        # save the global variables.
        $self->set_sample_count($sample_count);
        $self->set_total_seqs_count($total_seqs_count);
        $self->set_merged_seqs_count($merged_seqs_count);
        $self->set_match_count($match_count);
        $self->set_mismatch_count($mismatch_count);
        $self->set_mt_count($mt_count);
        $self->set_SRC_count($SRC_count);
        $self->set_samples_href(\%samples_hash);
        
        return 1;
    }
    
    sub print_summary_files {
        my ($self, $output_root) = @_;
        
        # print the summary text file
        my $out_file = "$output_root/all_summary.txt";
        open (my $OUT, ">", $out_file) or
            croak("Cannot open file: $out_file");
        
        print $OUT $self->to_string();
        
        close($OUT);
        
        # print the summary graph
        $self->_print_summary_graph($output_root);        
        
        return 1;
    }
    
    sub _print_summary_graph {
        my ($self, $output_root) = @_;
        
        print "output_root: $output_root\n";
        my $in_file = "$output_root/all_summary.txt";
        my $out_file = "$output_root/all_summary.ps";
        
        my $com = <<END;
        gnuplot << EOF

        reset
        
        set term postscript
        set output "$out_file"
        set title "Sample Summary Counts"
        set xlabel "Samples"
        set ylabel "Count"
        set xtics rotate by 270
        set xtics norangelimit
        set offsets .05,.05
        
        
        # plot total seqs
        plot "$in_file" i 1 u 3:xtic(1) title "total seqs", \\
        "" i 1 u 4:xtic(1) title "merged seqs", \\
        "" i 1 u 7:xtic(1) title "match seqs", \\
        "" i 1 u 9:xtic(1) title "MTs", \\
        "" i 1 u 10:xtic(1) title "SRCs"
END
        # the above line need to be at the very beginning because it signifies
        # the end of the heredoc
        
        # run the gnuplot command
        system($com);
        
        return 1;
    }
    
    sub _get_sample_name {
        my ($name) = @_;
        
        # the name is of the form: output_20130207/P1_CGTCGGTA/summary.txt
        my @values = split /\//, $name;
        
        return $values[-2];
    }
}


1;
__END__

#######
# POD #
#######
=head1 MTToolbox::RunSummary

MTToolbox::RunSummary - A class to store and print summary inforamtion
for an entire run

=head1 VERSION

This documentation refers to MTToolbox::RunSummary version 4.1.2.

=head1 Included Modules

MTToolbox::SampleSummary
MTToolbox::SampleSummary::SE
MTToolbox::SampleSummary::PE::NoOverlap
MTToolbox::SampleSummary::PE::Overlap
Chart::Gnuplot
Class::Std::Utils
Carp qw(carp croak)
version

=head1 Inherit

NA

=head1 SYNOPSIS
    
    use MTToolbox::RunSummary
    my $my_sum_info = MTToolbox::RunSummary->new();
    
    $my_sum_info->load_sample_summaries($output_root);
    
    # Printing Methods -- See also MTToolbox::SampleSummary
    $my_sum_info->print_summary_files($output_root);
    

=head1 DESCRIPTION

This class stores and outputs summary information about a Sample object
with PE reads.

=head1 METHODS

=over

    new
    set_sample_count
    set_total_seqs_count
    set_merged_seqs_count
    set_match_count
    set_mismatch_count
    set_mt_count
    set_SRC_count
    set_samples_href
    get_sample_count
    get_total_seqs_count
    get_merged_seqs_count
    get_not_merged_seqs_count
    get_merged_percent
    get_match_count
    get_mismatch_count
    get_mt_count
    get_SRC_count
    get_samples_href
    to_string
    load_sample_summaries
    print_summary_files
    _print_summary_graph
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: MTToolbox::RunSummary->new();
    Function: Creates a new MTToolbox::RunSummary object
    Returns: MTToolbox::RunSummary
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA

=head2 set_sample_count

    Title: set_sample_count
    Usage: $my_sum_info->set_sample_count($sample_count);
    Function: Sets the number of samples in the run
    Returns: 1 on successful completion
    Args: -sample_count => the number of samples
    Throws: NA
    Comments: NA
    See Also: NA

=head2 set_total_seqs_count

    Title: set_total_seqs_count
    Usage: $my_sum_info->set_total_seqs_count($total_seqs);
    Function: Sets the number of raw reads from the sequencing file
    Returns: 1 on successful completion
    Args: -total_seqs => the number of total raw reads
    Throws: NA
    Comments: NA
    See Also: NA

=head2 set_merged_seqs_count

    Title: set_merged_seqs_count
    Usage: $my_sum_info->set_merged_seqs_count($merged_seqs_count);
    Function: Sets the number of reads merged by FLASH (or other assembler)
    Returns: 1 on successful completion
    Args: -merged_seqs_count => the number of reads that successfully merge
    Throws: NA
    Comments: NA
    See Also: NA

=head2 set_match_count

    Title: set_match_count
    Usage: $my_sum_info->set_match_count($match_count);
    Function: Sets the number of reads that match the regex pattern we expect
    Returns: 1 on successful completion
    Args: -match_count => the number of reads that match the regex
    Throws: NA
    Comments: NA
    See Also: NA

=head2 set_mismatch_count

    Title: set_mismatch_count
    Usage: $my_sum_info->set_mismatch_count($mismatch_count);
    Function: Sets the number of reads that DO NOT match the regex pattern we
              expect
    Returns: 1 on successful completion
    Args: -mismatch_count => the number of reads that DO NOT match the regex
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 set_mt_count

    Title: set_mt_count
    Usage: $my_sum_info->set_mt_count($mt_count);
    Function: Sets the number of molecule tags that are found
    Returns: 1 on successful completion
    Args: -mt_count => the number of molecule tags
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 set_SRC_count

    Title: set_SRC_count
    Usage: $my_sum_info->set_SRC_count($SRC_count);
    Function: Sets the number of single read categories (SRC)
    Returns: 1 on successful completion
    Args: -SRC_count => the number of single read categories
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 set_samples_href

    Title: set_samples_href
    Usage: $my_sum_info->set_samples_href($href);
    Function: Sets the hash ref that the MTToolbox::SampleSummary objects are
              stored in
    Returns: 1 on successful completion
    Args: -href => a hash ref with MTToolbox::SampleSummary objects
    Throws: NA
    Comments: KEY => sample_id, VALUE => MTToolbox::SampleSummary object
    See Also: NA

=head2 get_sample_count

    Title: get_sample_count
    Usage: $my_sum_info->get_sample_count();
    Function: Gets the number of samples in the run
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_total_seqs_count

    Title: get_total_seqs_count
    Usage: $my_sum_info->get_total_seqs_count();
    Function: Gets the number of raw reads
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_merged_seqs_count

    Title: get_merged_seqs_count
    Usage: $my_sum_info->get_merged_seqs_count();
    Function: Gets the number of reads that were merged by FLASH
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_not_merged_seqs_count

    Title: get_not_merged_seqs_count
    Usage: $my_sum_info->get_not_merged_seqs_count();
    Function: Gets the number of reads that were NOT merged by FLASH
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_merged_percent

    Title: get_merged_percent
    Usage: $my_sum_info->get_merged_percent();
    Function: Gets the percentage of reads merged by FLASH
    Returns: String
    Args: None
    Throws: NA
    Comments: The return value is an int with the '%'
    See Also: NA
    
=head2 get_match_count

    Title: get_match_count
    Usage: $my_sum_info->get_match_count();
    Function: Gets the number of reads that match the expected pattern
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_mismatch_count

    Title: get_mismatch_count
    Usage: $my_sum_info->get_mismatch_count();
    Function: Gets the number of reads that do NOT match the expected pattern
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_mt_count

    Title: get_mt_count
    Usage: $my_sum_info->get_mt_count();
    Function: Gets the number of molecule tags
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_SRC_count

    Title: get_SRC_count
    Usage: $my_sum_info->get_SRC_count();
    Function: Gets the number of single read categories
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_samples_href

    Title: get_samples_href
    Usage: $my_sum_info->get_samples_href();
    Function: Gets the hash ref with the MTToolbox::SampleSummary objects
    Returns: href
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 load_sample_summaries

    Title: load_sample_summaries
    Usage: $my_sum_info->load_sample_summaries($output_root);
    Function: Saves all summary info for each sample based
    Returns: 1 on successful completion
    Args: -output_root => root of output file of MTToolbox analysis
    Throws: NA
    Comments: Sample summary files should be 2 directories underneath the given
              output_root directory and they should be named summary.txt.  
    See Also: NA
    
=head2 print_summary_files

    Title: print_summary_files
    Usage: $my_sum_info->print_summary_files($output_root);
    Function: Prints all run summary files inluding graphs
    Returns: 1 on successful completion
    Args: -output_root => root of output file of MTToolbox analysis
    Throws: NA
    Comments: NA
    See Also: NA

=head2 to_string

    Title: to_string
    Usage: $my_sum_info->to_string();
    Function: Gets the summary information stored in this object as a string
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 _print_summary_graph

    Title: _print_summary_graph
    Usage: _print_summary_graph(-output_root);
    Function: Prints the sample summary graph
    Returns: 1 on successful completion
    Args: -output_root => output root of MTToolbox output
    Throws: NA
    Comments: This method sends a system command to gnuplot through a command
              line call.  It outputs a postscript file in the output_root dir.
    See Also: NA

=head1 AUTHOR

Scott Yourstone     scott.yourstone81@gmail.com


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