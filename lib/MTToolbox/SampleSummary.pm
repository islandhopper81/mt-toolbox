package MTToolbox::SampleSummary;

# This object is a data structure that holds summary information for each sample.
# This is the parent object for MTToolbox::SampleSummary::SE.pm,
# MTToolbox::SampleSummary::PE::Overlap, and
# MTToolbox::SampleSummary::PE::NoOverlap.


use strict;
use warnings;

use Class::Std::Utils;
use Carp qw(carp croak);
use version; our $VERSION = qv('4.1.2');

{
    
    # Attributes #
    my %sample_name_of;
    my %sample_id_of;
    my %file_format_of;
    my %match_count_of;
    my %mismatch_count_of;
    my %seq_count_of;
    my %mt_count_of;
    my %SRC_count_of;

    # Setters #
    sub set_sample_name;
    sub set_sample_id;
    sub set_file_format;
    sub set_match_count;
    sub set_mismatch_count;
    sub set_seq_count;
    sub set_MT_count;
    sub set_SRC_count;
    
    # Getters #
    sub get_sample_name;
    sub get_sample_id;
    sub get_file_format;
    sub get_match_count;
    sub get_mismatch_count;
    sub get_seq_count;
    sub get_MT_count;
    sub get_SRC_count;

    # Others #
    sub print_summary;
    sub to_string;
    sub load_summary_file;
    
    
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
    sub set_sample_name {
        my ($self, $n) = @_;
        $sample_name_of{ident $self} = $n;
        return 1;
    }
    
    sub set_sample_id {
        my ($self, $id) = @_;
        $sample_id_of{ident $self} = $id;
        return 1;
    }
    
    sub set_file_format {
        my ($self, $f) = @_;
        $file_format_of{ident $self} = $f;
        return 1;
    }
    
    sub set_match_count {
        my ($self, $m) = @_;
        $match_count_of{ident $self} = $m;
        return 1;
    }
    
    sub set_mismatch_count {
        my ($self, $u) = @_;
        $mismatch_count_of{ident $self} = $u;
        return 1;
    }
    
    sub set_seq_count {
        my ($self, $t) = @_;
        $seq_count_of{ident $self} = $t;
        return 1;
    }
    
    sub set_MT_count {
        my ($self, $m) = @_;
        $mt_count_of{ident $self} = $m;
        return 1;
    }
    
    sub set_SRC_count {
        my ($self, $s) = @_;
        # SRC stands for single_read_categories
        $SRC_count_of{ident $self} = $s;
        return 1;
    }
    
    ###########
    # Getters #
    ###########
    sub get_sample_name {
        my ($self) = @_;
        
        my $sample_name = $sample_name_of{ident $self};
        if ( defined $sample_name ) {
            return $sample_name;
        }
        
        return "--";
    }
    
    sub get_sample_id {
        my ($self) = @_;
        
        my $sample_id = $sample_id_of{ident $self};
        if ( defined $sample_id ) {
            return $sample_id;
        }
        
        return "--";
    }
    
    sub get_file_format {
        my ($self) = @_;
        
        my $file_format = $file_format_of{ident $self};
        if ( defined $file_format ) {
            return $file_format;
        }
        
        return "--";  # In the case that file_format_of is not defined
    }
    
    sub get_match_count {
        my ($self) = @_;
        
        my $match_count = $match_count_of{ident $self};
        if ( defined $match_count ) {
            return $match_count;
        }
        
        return "--";
    }
    
    sub get_mismatch_count {
        my ($self) = @_;
        
        my $mismatch_count = $mismatch_count_of{ident $self};
        if ( defined $mismatch_count ) {
            return $mismatch_count;
        }
        
        return "--";
    }
    
    sub get_seq_count {
        my ($self) = @_;
        
        my $seq_count = $seq_count_of{ident $self};
        if ( defined $seq_count ) {
            return $seq_count;
        }
        
        return "--";
    }
    
    sub get_MT_count {
        my ($self) = @_;
        
        my $mt_count = $mt_count_of{ident $self};
        if ( defined $mt_count ) {
            return $mt_count;
        }
        
        return "--";
    }
    
    sub get_SRC_count {
        my ($self) = @_;
        
        #SRC stands for single_read_categories
        my $SRC_count = $SRC_count_of{ident $self};
        if ( defined $SRC_count ) {
            return $SRC_count;
        }
        
        return "--";
    }
    
    
    ##########
    # Others #
    ##########
    sub print_summary {
        my ($self, $outputFile) = @_;
        
        open (SUM, ">$outputFile") or croak("Cannot open $outputFile\nERROR: $!\n");
        
        print SUM $self->to_string();
        
        close(SUM);
        
        return 1;
    }
    
    sub to_string {
        my ($self) = @_;
        
        my $str = "Sample name: " . $self->get_sample_name() . "\n";
        $str .= "ID: " . $self->get_sample_id() . "\n";
        $str .= "File format: " . $self->get_file_format() . "\n";
        $str .= "Seq count: " . $self->get_seq_count() . "\n";
        $str .= "Match count: " . $self->get_match_count() . "\n";
        $str .= "Mismatch count: " . $self->get_mismatch_count() . "\n";
        $str .= "Molecule Tag (MT) count: " . $self->get_MT_count() . "\n";
        $str .= "Single read category count: " . $self->get_SRC_count() . "\n";
        
        return $str;                                                        
    }
    
    sub load_summary_file {
        my ($self, $file_path) = @_;
        
        open (my $IN, $file_path) or
            croak( "Cannot open file: $file_path\nERROR: $!" );
        
        foreach my $line ( <$IN> ) {
            chomp $line;
            if ( $line =~ m/^Sample name: (.*)/i ) {
                $self->set_sample_name($1);
            }
            elsif ( $line =~ m/^ID: (.*)/i ) {
                $self->set_sample_id($1);
            }
            elsif ( $line =~ m/^File format: (.*)/i ) {
                $self->set_file_format($1);
            }
            elsif ( $line =~ m/^Seq count: (\d+)/i ) {
                $self->set_seq_count($1);
            }
            elsif ( $line =~ m/^Match count: (\d+)/i ) {
                $self->set_match_count($1);
            }
            elsif ( $line =~ m/^Mismatch count: (\d+)/i ) {
                $self->set_mismatch_count($1);
            }
            elsif ( $line =~ m/^Molecule Tag \(MT\) count: (\d+)/i ) {
                $self->set_MT_count($1);
            }
            elsif ( $line =~ m/^Single read category count: (\d+)/i ) {
                $self->set_SRC_count($1);
            }
        }
        
        close($IN);
        
        return 1;
    }
}


1;
__END__

#######
# POD #
#######
=head1 MTToolbox::SampleSummary

MTToolbox::SampleSummary - A class to store and print summary inforamtion about
a Sample

=head1 VERSION

This documentation refers to MTToolbox::SampleSummary version 4.1.2.

=head1 Included Modules


use Class::Std::Utils;
use Carp qw(carp croak)
use version

=head1 Inherit

NA

=head1 SYNOPSIS
    
    use MTToolbox::SampleSummary;
    my $my_sum_info = MTToolbox::SampleSummary->new();
    
    # Attribute Setters
    $my_sum_info->set_sample_name($name);
    $my_sum_info->set_sample_id($id);
    $my_sum_info->set_file_format($format);
    $my_sum_info->set_match_count($match_count);
    $my_sum_info->set_mismatch_count($mismatch_count);
    $my_sum_info->set_seq_count($seqs_count);
    $my_sum_info->set_MT_count($mt_count);
    $my_sum_info->set_SRC_count($SRC_count);  #SRC - single_read_categories
    
    # Attribute Getters
    my $file_name = $my_sum_info->get_sample_name();
    my $id = $my_sum_info->get_sample_id();
    my $format = $my_sum_info->get_file_format();
    my $match_count = $my_sum_info->get_match_count();
    my $mismatch_count = $my_sum_info->get_mismatch_count();
    my $seq_count = $my_sum_info->get_seq_count();
    my $mt_count = $my_sum_info->get_MT_count();
    my $SRC_count = $my_sum_info->get_SRC_count();
    
    # Printing Methods
    $my_sum_info->print_summary($output_file);
    $my_sum_info->to_string();
    
    # Other Methods
    $my_sum_info->load_summary_file($file);
    

=head1 DESRCIPTION

This class stores and outputs summary information about a Sample object.

=head1 METHODS

=over

    new
    set_sample_name
    set_sample_id
    set_file_format
    set_match_count
    set_mismatch_count
    set_seq_count
    set_MT_count
    set_SRC_count
    get_sample_name
    get_sample_id
    get_file_format
    get_match_count
    get_mismatch_count
    get_seq_count
    get_MT_count
    get_SRC_count
    print_summary
    to_string
    load_summary_file
    
=back

=head1 METHODS DESRCIPTION

=head2 new
    
    Title: new
    Usage: MTToolbox::SampleSummary->new();
    Function: Creates a new MTToolbox::SampleSummary object
    Returns: MTToolbox::SampleSummary
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA

=head2 set_sample_name

    Title: set_sample_name
    Usage: $my_sum_info->set_sample_name($name);
    Function: Sets the name of the sample
    Returns: 1 on successful completion
    Args: -name => The name of the sample
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 set_sample_id

    Title: set_sample_id
    Usage: $my_sum_info->set_sample_id($id);
    Function: Sets the ID of the sample
    Returns: 1 on successful completion
    Args: -id => The ID of the sample.  Usually this will have a P-number
    Throws: NA
    Comments: NA
    See Also: NA

=head2 set_file_format

    Title: set_file_format
    Usage: $my_sum_info->set_file_format($format);
    Function: Sets the file format value
    Returns: 1 on successful completion
    Args: -format => technically this should only be FASTQ at this point
    Throws: NA
    Comments: Right now these scripts only support FASTQ format.  FASTA format is
              less widely used with sequenced data now so it was not included.
              Simple scripts can generate uniform quality scores for FASTA data
              so that it can be input to these scripts.
    See Also: NA

=head2 set_match_count

    Title: set_match_count
    Usage: $my_sum_info->set_match_count($match_count);
    Function: Sets the number of match objects in this sample
    Returns: 1 on successful completion
    Args: -match_count => the number of match objects in this sample
    Throws: NA
    Comments: NA
    See Also: NA

=head2 set_mismatch_count

    Title: set_mismatch_count
    Usage: $my_sum_info->set_mismatch_count($mismatch_count);
    Function: Sets the number of mismatch objects in this sample
    Returns: 1 on successful completion
    Args: -mismatch_count => the number of mismatch objects in this sample
    Throws: NA
    Comments: NA
    See Also: NA

=head2 set_seq_count

    Title: set_seq_count
    Usage: $my_sum_info->set_seq_count($seq_count);
    Function: Sets the number of raw reads in this sample
    Returns: 1 on successful completion
    Args: -seq_count => the number of raw reads in this sample
    Throws: NA
    Comments: NA
    See Also: NA

=head2 set_MT_count

    Title: set_MT_count
    Usage: $my_sum_info->set_MT_count($mt_count);
    Function: Sets the number of unique molecule tags in this sample
    Returns: 1 on successful completion
    Args: -mt_count => the number of unique molecule tags in this sample
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 set_SRC_count

    Title: set_SRC_count
    Usage: $my_sum_info->set_SRC_count($SRC_count);
    Function: Sets the number of single read categories in this sample
    Returns: 1 on successful completion
    Args: -SRC_count => the number of single catefory reads in this sample
    Throws: NA
    Comments: A single read categories matches the expected pattern but goes
              into a molecule tag category that only has a depth of 1 (i.e.
              no other reads have that molecule tag)
    See Also: NA

=head2 get_sample_name

    Title: get_sample_name
    Usage: $my_sum_info->get_sample_name();
    Function: Gets the name of the sample
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_sample_id

    Title: get_sample_id
    Usage: $my_sum_info->get_sample_id();
    Function: Gets the ID of the sample
    Returns: String
    Args: None
    Throws: NA
    Comments: As of right now the ID is some for of a P-number (ie P1_ATCGA)
    See Also: NA

=head2 get_file_format

    Title: get_file_format
    Usage: $my_sum_info->get_file_format();
    Function: Gets the file format value
    Returns: String
    Args: None
    Throws: NA
    Comments: Technically this can only be FASTQ right now.
    See Also: NA

=head2 get_match_count

    Title: get_match_count
    Usage: $my_sum_info->get_match_count();
    Function: Gets the number of match objects in this sample
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_mismatch_count

    Title: get_mismatch_count
    Usage: $my_sum_info->get_mismatch_count();
    Function: Gets the number of mismatch objects in this sample
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_seq_count

    Title: get_seq_count
    Usage: $my_sum_info->get_seq_count();
    Function: Gets the number of raw reads in this sample
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_MT_count

    Title: get_MT_count
    Usage: $my_sum_info->get_MT_count();
    Function: Gets the number of unique molecule tags in this sample
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_SRC_count

    Title: get_SRC_count
    Usage: $my_sum_info->get_SRC_count();
    Function: Gets the number of single cateogry reads in this sample
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 print_summary

    Title: print_summary
    Usage: $my_sum_info->print_summary($file);
    Function: Prints summary informtion to $file
    Returns: 1 on successful completion
    Args: -file => the file path where the summary inforation can be printed
    Throws: NA
    Comments: Calls the to_string() method
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

=head2 load_summary_file

    Title: load_summary_file
    Usage: $my_sum_info->load_summary_file($file);
    Function: Loads the values based on those in the summary file
    Returns: 1 on successful completion
    Args: None
    Throws: NA
    Comments: This erases any values previsouly stored in this object.
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