package MTToolbox::SampleSummary::PE::Overlap;

# This object is a data structure that holds summary information for
# a sample with overlapping reads. 


use strict;
use warnings;

use Class::Std::Utils;
use Carp qw(carp croak);
use version; our $VERSION = qv('4.1.2');
use MTToolbox::SampleSummary;
use base qw(MTToolbox::SampleSummary);  # inherits from MTToolbox::SampleSummary


{
    
    # Attributes #
    my %merged_count_of;
    my %not_merged_count_of;

    # Setters #
    sub set_merged_count;
    sub set_not_merged_count;
    
    # Getters #
    sub get_merged_count;
    sub get_not_merged_count;
    
    # Others #
    sub to_string;
    
    
    ###############
    # Constructor #
    ###############
    sub new {
        my ($class) = @_;

        # Croak if calling 'new' on an already blessed reference
        croak 'Constructor called on existing object instead of class'
            if ref $class;
        
        # Bless a scalar to instantiate an object
        my $new_obj = $class->SUPER::new();
        my $ident = ident($new_obj);
        
        # Initialize the objects attributes        
        
        return $new_obj;
    }
    
    
    ###########
    # Setters #
    ###########
    sub set_merged_count {
        my ($self, $m) = @_;
        $merged_count_of{ident $self} = $m;
        return 1;
    }
    
    sub set_not_merged_count {
        my ($self, $n) = @_;
        $not_merged_count_of{ident $self} = $n;
        return 1;
    }
    
    
    ###########
    # Getters #
    ###########
    sub get_merged_count {
        my ($self) = @_;
        
        my $merged_count = $merged_count_of{ident $self};
        if ( defined $merged_count ) {
            return $merged_count;
        }
        
        return "--";
    }
    
    sub get_not_merged_count {
        my ($self) = @_;
        
        my $not_merged_count = $not_merged_count_of{ident $self};
        if ( defined $not_merged_count ) {
            return $not_merged_count;
        }
        
        return "--";
    }
    
    sub get_perc_merged {
        my ($self) = @_;
        
        my $merged_count = $merged_count_of{ident $self};
        my $total_seqs = $self->get_seq_count();
        if ( defined $merged_count ) {
            return int($merged_count / $total_seqs * 100) . "%";
        }
        
        return "--";
    }
    
    
    
    ##########
    # Others #
    ##########
    sub to_string {
        my ($self) = @_;
        
        my $str = "Sample type: PE w/ Overlap\n";
        $str .= "Sample name: " . $self->get_sample_name() . "\n";
        $str .= "ID: " . $self->get_sample_id() . "\n";
        $str .= "File format: " . $self->get_file_format() . "\n";
        $str .= "Seq count: " . $self->get_seq_count() . "\n";
        $str .= "Merged count: " . $self->get_merged_count() . "\n";
        $str .= "Not merged count: " . $self->get_not_merged_count() . "\n";
        $str .= "Percent merged: " . $self->get_perc_merged() . "\n";
        $str .= "Match count: " . $self->get_match_count() . "\n";
        $str .= "Mismatch count: " . $self->get_mismatch_count() . "\n";
        $str .= "Molecule Tag (MT) count: " . $self->get_MT_count() . "\n";
        $str .= "Single read category count: " . $self->get_SRC_count() . "\n";
        
        return $str;                                                        
    }
    
    sub load_summary_file {
        my ($self, $file_path) = @_;
        
        open (my $IN, $file_path) or
            croak( "Cannot open file: $file_path" );
        
        foreach my $line ( <$IN> ) {
            chomp $line;
            if ( $line =~ m/^Sample name: (.*)/i ) {
                $self->set_sample_name($1);
            }
            elsif ( $line =~ m/ID: (.*)/i ) {
                $self->set_sample_id($1);
            }
            elsif ( $line =~ m/^File format: (.*)/i ) {
                $self->set_file_format($1);
            }
            elsif ( $line =~ m/^Seq count: (\d+)/i ) {
                $self->set_seq_count($1);
            }
            elsif ( $line =~ m/^Merged count: (\d+)/i ) {
                $self->set_merged_count($1);
            }
            elsif ( $line =~ m/^Not merged count: (\d+)/i ) {
                $self->set_not_merged_count($1);
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
=head1 MTToolbox::SampleSummary::PE::Overlap

MTToolbox::SampleSummary::PE::Overlap - A class to store and print summary
inforamtion about a Sample with overlapping paired-end reads

=head1 VERSION

This documentation refers to MTToolbox::SampleSummary::PE::Overlap version 4.1.2.

=head1 Included Modules

Class::Std::Utils;
Carp qw(carp croak)
version
MTToolbox::SampleSummary

=head1 Inherit

MTToolbox::SampleSummary

=head1 SYNOPSIS
    
    use MTToolbox::SampleSummary::PE::Overlap;
    my $my_sum_info = MTToolbox::SampleSummary::PE::Overlap->new();
    
    # Attribute Setters -- See also MTToolbox::SampleSummary
    $my_sum_info->set_merged_count($merged_count);
    $my_sum_info->set_not_merged_count($not_merged_count);
    
    # Attribute Getters -- See also MTToolbox::SampleSummary
    my $merged_count = $my_sum_info->get_merged_count();
    my $not_merged_count = $my_sum_info->get_not_merged_count();
    my $perc_merged = $my_sum_info->get_perc_merged();
    
    # Printing Methods -- See also MTToolbox::SampleSummary
    $my_sum_info->to_string();
    

=head1 DESCRIPTION

This class stores and outputs summary information about a Sample object
with PE reads.

=head1 METHODS

=over

    new
    set_merged_count
    set_not_merged_count
    get_merged_count
    get_not_merged_count
    get_perc_merged
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: MTToolbox::SampleSummary::PE::Overlap->new();
    Function: Creates a new MTToolbox::SampleSummary::PE::Overlap object
    Returns: MTToolbox::SampleSummary::PE::Overlap
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA

=head2 set_merged_count

    Title: set_merged_count
    Usage: $my_sum_info->set_merged_count($merged_count);
    Function: Sets the number of reads merged by FLASH
    Returns: 1 on successful completion
    Args: -merged_count => the number of merged reads
    Throws: NA
    Comments: NA
    See Also: NA

=head2 set_not_merged_count

    Title: set_not_merged_count
    Usage: $my_sum_info->set_not_merged_count($not_merged_count);
    Function: Sets the number of reads NOT merged by FLASH
    Returns: 1 on successful completion
    Args: -not_merged_count => the number of reads that failed to merge
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_merged_count

    Title: get_merged_count
    Usage: $my_sum_info->get_merged_count();
    Function: Gets the number of reads that were merged by FLASH
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_not_merged_count

    Title: get_not_merged_count
    Usage: $my_sum_info->get_not_merged_count();
    Function: Gets the number of reads that were NOT merged by FLASH
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_perc_merged

    Title: get_perc_merged
    Usage: $my_sum_info->get_perc_merged();
    Function: Gets the percentage of reads merged by FLASH
    Returns: String
    Args: None
    Throws: NA
    Comments: The return value is an int with the '%'
    See Also: NA
    
=head2 to_string

    Title: to_string
    Usage: $my_sum_info->to_string();
    Function: Gets the summary information stored in this object as a string
    Returns: String
    Args: None
    Throws: NA
    Comments: Overloaded from parent
    See Also: NA

=head2 load_summary_file

    Title: load_summary_file
    Usage: $my_sum_info->load_summary_file($file);
    Function: Loads the values based on those in the summary file
    Returns: 1 on successful completion
    Args: None
    Throws: NA
    Comments: This erases any values previsouly stored in this object.
              Overloaded from parent.
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