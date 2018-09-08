package MTToolbox::Mismatch;


use strict;
use warnings;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(croak);
use MyX::Generic 0.0.3;
use version; our $VERSION = qv('4.1.2');

{
    Readonly my $NEW_USAGE => q{ new( {id => ,
                                       desc => ,
                                       seq => ,
                                       quals_str => ,
                                       } ) };

    # Attributes #
    my %id_of;
    my %desc_of;
    my %seq_of;
    my %quals_str_of;
    
    # Setters #
    sub set_id;  # Id is the raw read id from the fastq file NOT the molecule tag
    sub set_desc;
    sub set_seq;
    sub set_quals_str;
    
    # Getters #
    sub get_id;
    sub get_desc;
    sub get_header;
    sub get_seq;
    sub get_quals_str;
    
    
    ##########
    # Others #
    ##########
    sub get_fastq_str;
    
    
    ###############
    # Constructor #
    ###############
    sub new {
        my ($class, $arg_href) = @_;
        
        # Croak if calling 'new' on an already blessed reference
        croak 'Constructor called on existing object instead of class'
            if ref $class;
        
        # Bless a scalar to instantiate an object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        # Make sure the required parameters are defined
        if ( any {!defined $_} $arg_href->{id},
                               $arg_href->{desc},
                               $arg_href->{seq},
                               $arg_href->{quals_str},
            ) {
            MyX::Generic::Undef::Param->throw(
                                              error => 'Undefined parameter value',
                                              usage => $NEW_USAGE,
                                              );
        }
        
        # Initialize the objects attributes
        $id_of{ident $new_obj} = $arg_href->{id};
        $desc_of{ident $new_obj} = $arg_href->{desc};
        $seq_of{ident $new_obj} = $arg_href->{seq};
        $quals_str_of{ident $new_obj} = $arg_href->{quals_str};
    
        return $new_obj;
    }
    
    
    ###########
    # Setters #
    ###########
    sub set_id {
        my ($self, $id) = @_;
        $id_of{ident $self} = $id;
        return 1;
    }
    
    sub set_desc {
        my ($self, $desc) = @_;
        $desc_of{ident $self} = $desc;
        return 1;
    }
    
    sub set_seq {
        my ($self, $seq) = @_;
        $seq_of{ident $self} = $seq;
        return 1;
    }
    
    sub set_quals_str {
        my ($self, $quals_str) = @_;
        $quals_str_of{ident $self} = $quals_str;
        return 1;
    }
    
    
    
    ###########
    # Getters #
    ###########
    sub get_id {
        my ($self) = @_;
        
        # check if the header is defined
        if ( ! defined $id_of{ident $self} ) {
            MyX::Generic::Undef::Attribute->throw(
                                                  error => 'Undefined ID'
                                                  );
        }
        
        return $id_of{ident $self};
    }
    
    sub get_desc {
        my ($self) = @_;
        
        if ( my $_desc = $desc_of{ident $self} ) {
            return $_desc;
        }
        
        return "";
    }
    
    sub get_header {
        my ($self) = @_;
        
        # Note this doesn't not attach a '@' or '>' to the beginning
        
        # only add the desc if there is something there
        if ( my $_desc = $self->get_desc() ) {
            return $self->get_id() . " " . $_desc;
        }
        
        return $self->get_id();
    }
    
    sub get_seq {
        my ($self) = @_;
        return $seq_of{ident $self};
    }
    
    sub get_quals_str {
        my ($self) = @_;
        return $quals_str_of{ident $self};
    }
    
    
    
    
    ##########
    # Others #
    ##########
    sub get_fastq_str {
        my ($self) = @_;
        
        my $str = "@" . $self->get_header() . "\n";
        $str .= $self->get_seq() . "\n";
        $str .= "+\n";
        $str .= $self->get_quals_str() . "\n";
        
        return $str;
    }
}


1;
__END__

#######
# POD #
#######
=head1 MTToolbox::Mismatch

MTToolbox::Mismatch - A data structure for storing information about raw reads
that don't match the expected sequence pattern defined in the _parseRead method
of one of the Sample classes

=head1 VERSION

This documentation refers to MTToolbox::Mismatch version 4.1.2.

=head1 Included Modules

    Class::Std::Utils
    List::MoreUtils qw(any)
    Readonly
    Carp qw(croak)
    MyX::Generic 1.0.0
    version

=head1 Inherit

    NA

=head1 SYNOPSIS
    
    use MTToolbox::Mismatch;
    my $my_mismatch = MTToolbox::Mismatch->new({
                            id => $id,
                            desc => $desc,
                            seq => $seq,
                            quals_str => $quals_str
                        });
    
    # Attribute Getters
    my $id = $my_mismatch->get_id();
    my $desc = $my_mismatch->get_desc();
    my $header = $my_mismatch->get_header();
    my $seq = $my_mismatch->get_seq();
    my $quals_str = $my_mismatch->get_quals_strStr();
    
    # Attribute Setters
    $my_mismatch->set_id("seq1");
    $my_mismatch->set_desc("arabidopsis root microbime sequence");
    $my_mismatch->set_seq("ATCG");
    $my_mismatch->set_quals_strStr('9999');
    
    # Others
    my $fastq_str = $my_mismatch->get_fastq_str();
    

=head1 DESCRIPTION

MTToolbox::Mismatch is a data structure for storing information about raw reads that don't
match the expected sequence pattern defined in the _parseRead method of one of
the Sample classes.  

=head1 METHODS

=over

    new
    set_id
    set_desc
    set_seq
    set_quals_str
    get_id
    get_desc
    get_header
    get_seq
    get_quals_str
    get_fastq_str
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: MTToolbox::Mismatch->new({
                id => $id,
                desc => $desc,
                seq => $seq,
                quals_str => $quals_str
            });
    Function: Creates a new MTToolbox::Mismatch object
    Returns: MTToolbox::Mismatch
    Args: -id => the read ID from the original FASTQ file
          -desc => a string representing the sequence description
          -seq => the sequence string
          -quals_str => the quality values as a string encoded in Illumina-1.8+
    Throws: MyX::Generic::Undef::Param
    Comments: NA
    See Also: NA

=head2 set_id

    Title: set_id
    Usage: $my_mismatch->set_id($id);
    Function: Sets the id value
    Returns: 1 on successful completion
    Args: -id => a sequence id
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 set_desc

    Title: set_desc
    Usage: $my_mismatch->set_desc($id);
    Function: Sets the description
    Returns: 1 on successful completion
    Args: -id => a string representing the sequence description
    Throws: NA
    Comments: The description can be anything.  When printing the sequences in
              common formats the description is followed by the ID and is
              seperate by a space.
    See Also: NA

=head2 set_seq

    Title: set_seq
    Usage: $my_mismatch->set_seq($seq);
    Function: Sets the seq value
    Returns: 1 on successful completion
    Args: -seq => a seqeunce string
    Throws: NA
    Comments: NA
    See Also: NA

=head2 set_quals_str

    Title: set_quals_str
    Usage: $my_mismatch->set_quals_str($quals_str);
    Function: Sets the quals_str value
    Returns: 1 on successful completion
    Args: -quals_str => a string of quality scores in Illumina-1.8+ encoding
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_id

    Title: get_id
    Usage: $my_mismatch->get_id();
    Function: Gets the ID
    Returns: String
    Args: None
    Throws: MyX::Generic::Undef::Attribute
    Comments: NA
    See Also: NA

=head2 get_desc
    
    Title: get_desc
    Usage: my $id = $my_mismatch->get_desc();
    Function: Gets the description
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_header
    
    Title: get_header
    Usage: my $id = $my_mismatch->get_header();
    Function: Gets the header
    Returns: String
    Args: None
    Throws: NA
    Comments: The header is not automatically preceded by a '@' or '>'
              character.  The header consists of the id followed by a space
              followed by the description (if one exists).
    See Also: NA

=head2 get_seq

    Title: get_seq
    Usage: $my_mismatch->get_seq();
    Function: Gets the sequence
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_quals_str

    Title: get_quals_str
    Usage: $my_mismatch->get_quals_str();
    Function: Gets the quality values as a string in Illumina-1.8+ encoding
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_fastq_str

    Title: get_fastq_str
    Usage: $my_mismatch->get_fastq_str();
    Function: Gets the mismatch sequence in FASTQ format
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
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
