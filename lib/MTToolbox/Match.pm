package MTToolbox::Match;


use strict;
use warnings;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(carp croak);
use BioUtils::FastqSeq 1.0.0;
use MyX::Generic 0.0.3;
use version; our $VERSION = qv('4.1.2');

{
    Readonly my $NEW_USAGE => q{ new( {id => ,
                                       desc => ,
                                       fwd_tag => ,
                                       fwd_linker => ,
                                       fwd_primer => ,
                                       fwd_amplicon => ,
                                       } ) };
    
    
    # Attributes #
    my %id_of;
    my %desc_of;
    my %fwd_tag_of;
    my %fwd_linker_of;
    my %fwd_primer_of;
    my %fwd_amplicon_of;

    # Setters #
    sub set_id;  # Remember the Id is actually the sequence id from the fastq file NOT the tag
    sub set_desc;
    sub set_fwd_tag;
    sub set_fwd_linker;
    sub set_fwd_primer;
    sub set_fwd_amplicon;
    
    # Getters #
    sub get_id;
    sub get_desc;
    sub get_header;
    sub get_fwd_tag;
    sub get_fwd_linker;
    sub get_fwd_primer;
    sub get_fwd_amplicon;
    sub get_fasta_str;           # Implement in child class
    sub get_fastq_str;           # Implement in child class
    sub get_qual_aref;           # Implement in child class
    
    # Others #
    
    
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
                               $arg_href->{fwd_tag},
                               $arg_href->{fwd_linker},
                               $arg_href->{fwd_primer},
                               $arg_href->{fwd_amplicon},
            ) {
            MyX::Generic::Undef::Param->throw(
                                              error => 'Undefined parameter value',
                                              usage => $NEW_USAGE,
                                              );
        }
        
        # Initialize the objects attributes
        $id_of{ident $new_obj} = $arg_href->{id};
        $desc_of{ident $new_obj} = $arg_href->{desc};
        $fwd_tag_of{ident $new_obj} = $arg_href->{fwd_tag};
        $fwd_linker_of{ident $new_obj} = $arg_href->{fwd_linker};
        $fwd_primer_of{ident $new_obj} = $arg_href->{fwd_primer};
        $fwd_amplicon_of{ident $new_obj} = $arg_href->{fwd_amplicon};
        
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
    
    sub set_fwd_tag {
        my ($self, $fwd_tag) = @_;
        $fwd_tag_of{ident $self} = $fwd_tag;
        return 1;
    }
    
    sub set_fwd_linker {
        my ($self, $fwd_linker) = @_;
        $fwd_linker_of{ident $self} = $fwd_linker;
        return 1;
    }
    
    sub set_fwd_primer {
        my ($self, $fwd_primer) = @_;
        $fwd_primer_of{ident $self} = $fwd_primer;
        return 1;
    }
    
    sub set_fwd_amplicon {
        my ($self, $fwd_amplicon) = @_;
        $fwd_amplicon_of{ident $self} = $fwd_amplicon;
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
    
    sub get_fwd_tag {
        my ($self) = @_;
        return $fwd_tag_of{ident $self};
    }
    
    sub get_fwd_linker {
        my ($self) = @_;
        return  $fwd_linker_of{ident $self};
    }
    
    sub get_fwd_primer {
        my ($self) = @_;
        return $fwd_primer_of{ident $self};
    }
    
    sub get_fwd_amplicon {
        my ($self) = @_;
        return $fwd_amplicon_of{ident $self};
    }
    
    # Implement these in child classes
    sub get_fasta_str { return -1; }
    sub get_fastq_str { return -1; }
    sub get_qual_aref { return -1; }
    
    
    ##########
    # Others #
    ##########
    
}

1;
__END__

#######
# POD #
#######
=head1 MTToolbox::Match

MTToolbox::Match - A class to store a sequence that has defined attributes
                   (eg. primer, amplicon, etc)

=head1 VERSION

This documentation refers to MTToolbox::Match version 4.1.2.

=head1 Included Modules

    Class::Std::Utils
    List::MoreUtils qw(any)
    Readonly
    Carp qw(croak carp)
    BioUtils::FastqSeq 1.0.0
    MyX::Generic 1.0.0
    version

=head1 Inherit

NA

=head1 SYNOPSIS
    
    use MTToolbxo::Match;
    my $my_match = MTToolbox::Match->new({ id => $id,
                        desc => $desc,
                        fwd_tag => $fwd_tag,
                        fwd_linker => $fwd_linker,
                        fwd_primer => $fwd_primer,
                        fwd_amplicon => $fwd_amplicon,
                    });
    
    # Attribute Getters
    my $id = $my_match->get_id();
    my $desc = $my_match->get_desc();
    my $header = $my_match->get_header();
    my $fwd_tag = $my_match->get_fwd_tag();
    my $fwd_linker = $my_match->get_fwd_linker();
    my $fwd_primer = $my_match->get_fwd_primer();
    my $fwd_amplicon = $my_match->get_fwd_amplicon();
    
    # Attribute Setters
    $my_match->set_id($id);
    $my_match->set_desc($desc);
    $my_match->set_fwd_tag($fwd_tag);
    $my_match->set_fwd_linker($fwd_linker);
    $my_match->set_fwd_primer($fwd_primer);
    $my_match->set_fwd_amplicon($fwd_amplicon);
    
    # Other Getter Methods
    $my_match->get_custalw_id();
    $my_match->get_fasta_str();            # Implemented in child class
    $my_match->get_fastq_str();            # Implemented in child class
    $my_match->get_qual_aref();            # Implemented in child class
    
    # Other Methods
    

=head1 DESCRIPTION

MTToolbox::Match is primarily a data structure to store the attributes of a
sequences that matches a defined pattern.  That pattern is defined in the Sample
object and the various child classes of Sample.  MTToolbox::Match is like an
abstract class.  The attributes that are stored in MTToolbox::Match, which will
be in all MTToolbox::Match objects, include id, fwd_tag, fwd_linker, fwd_primer,
and fwd_amplicon.

When clustalw is run (either from the commandline or from bioperl) it requires
that the names of the sequences are <= 30 characters long.  Originally the
headers where modified based on the given Illumina headers.  Now each sequence
is given a unique ID. This ID is the sample ID (ie P1) followed by an underscore
character followed by the current count of sequences seen in P1.  For example,
P1_1202 represents the 1202 sequence in sample P1.  

Another method to note is get_qual_aref.  This method is implemented in the
various child classes of MTToolbox::Match.  It is important because when I build
a consensus sequence from a set of MTToolbox::Match objects I need to access
their quality values.  Clustalw doesn't do anything with quality scores.  So
when I run Clustalw I only use the raw sequences.  Then when I make the
consensus I access the quality values by calling get_qual_aref.  I need the
quality values when making a consensus to resolve ties between two or more
aligned bases.

There are also various methods for obtainging the objects information in common
string formats (ie. fastq).  

=head1 METHODS

=over

    new
    set_id
    set_desc
    set_fwd_tag
    set_fwd_linker
    set_fwd_primer
    set_fwd_amplicon
    get_id
    get_desc
    get_header
    get_fwd_tag
    get_fwd_linker
    get_fwd_primer
    get_fwd_amplicon
    get_fasta_str           # Implemented in child class
    get_fastq_str           # Implemented in child class
    get_qual_aref           # Implemented in child class
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: MTToolbox::Match->new({
                id => $id,
                desc => $desc,
                fwd_tag => $fwd_tag,
                fwd_linker => $fwd_linker,
                fwd_primer => $fwd_primer,
                fwd_amplicon => $fwd_amplicon,
            });
    Function: Creates a new MTToolbox::Match object
    Returns: MTToolbox::Match
    Args: -id           => a string representing sequence id (NOT the molecule tag)
          -desc         => a string representing the sequence description
          -fwd_tag      => a FastqSeq representing the molecule tag
          -fwd_linker   => a FastqSeq representing the fwd_linker
          -fwd_primer   => a FastqSeq representing the forward primer
          -fwd_amplicon => a FastqSeq representing the forward amplicon
    Throws: MyX::Generic::Undef::Param
    Comments: NA
    See Also: BioUtils::FastqSeq

=head2 set_id

    Title: set_id
    Usage: $my_match->set_id($id);
    Function: Sets the id value
    Returns: 1 on successful completion
    Args: -id => a string representing the raw sequence id
    Throws: NA
    Comments: This NOT the molecule tag id
    See Also: NA

=head2 set_desc

    Title: set_desc
    Usage: $my_match->set_desc($id);
    Function: Sets the description
    Returns: 1 on successful completion
    Args: -id => a string representing the sequence description
    Throws: NA
    Comments: The description can be anything.  When printing the sequences in
              common formats the description is followed by the ID and is
              seperate by a space.
    See Also: NA

=head2 set_fwd_tag

    Title: set_fwd_tag
    Usage: $my_match->set_fwd_tag($fwd_tag);
    Function: Sets the fwd_tag value
    Returns: 1 on successful completion
    Args: -fwd_tag => a FastqSeq representing the forward tag
    Throws: NA
    Comments: NA
    See Also: BioUtils::FastqSeq

=head2 set_fwd_linker

    Title: set_fwd_linker
    Usage: $my_match->set_fwd_linker($fwd_linker);
    Function: Sets the fwd_linker value
    Returns: 1 on successful completion
    Args: -fwd_linker => a FastqSeq representing the forward linker
    Throws: NA
    Comments: NA
    See Also: BioUtils::FastqSeq

=head2 set_fwd_primer

    Title: set_fwd_primer
    Usage: $my_match->set_fwd_primer($fwd_primer);
    Function: Sets the fwd_primer value
    Returns: 1 on successful completion
    Args: -fwd_primer => a FastqSeq representing the forward primer
    Throws: NA
    Comments: NA
    See Also: BioUtils::FastqSeq

=head2 set_fwd_amplicon

    Title: set_fwd_amplicon
    Usage: $my_match->set_fwd_amplicon($fwd_amplicon);
    Function: Sets the fwd_amplicon value
    Returns: 1 on successful completion
    Args: -fwd_amplicon => a FastqSeq representing the forward amplicon
    Throws: NA
    Comments: NA
    See Also: BioUtils::FastqSeq
    
=head2 get_id
    
    Title: get_id
    Usage: my $id = $my_match->get_id();
    Function: Gets the id
    Returns: String
    Args: None
    Throws: MyX::Generic::Undef::Attribute
    Comments: This is the raw sequence ID and NOT the molecule tag ID
    See Also: NA
    
=head2 get_desc
    
    Title: get_desc
    Usage: my $id = $my_match->get_desc();
    Function: Gets the description
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_header
    
    Title: get_header
    Usage: my $id = $my_match->get_header();
    Function: Gets the header
    Returns: String
    Args: None
    Throws: NA
    Comments: The header is not automatically preceded by a '@' or '>'
              character.  The header consists of the id followed by a space
              followed by the description (if one exists).
    See Also: NA
    
=head2 get_fwd_tag
    
    Title: get_fwd_tag
    Usage: my $id = $my_match->get_fwd_tag();
    Function: Gets the forward tag
    Returns: BioUtils::FastqSeq
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_fwd_linker
    
    Title: get_fwd_linker
    Usage: my $id = $my_match->get_fwd_linker();
    Function: Gets the forward linker
    Returns: BioUtils::FastqSeq
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_fwd_primer
    
    Title: get_fwd_primer
    Usage: my $id = $my_match->get_fwd_primer();
    Function: Gets the forward primer
    Returns: BioUtils::FastqSeq
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_fwd_amplicon
    
    Title: get_fwd_amplicon
    Usage: my $id = $my_match->get_fwd_amplicon();
    Function: Gets the forward amplicon
    Returns: BioUtils::FastqSeq
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_fasta_str
    
    Title: get_fasta_str
    Usage: my $fasta_str = $my_match->get_fasta_str();
    Function: Get a string that can be printed to a file and subsiquently input
              into clustalw
    Returns: -1
    Args: None
    Throws: NA
    Comments: This method is implemented in the child classes of
              MTToolbox::Match. This method concatenates the primer and amplicon
              attribute sequences that this object stores.  To get additional or
              fewer parts of the sequence you can change the code of this method.
    See Also: NA
    
=head2 get_fastq_str
    
    Title: get_fastq_str
    Usage: my $fastqStr = $my_match->get_fastq_str();
    Function: Get a fastq formatted string that can be printed to a file
    Returns: -1
    Args: None
    Throws: NA
    Comments: This method is implemented in the child classes of
              MTToolbox::Match. This method concatenates the primer and amplicon
              attribute sequences and quality values that this object stores.
              To get additional or fewer parts of the sequence you can change
              the code of this method.
    See Also: NA
    
=head2 get_qual_aref
    
    Title: get_qual_aref
    Usage: my $quals_aref = $my_match->get_qual_aref();
    Function: Gets the quality values of the primer and amplicon attributes
    Returns: -1
    Args: None
    Throws: NA
    Comments: This method is implemented in the child classes of
              MTToolbox::Match. It concatenates the primer and amplicon
              attribute quality values that this object stores.  To get
              additional or fewer parts of the sequence you can change the code
              of this method.  Athough if you change those parts you will also
              want to change get_fastq_str and get_fasta_str.
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

