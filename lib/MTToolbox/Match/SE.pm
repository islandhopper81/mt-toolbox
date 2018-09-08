package MTToolbox::Match::SE;


use strict;
use warnings;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(croak);
use BioUtils::FastqSeq 1.0.0;
use MyX::Generic 1.0.0;
use version; our $VERSION = qv('4.1.2');

use MTToolbox::Match;
use base qw(MTToolbox::Match);

{
    Readonly my $NEW_USAGE => q{ new( {id => ,
                                       desc =>,
                                       fwd_tag => ,
                                       fwd_linker => ,
                                       fwd_primer => ,
                                       fwd_amplicon => ,
                                       fwd_tail => ,
                                       } ) };

    # Attributes #
    my %fwd_tail_of;
    
    # Setters #
    sub set_fwd_tail;
    
    # Getters #
    sub get_fwd_tail;
    sub get_fasta_str;
    sub get_fastq_str;
    sub get_quals_aref;
    
    
    # Others #
    
    
    ###############
    # Constructor #
    ###############
    sub new {
        my ($class, $arg_href) = @_;
        
        # Croak if calling 'new' on an already blessed reference
        croak 'Constructor called on existing object instead of class'
            if ref $class;
            
        # Make sure the required parameters are defined
        if ( any {!defined $_} $arg_href->{id},
                               $arg_href->{desc},
                               $arg_href->{fwd_tag},
                               $arg_href->{fwd_linker},
                               $arg_href->{fwd_primer},
                               $arg_href->{fwd_amplicon},
                               $arg_href->{fwd_tail},
            ) {
            MyX::Generic::Undef::Param->throw(
                                              error => 'Undefined parameter value',
                                              usage => $NEW_USAGE,
                                              );
        }
        
        # Bless a scalar to instantiate an object
        my $new_obj = $class->SUPER::new($arg_href);
        my $ident = ident($new_obj);
        
        $fwd_tail_of{ident $new_obj} = $arg_href->{fwd_tail};
        
        return $new_obj;
    }
    
    
    ###########
    # Setters #
    ###########
    sub set_fwd_tail {
        my ($self, $fwd_tail) = @_;
        $fwd_tail_of{ident $self} = $fwd_tail;
        return 1;
    }
    
    
    
    ###########
    # Getters #
    ###########
    sub get_fwd_tail {
        my ($self) = @_;
        return $fwd_tail_of{ident $self};
    }
    
    sub get_fasta_str {
        my ($self) = @_;
        
        my $str .= ">" . $self->get_header() . "\n";
        $str .= $self->get_fwd_amplicon()->get_seq() . "\n";
        
        return $str;
    }
    
    sub get_fastq_str {
        my ($self) = @_;
        
        my $str .= "@" . $self->get_header() . "\n";
        $str .= $self->get_fwd_amplicon()->get_seq() . "\n";
        $str .= "+\n";
        $str .= $self->get_fwd_amplicon()->get_quals_str() . "\n";
        
        return $str;
    }
    
    sub get_quals_aref {
        my ($self) = @_;
        
        my @quals = (@{$self->get_fwd_amplicon()->get_quals_aref()});
        
        return \@quals;
    }
    
    
    
    ##########
    # Others #
    ##########
}

1;
__END__

#######
# POD #
#######
=head1 MTToolbox::Match::SE

MTToolbox::Match::SE - A class to store a single end (SE) sequence that has
                       defined attributes (eg. primer, amplicon, etc).

=head1 VERSION

This documentation refers to MTToolbox::Match::SE version 4.1.2.

=head1 Included Modules

    Class::Std::Utils
    List::MoreUtils qw(any)
    Readonly
    Carp qw(carp croak)
    FastqSeq 1.0.0
    MyX::Generic 1.0.0
    version
    MTToolbox::Match

=head1 Inherit

    MTToolbox::Match

=head1 SYNOPSIS
    
    use MTToolbox::Match::SE;
    my $my_match = MTToolbox::Match::SE->new({
                        id => $id,
                        desc => $desc,
                        fwd_tag => $fwd_tag,
                        fwd_linker => $fwd_linker,
                        fwd_primer => $fwd_primer,
                        fwd_amplicon => $fwd_amplicon,
                        fwd_tail => $fwd_tail
                    });
    
    # Attribute Getters -- see also parent classes MTToolbox::Match
    # and MTToolbox::Match::PE
    my $fwd_tail = $my_match->get_fwd_tail();
    
    # Attribute Setters -- see also parent classes MTToolbox::Match
    # and MTToolbox::Match::PE
    $my_match->set_fwd_tail($fwd_tail);
    
    # Other Getter Methods implemented in this class
    $my_match->get_fasta_str();
    $my_match->get_fastq_str();
    $my_match->get_quals_aref();
    

=head1 DESCRIPTION

MTToolbox::Match::SE is primarily a data structure to store the attributes of
single end sequences that match a defined pattern.  That pattern is defined in
the Sample object and the various child classes of Sample.  MTToolbox::Match::SE
inherits from the parent class MTToolbox::Match.

A method to note is get_quals_aref.  It is important because when I build a
consensus sequence from a set of MTToolbox::Match objects I need to access their
quality values.  Clustalw doesn't do anything with quality scores.  So when I
run Clustalw I only use the raw sequences.  Then when I make the consensus I
access the quality values by calling get_quals_aref.  I need the quality values
when making a consensus to resolve ties between two or more aligned bases.

There are also various methods for obtainging the objects information in common
string formats (ie. fastq).  

=head1 METHODS

=over

    new
    set_fwd_tail
    get_fwd_tail
    get_fasta_str
    get_fastq_str
    get_quals_aref
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: MTToolbox::Match->new({ id => $id,
                desc => $desc,
                fwd_tag => $fwd_tag,
                fwd_linker => $fwd_linker,
                fwd_primer => $fwd_primer,
                fwd_amplicon => $fwd_amplicon,
                fwd_tail => $fwd_tail
            });
    Function: Creates a new MTToolbox::Match::SE object
    Returns: MTToolbox::Match::SE
    Args: -id           => a string representing sequence id (NOT the molecule tag)
          -desc         => a string representing the sequence description
          -fwd_tag      => a FastqSeq representing the molecule tag
          -fwd_linker   => a FastqSeq representing the fwd_linker
          -fwd_primer   => a FastqSeq representing the forward primer
          -fwd_amplicon => a FastqSeq representing the forward amplicon
          -fwd_tail     => a FastqSeq representing the forward tail
    Throws: MyX::Generic::Undef::Param
    Comments: NA
    See Also: BioUtils::BioUtils::FastqSeq
 
=head2 set_fwd_tail

    Title: set_fwd_tail
    Usage: $my_match->set_fwd_tail($fwd_tail);
    Function: Sets the fwd_tag value
    Returns: 1 on successful completion
    Args: -fwd_tail => a FastqSeq representing the forward tail
    Throws: NA
    Comments: A tail occurs in the event of an intended frame shift in the
              molecule tag.  It is advantageous in future clustering steps to
              have the sequences be the same length.  If there is no frame shift
              in the molecule tag the tail consists of two extra bases that are
              not included in the final sequence output and analysis.
    See Also: BioUtils::FastqSeq
    
=head2 get_fwd_tail
    
    Title: get_fwd_tail
    Usage: my $id = $my_match->get_fwd_tail();
    Function: Gets the forward tail
    Returns: BioUtils::FastqSeq
    Args: None
    Throws: NA
    Comments: A tail occurs in the event of an intended frame shift in the
              molecule tag.  It is advantageous in future clustering steps to
              have the sequences be the same length.  If there is no frame shift
              in the molecule tag the tail consists of two extra bases that are
              not included in the final sequence output and analysis.
    See Also: NA
   
=head2 get_fasta_str
    
    Title: get_fasta_str
    Usage: my $fasta_str = $my_match->get_fasta_str();
    Function: Get a string that can be printed to a file an subsiquently input
              into clustalw
    Returns: String
    Args: None
    Throws: NA
    Comments: This method concatenates only the amplicon attribute
              sequence that this object stores.  To get additional or fewer
              parts of the sequence you can change the code of this method.  If
              the attributes concatenated and returned by this method are
              changed, they should also be changed in get_fastq_str and
              get_quals_aref.
    See Also: get_fastq_str, get_quals_aref
    
=head2 get_fastq_str
    
    Title: get_fastq_str
    Usage: my $fastqStr = $my_match->get_fastq_str();
    Function: Get a fastq formatted string that can be printed to a file
    Returns: String
    Args: None
    Throws: NA
    Comments: This method concatenates only the amplicon attribute
              sequence and quality values that this object stores.  To get
              additional or fewer parts of the sequence you can change the code
              of this method.  If the attributes concatenated and returned by
              this method are changed, they should also be changed in
              get_fasta_str and get_quals_aref.
    See Also: get_fasta_str, get_quals_aref
    
=head2 get_quals_aref
    
    Title: get_quals_aref
    Usage: my $quals_aref = $my_match->get_quals_aref();
    Function: Gets the quality values of the primer and amplicon attributes
    Returns: Array reference
    Args: None
    Throws: NA
    Comments: This method is implemented in the child classes of MTToolbox::Match.
              It concatenates only the amplicon attribute quality values that
              this object stores.  To get additional or fewer parts of the
              sequence you can change the code of this method.  Athough if you
              change those parts you will also want to change get_fastq_str and
              get_fasta_str.
    See Also: get_fastq_str, get_fasta_str
    

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