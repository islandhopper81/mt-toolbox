package MTToolbox::Match::PE::NoOverlap;


use strict;
use warnings;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(carp croak);
use BioUtils::FastqSeq 1.0.0;
use MyX::Generic 0.0.3;
use version; our $VERSION = qv('4.1.2');

use MTToolbox::Match::PE;
use base qw(MTToolbox::Match::PE);  # inherits from PE (which also inherits from Match)


{
    Readonly my $NEW_USAGE => q{ new( {id => ,
                                       desc => ,
                                       fwd_tag => ,
                                       fwd_linker => ,
                                       fwd_primer => ,
                                       fwd_amplicon => ,
                                       fwd_tail => ,
                                       rev_tag => ,
                                       rev_linker => ,
                                       rev_primer => ,
                                       rev_amplicon => ,
                                       rev_tail => ,
                                       } ) }
                                       ;
    # Attributes #
    my %fwd_tail_of;
    my %rev_amplicon_of;
    my %rev_tail_of;
    
    # Setters #
    sub set_fwd_tail;
    sub set_rev_amplicon;
    sub set_rev_tail;
    
    # Getters #
    sub get_fwd_tail;
    sub get_rev_amplicon;
    sub get_rev_tail;
    sub get_fasta_str;
    sub get_fwd_fasta_str;
    sub get_rev_fasta_str;
    sub get_fastq_str;
    sub get_fwd_fastq_str;
    sub get_rev_fastq_str;
    sub get_quals_aref;
    sub get_fwd_quals_aref;
    sub get_rev_quals_aref;
    
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
                               $arg_href->{rev_tag},
                               $arg_href->{rev_primer},
                               $arg_href->{rev_linker},
                               $arg_href->{rev_amplicon},
                               $arg_href->{rev_tail},
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
        $rev_amplicon_of{$ident} = $arg_href->{rev_amplicon};
        $rev_tail_of{$ident} = $arg_href->{rev_tail};
        
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
    
    sub set_rev_amplicon {
        my ($self, $rev_amplicon) = @_;
        $rev_amplicon_of{ident $self} = $rev_amplicon;
        return 1;
    }
    
    sub set_rev_tail {
        my ($self, $rev_tail) = @_;
        $rev_tail_of{ident $self} = $rev_tail;
        return 1;
    }
    
    
    ###########
    # Getters #
    ###########
    sub get_fwd_tail {
        my ($self) = @_;
        return $fwd_tail_of{ident $self};
    }
    
    sub get_rev_amplicon() {
        my ($self) = @_;
        return $rev_amplicon_of{ident $self};
    }
    
    sub get_rev_tail() {
        my ($self) = @_;
        return $rev_tail_of{ident $self};
    }
    
    sub get_fasta_str {
        my ($self) = @_;
        
        # IMPORTANT: I should probably not use this method.  The main reason why
        # this method is not useful is that the id's on the fwd and rev seqeunces
        # are not unique.  At one point I tried to change this by putting "_fwd" and
        # "_rev on the end of the ids in get__fwd_fasta_str but that caused
        # problems in other sections of the code.
        
        my $str .= $self->get_fwd_fasta_str();
        $str .= $self->get_rev_fasta_str();
        
        return $str;
    }
    
    sub get_fwd_fasta_str {
        my ($self) = @_;
        
        my $str .= ">" . $self->get_header() . "\n";
        $str .= $self->get_fwd_amplicon()->get_seq() . "\n";
        
        return $str;
    }
    
    sub get_rev_fasta_str {
        my ($self) = @_;
        
        my $str .= ">" . $self->get_header() . "\n";
        $str .= $self->get_rev_amplicon()->get_seq() . "\n";
        
        return $str;
    }
    
    sub get_fastq_str {
        my ($self) = @_;
        
        my $str .= $self->get_fwd_fastq_str();
        $str .= $self->get_rev_fastq_str();
        
        return $str;
    }
    
    sub get_fwd_fastq_str {
        my ($self) = @_;
        
        # Store this so I don't have to re-get it for the quality values
        my $fwd_amplicon = $self->get_fwd_amplicon();
        
        my $str .= "@" . $self->get_header() . "\n";
        $str .= $fwd_amplicon->get_seq() . "\n"; 
        $str .= "+\n";
        $str .= $fwd_amplicon->get_quals_str() . "\n";
        
        return $str;
    }
    
    sub get_rev_fastq_str {
        my ($self) = @_;
        
        # Store this so I don't have to re-get it for the quality values
        my $revAmplicon = $self->get_rev_amplicon();
        
        my $str .= "@" . $self->get_header(). "\n";
        $str .= $revAmplicon->get_seq() . "\n"; 
        $str .= "+\n";
        $str .= $revAmplicon->get_quals_str() . "\n";
        
        return $str;
    }
    
    sub get_quals_aref {
        my ($self) = @_;
        # This method doesn't really make sense.  I can't imagine any time I would
        # want to get the concatenated quality values for both the fwd and rev
        # sequence quality scores.  I put it in for continuity.
        
        my @quals = ( 
                      @{ $self->get_fwd_amplicon()->get_quals_aref() },
                      @{ $self->get_rev_amplicon()->get_quals_aref() }
                    );
        
        return \@quals;
    }
    
    sub get_fwd_quals_aref {
        my ($self) = @_;
        
        my @quals = ( 
                      @{ $self->get_fwd_amplicon()->get_quals_aref() },
                    );
    
        return \@quals;
    }
    
    sub get_rev_quals_aref {
        my ($self) = @_;
        
        my @quals = ( 
                      @{ $self->get_rev_amplicon()->get_quals_aref() }
                    );
        
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
=head1 MTToolbox::Match::PE::NoOverlap

MTToolbox::Match::PE::NoOverlap - A class to store paired end (PE) sequence that
                                  do NOT have overlapping ends and has defined
                                  attributes (eg. primer, amplicon, etc).

=head1 VERSION

This documentation refers to MTToolbox::Match::PE::NoOverlap version 4.1.2.

=head1 Included Modules

    Class::Std::Utils
    List::MoreUtils qw(any)
    Readonly
    Carp qw(carp croak)
    BioUtils::FastqSeq 1.0.0
    MyX::Generic 1.0.0
    version
    Match::PE

=head1 Inherit

Match::PE

=head1 SYNOPSIS
    
    use MTToolbox::Match::PE::NoOverlap;
    my $my_match = MTToolbox::Match::PE::NoOverlap->new({
                        id => $id,
                        desc => $desc,
                        fwd_tag => $fwd_tag,
                        fwd_linker => $fwd_linker,
                        fwd_primer => $fwd_primer,
                        fwd_amplicon => $fwd_amplicon,
                        fwd_tail => $fwd_tail,
                        rev_tag => $rev_tag,
                        rev_linker => $rev_linker,
                        rev_primer => $rev_primer,
                        rev_amplicon => $rev_amplicon,
                        rev_tail => $rev_tail,
                     });
    
    # Attribute Getters -- see also parent classes Match and Match::PE
    my $fwd_tail = $my_match->get_fwd_tail();
    my $revAmplicon = $my_match->get_rev_amplicon();
    my $revTail = $my_match->get_rev_tail();
    
    # Attribute Setters -- see also parent classes Match and Match::PE
    $my_match->set_fwd_tail($fwd_tail);
    $my_match->set_rev_amplicon($rev_amplicon);
    $my_match->set_rev_tail($rev_tail);
    
    # Other Getter Methods
    $my_match->get_fasta_str()
    $my_match->get_fwd_fasta_str()
    $my_match->get_rev_fasta_str()
    $my_match->get_fastq_str()
    $my_match->get_fwd_fastq_str()
    $my_match->get_rev_fastq_str()
    $my_match->get_quals_aref()
    $my_match->get_fwd_quals_aref()
    $my_match->get_rev_quals_aref()
    

=head1 DESCRIPTION

MTToolbox::Match::PE::NoOverlap is primarily a data structure to store the
attributes of paired end (PE) sequences that do NOT have overlapping ends and
that match a defined pattern.  That pattern is defined in the Sample object and
the various child classes of Sample.  All attributes that are passed into the
constructor are handled by the parent classes Match and Match::PE.

Because the PE sequence do not overlap they are stored as seperate FastqSeq
objects.  If they are stored seperately that also means they are processed
seperately.  You will notice that there are nearly identical methods for
both the forward and reverse reads.  If one of the methods is changed it is
important to change the corresponding forward or reverse method.

Because MTToolbox::Match::PE::NoOverlap inherits from Match::PE there are
several methods that must be implemented.  Some of these methods (those that
don't seperate the forward and reverse reads. eg. get_fasta_str, get_fastq_str,
and get_quals_aref) don't make sense in the context of non-overlapping
PE reads.  Check the individual methods for notes on this, and use these methods
with extreme CAUTION.

The various get_quals_aref methods are important.  They are important because
when I build a consensus sequence from a set of Match objects I need to access
their quality values.  Clustalw doesn't do anything with quality scores.  So
when I run Clustalw I only use the raw sequences.  Then when I make the
consensus I access the quality values by calling get_quals_aref.  I need the
quality values when making a consensus to resolve ties between two or more
aligned bases.

=head1 METHODS

=over

    new
    set_fwd_tail
    set_rev_amplicon
    set_rev_tail
    get_fwd_tail
    get_rev_amplicon
    get_rev_tail
    get_fasta_str
    get_fwd_fasta_str
    get_rev_fasta_str
    get_fastq_str
    get_fwd_fastq_str
    get_rev_fastq_str
    get_quals_aref
    get_fwd_quals_aref
    get_rev_quals_aref
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: Match->new({ id => $id,
                        desc => $desc,
                        fwd_tag => $fwd_tag,
                        fwd_linker => $fwd_linker,
                        fwd_primer => $fwd_primer,
                        fwd_amplicon => $fwd_amplicon,
                        fwd_tail => $fwd_tail,
                        rev_tag => $rev_tag,
                        rev_linker => $rev_linker,
                        rev_primer => $rev_primer,
                        rev_amplicon => $rev_amplicon,
                        rev_tail = $rev_tail,
                    });
    Function: Creates a new MTToolbox::Match::PE::NoOverlap object
    Returns: MTToolbox::Match::PE::NoOverlap
    Args: -id           => a string representing sequence id (NOT the molecule tag)
          -desc         => a string representing the sequence description
          -fwd_tag      => a FastqSeq representing the forward molecule tag
          -fwd_linker   => a FastqSeq representing the forward linker
          -fwd_primer   => a FastqSeq representing the forward primer
          -fwd_amplicon => a FastqSeq representing the forward amplicon
          -fwd_tail     => a FastqSeq representing the forward tail
          -rev_tag      => a FastqSeq representing the reverse molecule tag
          -rev_linker   => a FastqSeq representing the reverse linker
          -rev_primer   => a FastqSeq representing the reverse primer
          -rev_amplicon => a FastqSeq representing the reverse amplicion
          -rev_tail     => a FastqSeq representing the reverse tail
    Throws: MyX::Generic::Undef::Param
    Comments: NA
    See Also: BioUtils::FastqSeq, Match::PEmatch

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

=head2 set_rev_amplicon

    Title: set_rev_amplicon
    Usage: $my_match->set_rev_amplicon($rev_amplicon);
    Function: Sets the revAmplicon value
    Returns: 1 on successful completion
    Args: -rev_amplicon => a FastqSeq representing the reverse amplicon
    Throws: NA
    Comments: NA
    See Also: BioUtils::FastqSeq

=head2 set_rev_tail

    Title: set_rev_tail
    Usage: $my_match->set_rev_tail($rev_tail);
    Function: Sets the rev_tag value
    Returns: 1 on successful completion
    Args: -rev_tail => a FastqSeq representing the reverse tail
    Throws: NA
    Comments: NA
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

=head2 get_rev_amplicon
    
    Title: get_rev_amplicon
    Usage: my $id = $my_match->get_rev_amplicon();
    Function: Gets the reverse amplicon
    Returns: BioUtils::FastqSeq
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_rev_tail
    
    Title: get_rev_tail
    Usage: my $id = $my_match->get_rev_tail();
    Function: Gets the reverse tail
    Returns: BioUtils::FastqSeq
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_fasta_str
    
    Title: get_fasta_str
    Usage: my $fasta_Str = $my_match->get_fasta_str();
    Function: Get a string formatted for clustalw input
    Returns: String
    Args: None
    Throws: NA
    Comments: This method doesn't make much sense in the context of paired end
              non-overlapping sequneces.  It is strongly recommended that the
              corresponding forward and reverse methods be used to get sequences
              in FASTQ format.  One drawback to this method is that it doesn't
              distinguish the forward and reverse sequences in the id.  The
              first sequence in the return string is the forward sequences and
              the second is the reverse sequences.
    See Also: get_fwd_fasta_str, get_rev_fasta_str

=head2 get_fwd_fasta_str
    
    Title: get_fwd_fasta_str
    Usage: my $fasta_str = $my_match->get_fwd_fasta_str();
    Function: Get a string formatted for clustalw input of the forward read
    Returns: String
    Args: None
    Throws: NA
    Comments: This method and the corresonding reverse method
              (get_rev_fasta_str) should be used to get the clustalw
              sequence input instead of get_fasta_str.
    See Also: NA
        
=head2 get_rev_fasta_str
    
    Title: get_rev_fasta_str
    Usage: my $fasta_str = $my_match->get_rev_fasta_str();
    Function: Get a string formatted for clustalw input of the reverse read
    Returns: String
    Args: None
    Throws: NA
    Comments: This method and the corresonding forward method
              (get_fwd_fasta_str) should be used to get the clustalw
              sequence input instead of get_fasta_str.
    See Also: NA

=head2 get_fastq_str
    
    Title: get_fastq_str
    Usage: my $fastqStr = $my_match->get_fastq_str();
    Function: Get a fastq formatted string that can be printed to a file
    Returns: String
    Args: None
    Throws: NA
    Comments: This method doesn't make much sense in the context of paired end
              non-overlapping sequneces.  It is strongly recommended that the
              corresponding forward and reverse methods be used to get sequences
              in FASTQ format.  One drawback to this method is that it doesn't
              distinguish the forward and reverse sequences in the id.  The
              first sequence in the return string is the forward sequences and
              the second is the reverse sequences.
    See Also: get_fwd_fastq_str, get_rev_fastq_str

=head2 get_fwd_fastq_str
    
    Title: get_fwd_fastq_str
    Usage: my $fastqStr = $my_match->get_fwd_fastq_str();
    Function: Get a fastq formatted forward string that can be printed to a file
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_rev_fastq_str
    
    Title: get_rev_fastq_str
    Usage: my $fastqStr = $my_match->get_rev_fastq_str();
    Function: Get a fastq formatted reverse string that can be printed to a file
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_quals_aref
    
    Title: get_quals_aref
    Usage: my $quals_aref = $my_match->get_quals_aref();
    Function: Gets the quality values of the primer and amplicon attributes
    Returns: Array reference
    Args: None
    Throws: NA
    Comments: This method doesn't make much sense in the context of paired end
              non-overlapping sequneces.  It is strongly recommended that the
              corresponding forward and reverse methods be used to get quality
              score array refs.  This method currently returns an array with
              quality scores for the forward primer and amplicon followed by
              the reverse primer and amplicon.
    See Also: get_fwd_quals_aref, get_rev_quals_aref

=head2 get_fwd_quals_aref
    
    Title: get_fwd_quals_aref
    Usage: my $quals_aref = $my_match->get_fwd_quals_aref();
    Function: Gets the quality values of the forward primer and amplicon attributes
    Returns: Array reference
    Args: None
    Throws: NA
    Comments: This method returns only the fwd amplicon attribute quality values
              that this object stores.  To get additional or fewer parts of the
              sequence you can change the code of this method.  Athough if you
              change those parts you will also want to change regular, forward,
              and reverse methods for get_fastq_str and get_fasta_str.
    See Also: NA

=head2 get_rev_quals_aref
    
    Title: get_rev_quals_aref
    Usage: my $quals_aref = $my_match->get_rev_quals_aref();
    Function: Gets the quality values of the reverse primer and amplicon attributes
    Returns: Array reference
    Args: None
    Throws: NA
    Comments: This method returns only the rev amplicon attribute quality values
              that this object stores.  To get additional or fewer parts of the
              sequence you can change the code of this method.  Athough if you
              change those parts you will also want to change regular, forward,
              and reverse methods for get_fastq_str and get_fasta_str.
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
