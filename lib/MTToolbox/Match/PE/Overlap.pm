package MTToolbox::Match::PE::Overlap;


use strict;
use warnings;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(carp croak);
use BioUtils::FastqSeq 1.0.0;
use MyX::Generic 1.0.0;
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
                                       rev_tag => ,
                                       rev_linker => ,
                                       rev_primer => ,
                                       rev_amplicon => ,
                                       } ) }
                                       ;

    # Attributes #
    
    # Setters #
    
    # Getters #
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
        if ( any {!defined $_}
                $arg_href->{id},
                $arg_href->{desc},
                $arg_href->{fwd_tag},
                $arg_href->{fwd_linker}, 
                $arg_href->{fwd_primer},
                $arg_href->{fwd_amplicon},
                $arg_href->{rev_tag},
                $arg_href->{rev_linker},
                $arg_href->{rev_primer},
            ) {
            MyX::Generic::Undef::Param->throw(
                                              error => 'Undefined parameter value',
                                              usage => $NEW_USAGE,
                                              );
        }
        
        # Bless a scalar to instantiate an object
        my $new_obj = $class->SUPER::new($arg_href);
        my $ident = ident($new_obj);
        
        return $new_obj;
    }
    
    
    ###########
    # Setters #
    ###########
    
    
    ###########
    # Getters #
    ###########
    sub get_fasta_str {
        my ($self) = @_;
        
        my $str .= ">" . $self->get_header() . "\n";
        $str .= $self->get_fwd_amplicon()->get_seq() .
                "\n";
        
        return $str;
    }
    
    sub get_fastq_str {
        my ($self) = @_;
        
        # The sequence id header
        my $str .= "@" . $self->get_header() . "\n";
        
        # The sequence
        $str .= $self->get_fwd_amplicon()->get_seq() . 
                "\n";
                
        # The quality value header  
        $str .= "+\n";
        
        # The quality values
        $str .= $self->get_fwd_amplicon()->get_quals_str() .
                "\n";
        
        return $str;
    }
    
    sub get_quals_aref {
        my ($self) = @_;
        
        my @quals = (
                     @{ $self->get_fwd_amplicon() -> get_quals_aref()   },
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
=head1 MTToolbox::Match::PE::Overlap

MTToolbox::Match::PE::Overlap - A class to store paired end (PE) sequence that
                                do NOT have overlapping ends and has defined
                                attributes (eg. primer, amplicon, etc).

=head1 VERSION

This documentation refers to MTToolbox::Match::PE::Overlap version 4.1.2.

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
    
    use MTToolbox::Match::PE::Overlap;
    my $my_match = MTToolbox::Match::PE::Overlap->new({
                        id => $id,
                        desc => $desc,
                        fwd_tag => $fwd_tag,
                        fwd_linker => $fwd_linker,
                        fwd_primer => $fwd_primer,
                        fwd_amplicon => $fwd_amplicon,
                        rev_tag => $rev_tag,
                        rev_linker => $rev_linker,
                        rev_primer => $rev_primer,
                        rev_amplicon => $rev_amplicon,
                    });
    
    # Attribute Getters are defined in parent classes Match and Match::PE
    
    # Attribute Setters are defined in parent classes Match and Match::PE
    
    # Other Getter Methods
    $my_match->get_fasta_str()
    $my_match->get_fastq_str()
    $my_match->get_quals_aref()

    

=head1 DESCRIPTION

MTToolbox::Match::PE::Overlap is primarily a data structure to store the attributes
of paired end (PE) sequences that DO have overlapping ends and that match a
defined pattern.  That pattern is defined in the Sample object and the various
child classes of Sample.  All attributes that are passed into the constructor
are handled by the parent classes Match and Match::PE.  By the time these
objects are created the PE overlapping reads have been merged.  So a PE
overlapping read is a single sequence that has attributes that are akin to PE
reads.

The get_quals_aref method is important because when the consensus sequence is
built from a set of Match objects I need to access their quality values.
Clustalw doesn't do anything with quality scores.  So when I run Clustalw I only
use the raw sequences.  Then when I make the consensus I access the quality
values by calling get_quals_aref.  I need the quality values when making a
consensus to resolve ties between two or more aligned bases.

=head1 METHODS

=over

    new
    get_fasta_str
    get_fastq_str
    get_quals_aref
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: MTToolbox::Match::PE::Overlap->new({
                id => $id,
                desc => $desc,
                fwd_tag => $fwd_tag,
                fwd_linker => $fwd_linker,
                fwd_primer => $fwd_primer,
                fwd_amplicon => $fwd_amplicon,
                rev_tag => $rev_tag,
                rev_linker => $rev_linker,
                rev_primer => $rev_primer,
                rev_amplicon => $rev_amplicon,
            });
    Function: Creates a new MTToolbox::Match::PE::Overlap object
    Returns: MTToolbox::Match::PE::Overlap
    Args: -id           => a string representing the sequence id (NOT the molecule tag)
          -desc         => a string representing the sequence description
          -fwd_tag      => a FastqSeq representing the forward molecule tag
          -fwd_linker   => a FastqSeq representing the forward linker
          -fwd_primer   => a FastqSeq representing the forward primer
          -fwd_amplicon => a FastqSeq representing the forward amplicon
          -rev_tag      => a FastqSeq representing the reverse molecule tag
          -rev_linker   => a FastqSeq representing the reverse linker
          -rev_primer   => a FastqSeq representing the reverse primer
          -rev_amplicon => a FastqSeq representing the reverse amplicion
    Throws: MyX::Generic::Undef::Param
    Comments: NA
    See Also: BioUtils::FastqSeq, Match::PE
    
=head2 get_fasta_str
    
    Title: get_fasta_str
    Usage: my $fasta_Str = $my_match->get_fasta_str();
    Function: Get a string formatted for clustalw input
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: get_fwd_fasta_str, get_rev_fasta_str

=head2 get_fastq_str
    
    Title: get_fastq_str
    Usage: my $fastqStr = $my_match->get_fastq_str();
    Function: Get a fastq formatted string that can be printed to a file
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: get_fwd_fastq_str, get_rev_fastq_str

=head2 get_quals_aref
    
    Title: get_quals_aref
    Usage: my $quals_aref = $my_match->get_quals_aref();
    Function: Gets the quality values of the amplicon attribute
    Returns: Array reference
    Args: None
    Throws: NA
    Comments: NA
    See Also: get_fwd_quals_aref, get_rev_quals_aref
    

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
