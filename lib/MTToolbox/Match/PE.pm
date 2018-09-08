package MTToolbox::Match::PE;

use strict;
use warnings;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(croak);
use BioUtils::FastqSeq 1.0.0;
use MyX::Generic 0.0.3;
use version; our $VERSION = qv('4.1.2');

use MTToolbox::Match;
use base qw(MTToolbox::Match);


{
    Readonly my $NEW_USAGE => q{ new( {id => ,
                                       desc => ,
                                       fwd_tag => ,
                                       fwd_linker => ,
                                       fwd_primer => ,
                                       rev_tag => ,
                                       rev_linker => ,
                                       rev_primer => ,
                                       } ) }
                                       ;
    
    # Attributes #
    my %rev_tag_of;
    my %rev_linker_of;
    my %rev_primer_of;

    # Setters #
    sub set_rev_tag;
    sub set_rev_linker;
    sub set_rev_primer;
    
    # Getters #
    sub get_rev_tag;
    sub get_rev_linker;
    sub get_rev_primer;
    
    
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
        
        # Initialize attributes specific to this class
        $rev_tag_of{$ident} = $arg_href->{rev_tag};
        $rev_linker_of{$ident} = $arg_href->{rev_linker};
        $rev_primer_of{$ident} = $arg_href->{rev_primer};
        
        return $new_obj;
    }
    
    
    
    ###########
    # Setters #
    ###########
    sub set_rev_tag {
        my ($self, $rev_tag) = @_;
        $rev_tag_of{ident $self} = $rev_tag;
        return 1;
    }
    
    sub set_rev_linker {
        my ($self, $rev_linker) = @_;
        $rev_linker_of{ident $self} = $rev_linker;
        return 1;
    }
    
    sub set_rev_primer {
        my ($self, $rev_primer) = @_;
        $rev_primer_of{ident $self} = $rev_primer;
        return 1;
    }
    
    
    ###########
    # Getters #
    ###########
    sub get_rev_tag {
        my ($self) = @_;
        return $rev_tag_of{ident $self};
    }
    
    sub get_rev_linker {
        my ($self) = @_;
        return $rev_linker_of{ident $self};
    }
    
    sub get_rev_primer() {
        my ($self) = @_;
        return $rev_primer_of{ident $self};
    }
}


1;
__END__


#######
# POD #
#######
=head1 MTToolbox::Match::PE

MTToolbox::Match::PE - An abstract class to store paired end (PE) sequence that
                       have defined attributes (eg. primer, amplicon, etc).

=head1 VERSION

This documentation refers to MTToolbox::Match::PE version 4.1.2.

=head1 Included Modules

    Class::Std::Utils
    List::MoreUtils qw(any)
    Readonly
    Carp qw(carp croak)
    BioUtils::FastqSeq 1.0.0
    MyX::Generic 1.0.0
    version
    MTToolbox::Match

=head1 Inherit

    MTToolbox::Match

=head1 SYNOPSIS
    
    use MTToolbox::Match::PE;
    my $my_match = MTToolbox::Match::PE->new({
                        id => $id,
                        desc => $desc,
                        fwd_tag => $fwd_tag,
                        fwd_linker => $fwd_linker,
                        fwd_primer => $fwd_primer,
                        fwd_amplicon => $fwd_amplicon,
                        rev_tag => $rev_tag,
                        rev_linker => $rev_linker,
                        rev_primer => $rev_primer,
                    });
    
    # Attribute Getters (see also parent class MTToolbox::Match)
    my $rev_tag = $my_match->get_rev_tag();
    my $rev_linker = $my_match->get_rev_linker();
    my $revPrimer = $my_match->get_rev_primer();

    
    # Attribute Setters (see also parent class MTToolbox::Match)
    $my_match->set_rev_tag($rev_tag);
    $my_match->set_rev_spcer($rev_linker);
    $my_match->set_rev_primer($rev_primer);
    

=head1 DESCRIPTION

MTToolbox::Match::PE is primarily a data structure to store the attributes of
paired end (PE) sequences that match a defined pattern.  That pattern is defined
in the Sample object and the various child classes of Sample.
MTToolbox::Match::PE is like an abstract class.  The attributes that are stored
in MTToolbox::Match::PE, which will be in all MTToolbox::Match::PE objects,
include rev_tag, rev_primer and all other attributes stored in the parent class
object, MTToolbox::Match.

NOTE: I remvoed the rev_amplicon and rev_tail attributes from this class because
MTToolbox::Match::PE::Overlap objects do not have these attributes.  So the
rev_amplicon and rev_tail attributes are defined in the
MTToolbox::Match::PE::NoOverlap class.

=head1 METHODS

=over

    new
    set_rev_tag
    set_rev_linker
    set_rev_primer
    get_rev_tag
    get_rev_linker
    get_rev_primer
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: MTToolbox::Match::PE->new({
                id => $id,
                desc => $desc,
                fwd_tag => $fwd_tag,
                fwd_linker => $fwd_linker,
                fwd_primer => $fwd_primer,
                fwd_amplicon => $fwd_amplicon,
                rev_tag => $rev_tag,
                rev_linker => $rev_linker,
                rev_primer => $rev_primer,
            });
    Function: Creates a new MTToolbox::Match::PE object
    Returns: MTToolbox::Match::PE
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
    See Also: BioUtils::FastqSeq

=head2 set_rev_tag

    Title: set_rev_tag
    Usage: $my_match->set_rev_tag($rev_tag);
    Function: Sets the rev_tag value
    Returns: 1 on successful completion
    Args: -rev_tag => a FastqSeq representing the reverse tag
    Throws: NA
    Comments: NA
    See Also: BioUtils::FastqSeq
    
=head2 set_rev_linker

    Title: set_rev_linker
    Usage: $my_match->set_rev_linker($rev_linker);
    Function: Sets the rev_spcare value
    Returns: 1 on successful completion
    Args: -rev_linker => a FastqSeq representing the reverse linker
    Throws: NA
    Comments: NA
    See Also: BioUtils::FastqSeq

=head2 set_rev_primer

    Title: set_rev_primer
    Usage: $my_match->set_rev_primer($rev_primer);
    Function: Sets the revPrimer value
    Returns: 1 on successful completion
    Args: -rev_primer => a FastqSeq representing the reverse primer
    Throws: NA
    Comments: NA
    See Also: BioUtils::FastqSeq
  
=head2 get_rev_tag
    
    Title: get_rev_tag
    Usage: my $id = $my_match->get_rev_tag();
    Function: Gets the reverse tag
    Returns: BioUtils::FastqSeq
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_rev_linker
    
    Title: get_rev_linker
    Usage: my $id = $my_match->get_rev_linker();
    Function: Gets the reverse linker
    Returns: BioUtils::FastqSeq
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_rev_primer
    
    Title: get_rev_primer
    Usage: my $id = $my_match->get_rev_primer();
    Function: Gets the reverse primer
    Returns: BioUtils::FastqSeq
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



