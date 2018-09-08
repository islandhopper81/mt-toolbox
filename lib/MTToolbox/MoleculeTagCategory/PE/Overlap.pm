package MTToolbox::MoleculeTagCategory::PE::Overlap;

use strict;
use warnings;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(carp croak);
use MyX::Generic 1.0.0;
use version; our $VERSION = qv('4.1.2');
use BioUtils::ConsensusBuilder::ConsensusBuilder 1.0.7 qw(build_consensus);

use MTToolbox::MoleculeTagCategory;
use base qw(MTToolbox::MoleculeTagCategory);

{
    Readonly my $NEW_USAGE => q{ new( {tag => ,
                                       } ) };

    # Attributes #
    
    # Setters #
    
    # Getters #
    
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
        if ( any {!defined $_} $arg_href->{tag},
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
   
    
    
    ##########
    # Others #
    ##########
}

1;
__END__

#######
# POD #
#######
=head1 MTToolbox::MoleculeTagCategory::PE::Overlap

MTToolbox::MoleculeTagCategory::PE::Overlap - A class to store and operate on a
set of paired-end (PE) reads that have the same molecule tag

=head1 VERSION

This documentation refers to MTToolbox::MoleculeTagCategory::PE::Overlap version 4.1.2.

=head1 Included Modules

Class::Std::Utils
List::MoreUtils qw(any)
Readonly
Carp qw(carp croak)
MyX::Generic 1.0.0
version
MTToolbox::MoleculeTagCategory

=head1 Inherit

MTToolbox::MoleculeTagCategory

=head1 SYNOPSIS
    
    use MTToolbox::MoleculeTagCategory::PE::Overlap;
    my $my_mt_bin = MTToolbox::MoleculeTagCategory::PE::Overlap->new(
                        "ATTGATGAG_ATGG"
                    );
    
    # See documentation for parent class

=head1 DESCRIPTION

MTToolbox::MoleculeTagCategory::PE::Overlap is a child class of
MoleculeTagCategory.  It represents a set of paired-end (PE) reads that have the
same molecule tag.  

=head1 METHODS

=over

    new
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: MTToolbox::MoleculeTagCategory::PE::Overlap->new($tag);
    Function: Creates a new MTToolbox::MoleculeTagCategory::PE::Overlap object
    Returns: MTToolbox::MoleculeTagCategory::PE::Overlap
    Args: -tag => a unique identifier representing the molecule tag(s) of a read
    Throws: NA
    Comments: NA
    See Also: MTToolbox::MoleculeTagCategory

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