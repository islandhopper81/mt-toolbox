package MTToolbox::SampleSummary::SE;

# This object is a data structure that holds summary information for
# a sample with single-end reads. 


use strict;
use warnings;

use Class::Std::Utils;
use Carp qw(carp croak);
use version; our $VERSION = qv('4.1.2');
use MTToolbox::SampleSummary;
use base qw(MTToolbox::SampleSummary);  # inherits from MTToolbox::SampleSummary


{
    
    # Attributes #

    # Setters #
    
    # Getters #
    
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
    
    
    ###########
    # Getters #
    ###########
    
    
    ##########
    # Others #
    ##########
    sub to_string() {
        my ($self) = @_;
        
        my $str = "Sample type: SE\n";
        $str .= $self->SUPER::to_string();
        
        return $str;                                                        
    }
}


1;
__END__

#######
# POD #
#######
=head1 MTToolbox::SampleSummary::SE

MTToolbox::SampleSummary::SE - A class to store and print summary inforamtion
about a Sample with single-end reads.

=head1 VERSION

This documentation refers to MTToolbox::SampleSummary::SE version 4.1.2.

=head1 Included Modules

Class::Std::Utils;
Carp qw(carp croak)
version
MTToolbox::SampleSummary

=head1 Inherit

MTToolbox::SampleSummary

=head1 SYNOPSIS
    
    use MTToolbox::SampleSummary::SE;
    my $my_sum_info = MTToolbox::SampleSummary::SE->new();
    
    # Attribute Setters -- See also MTToolbox::SampleSummary
    
    # Attribute Getters -- See also MTToolbox::SampleSummary
    
    # Printing Methods -- See also MTToolbox::SampleSummary
    $my_sum_info->to_string();
    

=head1 DESCRIPTION

This class stores and outputs summary information about a Sample object
with SE reads.

=head1 METHODS

=over

    new
    to_string
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: MTToolbox::SampleSummary::SE->new();
    Function: Creates a new MTToolbox::SampleSummary::SE object
    Returns: MTToolbox::SampleSummary::SE
    Args: NA
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
    Comments: Overwritten from parent class
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
