package MTToolbox::MyX::MTParams;

use version; our $VERSION = qv('4.1.2');

use Exception::Class (
    'MTToolbox::MyX::MTParams' => {
    },
    
    'MTToolbox::MyX::MTParams::MissingRequiredTag' => {
        isa => 'MTToolbox::MyX::MTParams',
        fields => ['tag_name'],
    },
    
    'MTToolbox::MyX::MTParams::MissingTagValue' => {
        isa => 'MTToolbox::MyX::MTParams',
        fields => ['tag_name'],
    },
    
    'MTToolbox::MyX::MTParams::UnrecognizedTag' => {
        isa => 'MTToolbox::MyX::MTParams',
        fields => ['tag_name'],
    },
    
    'MTToolbox::MyX::MTParams::MissingExec' => {
        isa => 'MTToolbox::MyX::MTParams',
        fields => ['exec_name'],
    },
    
    'MTToolbox::MyX::MTParams::LongLinker' => {
        isa => 'MTToolbox::MyX::MTParams',
        fields => ['orientation', 'length'],
    },
    
    'MTToolbox::MyX::MTParams::UnrecognizedChar' => {
        isa => 'MTToolbox::MyX::MTParams',
        fields => ['char'],
    },
    
    'MTToolbox::MyX::MTParams::BadSeqType' => {
        isa => 'MTToolbox::MyX::MTParams',
        fields => ['value'],
    },
);

MTToolbox::MyX::MTParams->Trace(1);

1;
__END__


#######
# POD #
#######
=head1 NAME

MTToolbox::MyX::MTParams - A hierarchy of exceptions that can be used in
                           MTToolbox::MTParams

=head1 VERSION

This documentation refers to MTToolbox::MyX::MTParams version 4.1.2.

=head1 Included Modules

    Exception::Class

=head1 Inherit

    NA

=head1 SYNOPSIS

    # Throw a MTToolbox::MyX::MTParams exception
    use MTToolbox::MyX::MTParams;
    if ( ... ) {   # Some code looking for an error
        MTToolbox::MyX::MTParams->throw(
                            error => 'A MTParams execption'
                            );
    }
    
    # In caller catch the MTToolbox::MYX::MTParams exception
    eval { ... };
    if ( my $e = MTToolbox::MyX::MTParams->caught() ) {
        # Do something to handle the exception like print an error message
        print $e->error(), " via package ", $e->package(), " at ", $e->file,
            " line ", $e->line();
    }
    

=head1 DESCRIPTION

MTToolbox::MyX::MTParams holds a hierarchy of exception classes that can be used
in MTToolbox::MTParams objects

For more information what can be done when throwing and catching an exception
see Exception::Class and Exception::Class::Base.

=head1 CLASSES

=over

    MTToolbox::MyX::MTParams
    MTToolbox::MyX::MTParams::MissingRequiredTag
    MTToolbox::MyX::MTParams::MissingTagValue
    MTToolbox::MyX::MTParams::UnrecognizedTag
    MTToolbox::MyX::MTParams::MissingExec
    MTToolbox::MyX::MTParams::LongLinker
    MTToolbox::MyX::MTParams::UnrecognizedChar
    MTToolbox::MyX::MTParams::BadSeqType
        
    
=back

=head1 CLASSES DESCRIPTION

=head2 MTToolbox::MyX::MTParams
    
    Title: MTToolbox::MyX::MTParams
    Throw Usage: MTToolbox::MyX::MTParams->throw(
                    error => 'Any MTParams error message'
                );
    Catch Usage: if ( my $e = MTToolbox::MyX::MTParams->caught() ) { ... }
    Function: Throw/Catch a MTToolbox::MyX::MTParams exception
    Fields: error => an error message
    Inherits: NA
    Comments: NA
    See Also: NA

=head2 MTToolbox::MyX::MTParams::MissingRequiredTag

    Title: MTToolbox::MyX::MTParams::MissingRequiredTag
    Throw Usage: MTToolbox::MyX::MTParams::MissingRequiredTag->throw(
                    tag_name => $tag,
                );
    Catch Usage: if ( my $e = MTToolbox::MyX::MTParams::MissingRequiredTag->caught() ) { ... }
    Function: Throw/Catch a MTToolbox::MyX::MTParams::MissingRequiredTag
              exception when a required tag is missing.
    Fields: tag_name => The name of the missing tag
    Inherits: MTToolbox::MyX::MTParams
    Comments: NA
    See Also: NA
    
=head2 MTToolbox::MyX::MTParams::MissingTagValue

    Title: MTToolbox::MyX::MTParams::MissingTagValue
    Throw Usage: MTToolbox::MyX::MTParams::MissingTagValue->throw(
                    tag_name => $tag,
                );
    Catch Usage: if ( my $e = MTToolbox::MyX::MTParams::MissingTagValue->caught() ) { ... }
    Function: Throw/Catch a MTToolbox::MyX::MTParams::MissingTagValue exception
              when a tag is missing a value (i.e. it is undefined).
    Fields: tag_name => The name of the missing tag missing a value
    Inherits: MTToolbox::MyX::MTParams
    Comments: NA
    See Also: NA
    
=head2 MTToolbox::MyX::MTParams::UnrecognizedTag

    Title: MTToolbox::MyX::MTParams::UnrecognizedTag
    Throw Usage: MTToolbox::MyX::MTParams::UnrecognizedTag->throw(
                    tag_name => $tag,
                );
    Catch Usage: if ( my $e = MTToolbox::MyX::MTParams::UnrecognizedTag->caught() ) { ... }
    Function: Throw/Catch a MTToolbox::MyX::MTParams::UnrecognizedTag exception
              when a unrecognized tag is found.
    Fields: tag_name => The name of the unrecognized tag
    Inherits: MTToolbox::MyX::MTParams
    Comments: NA
    See Also: NA
    
=head2 MTToolbox::MyX::MTParams::MissingExec

    Title: MTToolbox::MyX::MTParams::MissingExec
    Throw Usage: MTToolbox::MyX::MTParams::MissingExec->throw(
                    exec_name => $exec,
                );
    Catch Usage: if ( my $e = MTToolbox::MyX::MTParams::MissingExec->caught() ) { ... }
    Function: Throw/Catch a MTToolbox::MyX::MTParams::MissingExec exception when
              an executable file cannot be found.
    Fields: exec_name => The name of the missing executable
    Inherits: MTToolbox::MyX::MTParams
    Comments: NA
    See Also: NA
    
=head2 MTToolbox::MyX::MTParams::LongLinker

    Title: MTToolbox::MyX::MTParams::LongLinker
    Throw Usage: MTToolbox::MyX::MTParams::LongLinker->throw(
                    orientation => 'fwd',
                    length => 25,
                );
    Catch Usage: if ( my $e = MTToolbox::MyX::MTParams::LongLinker->caught() ) { ... }
    Function: Throw/Catch a MTToolbox::MyX::MTParams::LongLinker exception when
              an linker is very long.  Perhaps the user made a mistake because
              linkers are generally very short (i.e. 2-5 bp).
    Fields: orientation => The fwd or rev linker
            length => The length of the linker
    Inherits: MTToolbox::MyX::MTParams
    Comments: NA
    See Also: NA
    
=head2 MTToolbox::MyX::MTParams::UnrecognizedChar

    Title: MTToolbox::MyX::MTParams::UnrecognizedChar
    Throw Usage: MTToolbox::MyX::MTParams::UnrecognizedChar->throw(
                    char => $char,
                );
    Catch Usage: if ( my $e = MTToolbox::MyX::MTParams::UnrecognizedChar->caught() ) { ... }
    Function: Throw/Catch a MTToolbox::MyX::MTParams::UnrecognizedChar exception
              when an character is unrecognized or not allowed in that context.
    Fields: char => The character that is unrecognized or not allowed
    Inherits: MTToolbox::MyX::MTParams
    Comments: NA
    See Also: NA
    
=head2 MTToolbox::MyX::MTParams::BadSeqType

    Title: MTToolbox::MyX::MTParams::BadSeqType
    Throw Usage: MTToolbox::MyX::MTParams::BadSeqType->throw(
                    value => $seq_type,
                );
    Catch Usage: if ( my $e = MTToolbox::MyX::MTParams::BadSeqType->caught() ) { ... }
    Function: Throw/Catch a MTToolbox::MyX::MTParams::BadSeqType exception when
              an sequence type is unrecognized.  Seq type can be set as either
              1, 2, or 3 for SE, PE w/ overlap, or PE w/o overlap respectively
    Fields: value => The unaccepted seq type value
    Inherits: MTToolbox::MyX::MTParams
    Comments: NA
    See Also: NA
        
=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to the author


=head1 AUTHOR

Scott Yourstone  C<< <scott.yourstone81@gmail.com> >>


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

