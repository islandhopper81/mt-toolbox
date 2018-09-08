package MTToolbox::MTEnv;

use warnings;
use strict;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(carp croak);
use IPC::Cmd qw( can_run );
use Exception::Class;
use MyX::Generic 0.0.3;
use MTToolbox::MyX::MTParams 4.1.2;
use version; our $VERSION = qv('4.1.2');

{
	# Usage statement
	Readonly my $NEW_USAGE => q{ new() };

    # Subroutines #
	sub check_env;
	
    sub check_cat;
    sub check_rm;
    sub check_flash;
    sub check_gnuplot;
    sub check_clustalw;
    
    
    ###############
    # Constructor #
    ###############
    sub new {
        my ($class) = @_;
        
        # Croak if calling 'new' on an already blessed reference
        croak 'Constructor called on existing object instead of class'
            if ref $class;
        
        # Bless a scalar to instantiate an object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        return $new_obj;
    }
    
    ###############
	# Subroutines #
	###############
    sub check_env {
		my ($self) = @_;
		
		$self->check_cat();
		$self->check_rm();
		$self->check_flash();
		$self->check_gnuplot();
		$self->check_clustalw();
		
		return 1;
	}
	
    sub check_cat {
		my ($self) = @_;
		
		if ( ! can_run('cat') ) {
			MTToolbox::MyX::MTParams::MissingExec->throw(
                    error => "Cannot find cat system command",
                    exec_name => 'cat',
                );
		}
		
		return 1;
	}
	
    sub check_rm {
		my ($self) = @_;
		
		if ( ! can_run('rm') ) {
			MTToolbox::MyX::MTParams::MissingExec->throw(
                    error => "Cannot find rm system command",
                    exec_name => 'rm',
                );
		}
		
		return 1;
	}
	
    sub check_flash {
		my ($self) = @_;
		
		if ( ! can_run('flash') ) {
			MTToolbox::MyX::MTParams::MissingExec->throw(
                    error => "Cannot find flash command",
                    exec_name => 'flash',
                );
		}
		
		return 1;
	}
	
    sub check_gnuplot {
		my ($self) = @_;
		
		if ( ! can_run('gnuplot') ) {
			MTToolbox::MyX::MTParams::MissingExec->throw(
                    error => "Cannot find gnuplot command",
                    exec_name => 'gnuplot',
                );
		}
		
		return 1;
	}
	
    sub check_clustalw {
		my ($self) = @_;
		
		if ( ! can_run('clustalw') ) {
			MTToolbox::MyX::MTParams::MissingExec->throw(
                    error => "Cannot find clustalw command",
                    exec_name => 'clustalw',
                );
		}
		
		return 1;
	}
}



1; # Magic true value required at end of module
__END__

=head1 NAME

MTToolbox::MTEnv - A module for checking the working environment for
				   required executables.


=head1 VERSION

This document describes MTToolbox::MTEnv version 4.1.2


=head1 SYNOPSIS

    use MTEnv;
	
    # Build the object
	my $mt_env = MTToolbox::MTEnv->new();
    
    # Run all checks
	$mt_env->check_env();
	

  
=head1 DESRCIPTION

MTEnv checks for all required excutables.


=head1 CONFIGURATION AND ENVIRONMENT

    None


=head1 DEPENDENCIES
	
	Class::Std::Utils
	List::MoreUtils qw(any)
	Readonly
	Carp qw(carp croak)
	IPC::Cmd qw( can_run )
	Exception::Class
	MyX::Generic 1.0.7
	MTToolbox::MyX::MTParams 4.1.2
	version our $VERSION = qv('4.1.2')


=head1 INCOMPATIBILITIES

    None reported.

=head1 METHODS

	check_env
	
    check_cat
    check_rm
    check_flash
    check_gnuplot
    check_clustalw
    
=head1 METHODS DESRCIPTION

=head2 new

	Title: new
	Usage: my $mt_env = MTToolbox::MTEnv->new();
	Function: Create a new MTToolbox::MTEnv object
	Returns: Reference to a MTToolbox::MTEnv object
	Args: NA
	Throws: croak 'Constructor called on existing object instead of class'
	Comments: NA
	See Also: NA

=head2 check_env

	Title: check_env
	Usage: $mt_env->check_env()
	Function: Runs all check subroutines
	Returns: 1 on successful completion
	Args: NA
	Throws: NA
	Comments: If any of the the called subroutines has an error with throw
			  MTToolbox::MyX::MTParams::MissingExec.
	See Also: NA
	
=head2 check_cat

	Title: check_cat
	Usage: $mt_env->check_cat()
	Function: Checks for cat system command
	Returns: 1 on successful completion
	Args: NA
	Throws: MTToolbox::MyX::MTParams::MissingExec.
	Comments: NA
	See Also: NA
	
=head2 check_rm

	Title: check_rm
	Usage: $mt_env->check_rm()
	Function: Checks for rm system command
	Returns: 1 on successful completion
	Args: NA
	Throws: MTToolbox::MyX::MTParams::MissingExec.
	Comments: NA
	See Also: NA
	
=head2 check_flash

	Title: check_flash
	Usage: $mt_env->check_flash()
	Function: Checks for flash command
	Returns: 1 on successful completion
	Args: NA
	Throws: MTToolbox::MyX::MTParams::MissingExec.
	Comments: NA
	See Also: NA
	
=head2 check_gnuplot

	Title: check_gnuplot
	Usage: $mt_env->check_gnuplot()
	Function: Checks for gnuplot command
	Returns: 1 on successful completion
	Args: NA
	Throws: MTToolbox::MyX::MTParams::MissingExec.
	Comments: NA
	See Also: NA
	
=head2 check_clustalw

	Title: check_clustalw
	Usage: $mt_env->check_clustalw()
	Function: Checks for clustalw command
	Returns: 1 on successful completion
	Args: NA
	Throws: MTToolbox::MyX::MTParams::MissingExec.
	Comments: NA
	See Also: NA


=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to author


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

