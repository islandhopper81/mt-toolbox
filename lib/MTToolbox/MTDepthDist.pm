package MTToolbox::MTDepthDist;

use strict;
use warnings;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(carp croak);
use Chart::Gnuplot;
use MyX::Generic 0.0.3;
use version; our $VERSION = qv('4.1.2');


{
    Readonly my $NEW_USAGE => q{ new( {dist_file => } ) };

    # Attributes #
    my %dist_file_of;
    my %hist_href_of;
    
    # Setters #
    sub set_dist_file;
    sub set_hist_href;
    
    # Getters #
    sub get_dist_file;
    sub set_hist_href;
    
    # Others #
    sub build_hist;
    sub print_hist;
    sub print_hist_graph;
    sub hist_to_string;
    
    
    
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
        
        # Initialize the objects attributes
        if ( defined $arg_href->{dist_file} ) {
            $dist_file_of{ident $new_obj} = $arg_href->{dist_file};
        }
        
        return $new_obj;
    }
    
    
    ###########
    # Setters #
    ###########
    sub set_dist_file {
        my ($self, $file) = @_;
        $dist_file_of{ident $self} = $file;
        return 1;
    }
    
    sub set_hist_href {
        my ($self, $href) = @_;
        $hist_href_of{ident $self} = $href;
        return 1;
    }

    
    ###########
    # Getters #
    ###########
    sub get_dist_file {
        my ($self) = @_;
        return $dist_file_of{ident $self};
    }
    
    sub get_hist_href {
        my ($self) = @_;
        return $hist_href_of{ident $self};
    }
    
    
    
    ##########
    # Others #
    ##########
    sub build_dist {
        my ($self, $mt_objs_href, $sample_id) = @_;
        
        my $file = $self->get_dist_file();
        open my $OUT, ">", $file or
            croak("Cannot open file: $file");
        
        foreach my $mt ( keys %{$mt_objs_href} ) {
            print $OUT $sample_id . "_" .
                       $mt . "\t" .
                       $mt_objs_href->{$mt}->get_seq_count() .
                       "\n";
        }
        
        close($OUT);
        
        return 1;
    }
    
    sub build_hist {
        my ($self) = @_;
        
        my %hist = ();
        my $file = $self->get_dist_file();
        
        open my $IN, "<", $file or
            croak("Cannot open hist file: $file");
        
        foreach my $line ( <$IN> ) {
            chomp $line;
            
            my @values = split /\s+/, $line;
            
            if ( defined $hist{$values[1]} ) {
                $hist{$values[1]}++;
            }
            else {
                $hist{$values[1]} = 1;
            }
        }
        
        # set the hist href
        $hist_href_of{ident $self} = \%hist;
        
        return \%hist;
    }

    sub print_hist {
        my ($self, $file) = @_;
        
        my $hist = $self->get_hist_href();
        
        open my $OUT, ">", $file or
            croak("Cannot open file: $file");
        
        print $OUT $self->hist_to_string();
        
        close($OUT);
        
        return 1;
    }
    
    sub print_hist_graph {
        my ($self, $file, $title) = @_;
        
        # variables for graph step
        my @x_vals = (); 
        my @y_vals = ();
        my $chart = Chart::Gnuplot->new(
            output => $file,
            title => $title,
            xlabel => "MT Depth",
            ylabel => "Count",
        );
        $chart->command("set boxwidth .9");
        $chart->command("set style fill solid");
        
        
        # get the x and y values
        my $href = $self->get_hist_href();
        foreach my $key ( sort {$a<=>$b} keys %{$href} ) {
            push @x_vals, $key;
            push @y_vals, $href->{$key};
        }
        
        # print hist graph
        my $data_set = Chart::Gnuplot::DataSet->new(
            xdata => \@x_vals,
            ydata => \@y_vals,
            style => "boxes",
        );
        $chart->plot2d($data_set);
        
        return 1;
    }
    
    sub hist_to_string {
        my ($self) = @_;
        
        my $hist = $self->get_hist_href();
        my $str;
        
        foreach my $key ( sort {$a <=> $b} keys %{$hist} ) {
            $str .= $key . "\t" . $hist->{$key} . "\n";
        }
        
        return $str;
    }
}

1;
__END__

#######
# POD #
#######
=head1 NAME

MTToolbox::MTDepthDist - An object for handling molecule tag depth distributions

=head1 VERSION

This documentation refers to MTToolbox::MTDepthDist version 4.1.2.

=head1 Included Modules

    Class::Std::Utils
    List::MoreUtils qw(any)
    Readonly
    Carp qw(carp croak)
    Chart::Gnuplot
    version our $VERSION = qv('4.1.2')
    MyX::Generic 1.0.0

=head1 Inherit

    NA

=head1 SYNOPSIS

    use MTToolbox::MTDepthDist;
    my $mt_depth_dist = MTToolbox::MTDepthDist->new({dist_file => $file});
    
    $mt_depth_dist->build_hist();
    $mt_depth_dist->print_hist("seqs_per_mt_hist.txt");
    $mt_depth_dist->print_hist_graph("seqs_per_mt_hist.ps", "MT Depth");
    

=head1 DESCRIPTION

MTToolbox::MTDepthDist is an object for handling molecule tag depth
distributions. Specifically, it can create the distributions from one
seqs_per_mt_dist.txt files that are output for each sample and all samples
combined when MTToolbox is run.  It also used gnuplot to create graphs to show
the molecule tag depth distribution histogram.

=head1 METHODS

=over

    new
    set_dist_file
    set_hist_href
    get_dist_file
    get_hist_href
    build_dist
    build_hist
    print_hist
    print_hist_graph
    hist_to_string
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: MTToolbox::MTDepthDist->new({dist_file => $file});
    Function: Creates a new Sample object
    Returns: MTToolbox::MTDepthDist
    Args: -dist_file => file path to the seqs_per_mt_dist.txt file

=head2 set_dist_file

    Title: set_dist_file
    Usage: $depth_dist->set_dist_file($file);
    Function: Sets the file path to the seqs_per_mt_dist.txt file
    Returns: 1 on successful completion
    Args: -file => file path to the seqs_per_mt_dist.txt file
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 set_hist_href

    Title: set_hist_href
    Usage: $depth_dist->set_hist_href($href);
    Function: Sets the histogram hash ref
    Returns: 1 on successful completion
    Args: -href => a hash ref to a histogram
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_dist_file

    Title: get_dist_file
    Usage: my $id = $depth_dist->get_dist_file();
    Function: Gets the file path to the seqs_per_mt_dist.txt file
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_hist_href

    Title: get_hist_href
    Usage: my $id = $depth_dist->get_hist_href();
    Function: Gets the hash ref to the histogram
    Returns: Href
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 build_dist

    Title: build_dist
    Usage: $depth_dist->build_dist($mt_objs_href, $sample_id);
    Function: Prints the dist file given the mt objects and sample id
    Returns: 1 on successful completion
    Args: -mt_obj_href => a hash ref with MoleculeTagCategory objects
          -sample_id => the sample id to be prepended to the sample name
    Throws: NA
    Comments: the dist is not saved in this object as the hist is
    See Also: NA
    
=head2 build_hist

    Title: build_hist
    Usage: $depth_dist->build_hist();
    Function: Builds and stores the histogram hash ref
    Returns: Href
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 print_hist

    Title: print_hist
    Usage: $depth_dist->print_hist($file);
    Function: Prints the histogram stored in the object
    Returns: 1 on successful completion
    Args: -file => file path to print the histogram
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 print_hist_graph

    Title: print_hist_graph
    Usage: $depth_dist->print_hist_graph($file, $title);
    Function: Print the histogram
    Returns: 1 on successful completion
    Args: -file  => path to the file to print to
          -title => a title for the graph
    Throws: NA
    Comments: Uses gnuplot
    See Also: NA
    
=head2 hist_to_string

    Title: hist_to_string
    Usage: $depth_dist->hist_to_string();
    Function: Get the histogram in string format
    Returns: String
    Args: NA
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
