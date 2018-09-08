package MTToolbox::MoleculeTagCategory::PE::NoOverlap;

use strict;
use warnings;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(carp croak);
use MyX::Generic 0.0.3;
use BioUtils::ConsensusBuilder::ConsensusBuilder 1.0.0 qw(build_con_from_file build_from_clustalw_file build_consensus);
use version; our $VERSION = qv('4.1.2');

use MTToolbox::MoleculeTagCategory;
use base qw(MTToolbox::MoleculeTagCategory);

{
    Readonly my $NEW_USAGE => q{ new( {tag => ,
                                       } ) };

    # Attributes #
    my %fwd_simple_align_of;
    my %rev_simple_align_of;
    my %fwd_fastq_con_of;
    my %rev_fastq_con_of;
    
    # Setters #
    sub set_fwd_fastq_con;
    sub set_rev_fastq_con;
    
    # Getters #
    sub get_unaligned_seqs_fasta_str;
    sub get_unaligned_seqs_fastq_str;
    sub get_unaligned_quals_href;
    sub get_fwd_fastq_con;
    sub get_rev_fastq_con;
    sub get_con_fastq_str;
    sub _get_fwd_con_fastq_str;
    sub _get_rev_con_fastq_str;
    
    # Others #
    sub build_con;
    sub run_command_line_muscle;
    sub run_command_line_clustalw;
    sub build_con_from_fasta_file;
    sub _build_fwd_con_from_fasta_file;
    sub _build_rev_con_from_fasta_file;
    sub build_con_from_clustalw_file;
    sub _build_fwd_con_from_clustalw_file;
    sub _build_rev_con_from_clustalw_file;
    sub print_unaligned_seqs_fasta;
    sub print_unaligned_seqs_fastq;
    sub print_con_fastq_str;
    sub _check_type_param;
    
    
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
        
        # Initialize the objects attributes
        $fwd_fastq_con_of{ident $new_obj} = undef;
        $rev_fastq_con_of{ident $new_obj} = undef;
        
        return $new_obj;
    }
    
    
    ###########
    # Setters #
    ###########
    sub set_fwd_fastq_con {
        my ($self, $con) = @_;
        
        $fwd_fastq_con_of{ident $self} = $con;
        return 1;
    }
    
    sub set_rev_fastq_con {
        my ($self, $con) = @_;
        
        $rev_fastq_con_of{ident $self} = $con;
        return 1;
    }
    
    
    ###########
    # Getters #
    ###########
    sub get_unaligned_seqs_fasta_str {
        my ($self, $type) = @_;
        
        # Make sure that $type is passed in correctly
        $type = lc $type;
        $self->_check_type_param($type, 'get_unaligned_seqs_fasta_str');
        
        my $str = '';
        my $temp_match_objs_href = $self->get_match_objs_href();
        
        # Loop throught the match objects of the given type and get their
        #   clustalw fasta str.
        if ( $type eq 'rev' ) {
            foreach my $id (keys %{$temp_match_objs_href} ) {
                $str .= $temp_match_objs_href->{$id}->get_rev_fasta_str();
            }
        }
        else {
            foreach  my $id (keys %{$temp_match_objs_href} ) {
                $str .= $temp_match_objs_href->{$id}->get_fwd_fasta_str();
            }
        }
        
        return $str;
    }
    
    sub get_unaligned_seqs_fastq_str {
        my ($self, $type) = @_;
        
        # Make sure that $type is passed in correctly
        $type = lc $type;
        $self->_check_type_param($type, 'get_unaligned_seqs_fastq_str');
        
        my $str = '';
        my $temp_match_objs_href = $self->get_match_objs_href();
        
        # Loop throught the match objects of the given type and get their
        #   fastq str.
        if ( $type eq 'rev' ){
            foreach my $id (keys %{$temp_match_objs_href} ) {
                $str .= $temp_match_objs_href->{$id}->get_rev_fastq_str();
            }
        }
        else {
            foreach my $id (keys %{$temp_match_objs_href} ) {
                $str .= $temp_match_objs_href->{$id}->get_fwd_fastq_str();
            }
        }
        
        return $str;
    }
    
    sub get_unaligned_quals_href {
        my ($self, $type) = @_;
        # I use this method when I parse the clustalw output.  As I am building the
        # consensus sequence I need to check the quality values in the event that
        # there is a tie between two bases.  If there is a tie I take the base with
        # the highest average quality score.
        
        # Make sure that $type is passed in correctly
        $type = lc $type;
        $self->_check_type_param($type, 'get_unaligned_quals_href');
        
        my $temp_match_objs_href = $self->get_match_objs_href();
        my %qualsHash = ();
        
        if ( $type eq 'rev' ) {
            foreach my $id ( keys %{$temp_match_objs_href} ) {
                $qualsHash{$id} = $temp_match_objs_href->{$id}->get_rev_quals_aref();
            }
        }
        else {
            foreach my $id ( keys %{$temp_match_objs_href} ) {
                $qualsHash{$id} = $temp_match_objs_href->{$id}->get_fwd_quals_aref();
            }
        }
        
        return \%qualsHash;
    }
    
    sub get_fwd_fastq_con {
        my ($self) = @_;
        
        return $fwd_fastq_con_of{ident $self};
    }
    
    sub get_rev_fastq_con {
        my ($self) = @_;
        
        return $rev_fastq_con_of{ident $self};
    }
    
    sub get_con_fastq_str {
        my ($self, $id, $type) = @_;
        
        # Make sure that $type is passed in correctly
        $type = lc $type;
        $self->_check_type_param($type, 'get_con_fastq_str');
        
        my $str;
        if ( $type eq 'fwd' ) {
                $str = $self->_get_fwd_con_fastq_str($id);
        }
        else {
            $str = $self->_get_rev_con_fastq_str($id);
        }
        
        return $str;
    }
    
    sub _get_fwd_con_fastq_str {
        my ($self, $id) = @_;
        
        my $temp_fwd_fastq_con = $fwd_fastq_con_of{ident $self};
        
        if ( ! defined $temp_fwd_fastq_con ) {
            MyX::Generic::Undef::Attribute->throw(
                                                  error => 'Undefined attribute value',
                                                  att_name => 'fwd_fastq_con_of',
                                                  );
        }
        if ( ! defined $temp_fwd_fastq_con->get_seq() ) {
            MyX::Generic::Undef->throw(
                                       error => 'Undefined rev sequence value',
                                       );
        }
        if ( ! defined $temp_fwd_fastq_con->get_quals_str() ) {
            MyX::Generic::Undef->throw(
                                       error => 'Undefined fwd sequence value',
                                       );
        }
        
        my $str = "@" . $id . '_' . $self->get_tag() . ";";
        $str .= "barcodelabel=$id;";
        $str .= "c_score=" . $temp_fwd_fastq_con->get_c_score() . ";";
        $str .= "depth=" . $self->get_seq_count() . "\n";
        $str .= $temp_fwd_fastq_con->get_seq() . "\n";
        $str .= "+\n";
        $str .= $temp_fwd_fastq_con->get_quals_str() . "\n";
        
        return $str;
    }
    
    sub _get_rev_con_fastq_str {
        my ($self, $id) = @_;
        
        my $temp_rev_fastq_con = $rev_fastq_con_of{ident $self};
        
        if ( ! defined $temp_rev_fastq_con ) {
            MyX::Generic::Undef::Attribute->throw(
                                                  error => 'Undefined attribute value',
                                                  att_name => 'rev_fastq_con_of',
                                                  );
        }
        if ( ! defined $temp_rev_fastq_con->get_seq() ) {
            MyX::Generic::Undef->throw(
                                       error => 'Undefined rev sequence value',
                                       );
        }
        if ( ! defined $temp_rev_fastq_con->get_quals_str() ) {
            MyX::Generic::Undef->throw(
                                       error => 'Undefined fwd sequence value',
                                       );
        }
        
        my $str = "@" . $id . '_' . $self->get_tag(). ";";
        $str .= "barcodelabel=$id;";
        $str .= "c_score=" . $temp_rev_fastq_con->get_c_score() . ";";
        $str .= "depth=" . $self->get_seq_count() . "\n";
        $str .= $temp_rev_fastq_con->get_seq() . "\n";
        $str .= "+\n";
        $str .= $temp_rev_fastq_con->get_quals_str() . "\n";
        
        return $str;
    }
    
    
    ##########
    # Others #
    ##########
    sub build_con {
        my ($self) = @_;
        
        my @fwd_seqs;
        my @rev_seqs;
        
        foreach my $id (keys %{$self->get_match_objs_href()} ) {
            push @fwd_seqs, $self->get_match_objs_href()->{$id}->get_fwd_amplicon();
        }
        foreach my $id (keys %{$self->get_match_objs_href()} ) {
            push @rev_seqs, $self->get_match_objs_href()->{$id}->get_rev_amplicon();
        }
        
        # Build and save the fastq consensus sequences
        my $fwd_fastq_con = build_consensus(\@fwd_seqs);
        $self->set_fwd_fastq_con($fwd_fastq_con);
        my $rev_fastq_con = build_consensus(\@rev_seqs);
        $self->set_rev_fastq_con($rev_fastq_con);
        
        # Return the fastq consensus sequence
        return ($fwd_fastq_con, $rev_fastq_con);
    }
    
    sub run_command_line_muscle {
        my ($self, $fasta_dir, $aln_dir) = @_;  # Fix the parameters
        
        # Make sure the fasta and aln directories already exist
        if ( ! -d $fasta_dir ) {
            MyX::Generic::DoesNotExist::Dir->throw(
                                                   error => 'Directory does not exist',
                                                   dir_name => $fasta_dir,
                                                   );
        }
        if ( ! -d $aln_dir ) {
            MyX::Generic::DoesNotExist::Dir->throw(
                                                   error => 'Directory does not exist',
                                                   dir_name => $aln_dir,
                                                   );
        }
        
        # Create fwd and rev fasta files to print the fwd and rev sequences to
        my $fwd_fasta_file = $fasta_dir . "/" . $self->get_tag() . "_fwd.fasta";
        my $rev_fasta_file = $fasta_dir . "/" . $self->get_tag() . "_rev.fasta";
        
        # Print the fwd and rev fasta files
        $self->print_unaligned_seqs_fasta($fwd_fasta_file, 'fwd');
        $self->print_unaligned_seqs_fasta($rev_fasta_file, 'rev');
        
        # Create fwd and rev output files for muscle to print the output to
        my $fwd_output_file = $aln_dir . "/" . $self->get_tag() . "_fwd.fasta";
        my $rev_output_file = $aln_dir . "/" . $self->get_tag() . "_rev.fasta";
        
        # Run clustalw through a system command on both fwd and rev fasta files
        my $fwd_command =   "muscle " .
                            "-in $fwd_fasta_file " .
                            "-out $fwd_output_file " .
                            "-diags " .
                            "-quiet ";
        my $rev_command =   "muscle " .
                            "-in $rev_fasta_file " .
                            "-out $rev_output_file " .
                            "-diags " .
                            "-quiet ";
        system($fwd_command);
        system($rev_command);
        
        return 1;
    }
    
    sub run_command_line_clustalw {
        my ($self, $fasta_dir, $aln_dir) = @_;  # Fix the parameters
        
        # Make sure the fasta and aln directories already exist
        if ( ! -d $fasta_dir ) {
            MyX::Generic::DoesNotExist::Dir->throw(
                                                   error => 'Directory does not exist',
                                                   dir_name => $fasta_dir,
                                                   );
        }
        if ( ! -d $aln_dir ) {
            MyX::Generic::DoesNotExist::Dir->throw(
                                                   error => 'Directory does not exist',
                                                   dir_name => $aln_dir,
                                                   );
        }
        
        # Create fwd and rev fasta files to print the fwd and rev sequences to
        my $fwd_fasta_file = $fasta_dir . "/" . $self->get_tag() . "_fwd.fasta";
        my $rev_fasta_file = $fasta_dir . "/" . $self->get_tag() . "_rev.fasta";
        
        # Print the fwd and rev fasta files
        $self->print_unaligned_seqs_fasta($fwd_fasta_file, 'fwd');
        $self->print_unaligned_seqs_fasta($rev_fasta_file, 'rev');
        
        # Create fwd and rev output files for clustalw to print the output to
        my $fwd_output_file = $aln_dir . "/" . $self->get_tag() . "_fwd.gde";
        my $rev_output_file = $aln_dir . "/" . $self->get_tag() . "_rev.gde";
        
        # Run clustalw through a system command on both fwd and rev fasta files
        my $fwd_command = "clustalw -infile=$fwd_fasta_file -output=gde -outorder=input -case=upper -outfile=$fwd_output_file -quiet -quicktree";
        my $rev_command = "clustalw -infile=$rev_fasta_file -output=gde -outorder=input -case=upper -outfile=$rev_output_file -quiet -quicktree";
        system($fwd_command);
        system($rev_command);
        
        return 1;
    }
    
    sub build_con_from_fasta_file {
        my ($self, $aln_dir) = @_;
        # This is a wrapper method to run both the fwd and rev consensus builders.
        
        # find the correct fasta alignment file
        my $aligned_seqs_fwd_file = $aln_dir . "/" . $self->get_tag() . "_fwd.fasta";
        my $aligned_seqs_rev_file = $aln_dir . "/" . $self->get_tag() . "_rev.fasta";
        
        # Note if the file does not exist errors are thrown from the ConsensusBuilder
        
        my $fwd_fastq_con = $self->_build_fwd_con_from_fasta_file(
            $aligned_seqs_fwd_file
        );
        my $rev_fastq_con = $self->_build_rev_con_from_fasta_file(
            $aligned_seqs_rev_file
        );
        
        return ($fwd_fastq_con, $rev_fastq_con);
    }
    
    sub _build_fwd_con_from_fasta_file {
        my ($self, $aligned_seqs_file) = @_;
        
        my $fwd_fastq_con = build_con_from_file(
                                $aligned_seqs_file,
                                "fasta",
                                $self->get_unaligned_quals_href('fwd')
                            );
        
        # Save the fwd fastq consensus sequence object
        $fwd_fastq_con_of{ident $self} = $fwd_fastq_con;
        
        # Return the fwd fastq consensus sequence object
        return $fwd_fastq_con;
    }
    
    sub _build_rev_con_from_fasta_file {
        my ($self, $aligned_seqs_file) = @_;
        
        my $rev_fastq_con = build_con_from_file(
                                $aligned_seqs_file,
                                "fasta",
                                $self->get_unaligned_quals_href('rev')
                            );
        
        # Save the rev fastq consensus sequence object
        $rev_fastq_con_of{ident $self} = $rev_fastq_con;
        
        # Return the rev fastq consensus sequence object
        return $rev_fastq_con;
    }
    
    sub build_con_from_clustalw_file {
        my ($self, $aln_dir) = @_;
        # This is a wrapper method to run both the fwd and rev consensus builders.
        
        # find the correct clustalw alignment file
        my $aligned_seqs_fwd_file = $aln_dir . "/" . $self->get_tag() . "_fwd.gde";
        my $aligned_seqs_rev_file = $aln_dir . "/" . $self->get_tag() . "_rev.gde";
        
        # Note if the file does not exist errors are thrown from the ConsensusBuilder
        
        my $fwd_fastq_con = $self->_build_fwd_con_from_clustalw_file($aligned_seqs_fwd_file);
        my $rev_fastq_con = $self->_build_rev_con_from_clustalw_file($aligned_seqs_rev_file);
        
        return ($fwd_fastq_con, $rev_fastq_con);
    }
    
    sub _build_fwd_con_from_clustalw_file {
        my ($self, $aligned_seqs_file) = @_;
        
        my $fwd_fastq_con = build_con_from_file(
                                $aligned_seqs_file,
                                "gde",
                                $self->get_unaligned_quals_href('fwd')
                            );
        
        # Save the fwd fastq consensus sequence object
        $fwd_fastq_con_of{ident $self} = $fwd_fastq_con;
        
        # Return the fwd fastq consensus sequence object
        return $fwd_fastq_con;
    }
    
    sub _build_rev_con_from_clustalw_file {
        my ($self, $aligned_seqs_file) = @_;
        
        my $rev_fastq_con = build_con_from_file(
                                $aligned_seqs_file,
                                "gde",
                                $self->get_unaligned_quals_href('rev')
                            );
        
        # Save the rev fastq consensus sequence object
        $rev_fastq_con_of{ident $self} = $rev_fastq_con;
        
        # Return the rev fastq consensus sequence object
        return $rev_fastq_con;
    }
    
    sub print_unaligned_seqs_fasta {
        my ($self, $file, $type) = @_;
        # NOTE: The user must be a little careful here.  They have to ensure
        # that the file type they pass in matches the file type of the file
        # name that they pass in.
        
        # Make sure that $type is passed in correctly
        $type = lc $type;
        $self->_check_type_param($type, 'print _unaligned_seqs_fasta');
        
        open (OUT, ">>$file") or die "Cannot open $file\nERROR: $!\n";
        
        if ( $type eq 'rev' ) {
            print OUT $self->get_unaligned_seqs_fasta_str('rev'); 
        }
        else {
            print OUT $self->get_unaligned_seqs_fasta_str('fwd'); 
        }
        
        close (OUT);
        
        return 1;
    }
    
    sub print_unaligned_seqs_fastq {
        my ($self, $file, $type) = @_;
        # NOTE: The user must be a little careful here.  They have to ensure
        # that the file type they pass in matches the file type of the file
        # name that they pass in.
        
        # Make sure that $type is passed in correctly
        $type = lc $type;
        $self->_check_type_param($type, 'print_unaligned_seqs_fastq');
        
        open (OUT, ">>$file") or die "Cannot open $file\nERROR: $!\n";
        
        if ( $type eq 'rev' ) {
            print OUT $self->get_unaligned_seqs_fastq_str('rev');
        }
        else {
            print OUT $self->get_unaligned_seqs_fastq_str('fwd');
        }
        
        close (OUT);
        
        return 1;
    }
    
    sub print_con_fastq_str {
        my ($self, $file, $id, $type) = @_;
        # NOTE: The user must be a little careful here.  They have to ensure
        # that the file type they pass in matches the file type of the file
        # name that they pass in.
        
        # Make sure that $type is passed in correctly
        $type = lc $type;
        $self->_check_type_param($type, 'print_con_fastq_str');
        
        open(OUT, ">>$file") or die "Cannot open $file\nERROR: $!\n";
        
        if ( $type eq 'rev' ) {
            print OUT $self->_get_con_fastq_str($id, 'rev');
        }
        else {
            print OUT $self->_get_con_fastq_str($id, 'fwd');
        }
        
        close(OUT);
        
        return 1;
    }

    sub _check_type_param {
        my ($self, $type, $method) = @_;
        
        # Make sure type is defined
        if ( ! defined $type ) {
            croak "type parameter not defined when calling $method";
        }
        
        # Make sure it is either 'fwd' or 'rev'
        $type = lc $type;
        my %allowed_types = (fwd => 1, rev => 1);
        if ( ! defined $allowed_types{$type} ) {
            croak "type paramter ($type) must be either \'fwd\' or \'rev\' when calling $method";
        }
        
        return 1;
    }
}


1;
__END__


#######
# POD #
#######
=head1 MTToolbox::MoleculeTagCategory::PE::NoOverlap

MTToolbox::MoleculeTagCategory::PE::NoOverlap - A class to store and operate on a set of paired end
(PE) sequences that do NOT have overlapping ends and that have the same molecule
tag

=head1 VERSION

This documentation refers to MTToolbox::MoleculeTagCategory::PE::NoOverlap version 4.1.2.

=head1 Included Modules

    Class::Std::Utils
    List::MoreUtils qw(any)
    Readonly
    Carp qw(carp croak)
    MyX::Generic 1.0.0
    BioUtils::ConsensusBuilder::ConsensusBuilder 1.0.0
    version
    MTToolbox::MoleculeTagCategory

=head1 Inherit

    MTToolbox::MoleculeTagCategory

=head1 SYNOPSIS
    
    use MTToolbox::MoleculeTagCategory::PE::NoOverlap;
    my $my_mt_category = MTToolbox::MoleculeTagCategory::PE::NoOverlap;->new(
                            "ATTGATGAG_ATGG"
                        );
    
    # Most important methods
    $my_mt_category->run_bio_perl_clustalw($factory);
    $my_mt_category->get_con_fastq_str("P1", 'fwd');
    $my_mt_category->get_con_fastq_str("P1", 'rev');

=head1 DESCRIPTION

MTToolbox::MoleculeTagCategory::PE::NoOverlap is a child class of
MTToolbox::MoleculeTagCategory.  It represents a set of sequences that do NOT
have overlapping ends and that have the same molecule tag.  

IMPORTANT: Because this object stores the forward and reverse sequences
serperately they are also processed seperately.  This required the main methods
in the parent class, MTToolbox::MoleculeTagCategory, to be overridden in this
class.  Many of those methods are wrapper methods that actually process the
forward and reverse reads seperately or processes only the type of read
specified by the $type parameter.

=head1 METHODS

=over

    new
    set_fwd_fastq_con
    set_rev_fastq_con
    get_unaligned_seqs_fasta_str
    get_unaligned_seqs_fastq_str
    get_unaligned_quals_href
    get_fwd_fastq_con
    get_rev_fastq_con
    get_con_fastq_str
    _get_fwd_con_fastq_str
    _get_rev_con_fastq_str
    build_con
    run_command_line_muscle
    run_command_line_clustalw
    build_con_from_fasta_file
    _build_fwd_con_from_fasta_file
    _build_rev_con_from_fasta_file
    build_con_from_clustalw_file
    _build_fwd_con_from_clustalw_file
    _build_rev_con_from_clustalw_file
    print_unaligned_seqs_fasta
    print_unaligned_seqs_fastq
    print_con_fastq_str
    _check_type_param
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: MTToolbox::MoleculeTagCategory::PE::NoOverlap->new($tag);
    Function: Creates a new MTToolbox::MoleculeTagCategory::PE::NoOverlap; object
    Returns: MTToolbox::MoleculeTagCategory::PE::::NoOverlap
    Args: -tag => a unique identifier representing the molecule tag(s) of a read
    Throws: MyX::Generic::Undef::Param
    Comments: NA
    See Also: MTToolbox::MoleculeTagCategory

=head2 set_fwd_fastq_con

    Title: set_fwd_fastq_con
    Usage: $my_mt_category->set_fwd_fastq_con($fwd_con);
    Function: Sets the fwd consensus sequence
    Returns: 1 on success
    Args: -con => BioUtils::ConsensusBuilder::FastqConsensus
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 set_rev_fastq_con

    Title: set_rev_fastq_con
    Usage: $my_mt_category->set_rev_fastq_con($rev_con);
    Function: Sets the rev consensus sequence
    Returns: 1 on success
    Args: -con => BioUtils::ConsensusBuilder::FastqConsensus
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_fwd_fastq_con

    Title: get_fwd_fastq_con
    Usage: my $con = $my_mt_category->get_fwd_fastq_con();
    Function: Gets the fwd fastq consensus sequence object
    Returns: BioUtils::ConsensusBuilder::FastqConsensus
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_rev_fastq_con

    Title: get_rev_fastq_con
    Usage: my $con = $my_mt_category->get_rev_fastq_con();
    Function: Gets the rev fastq consensus sequence object
    Returns: BioUtils::ConsensusBuilder::FastqConsensus
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_unaligned_seqs_fasta_str

    Title: get_unaligned_seqs_fasta_str
    Usage: $my_mt_category->get_unaligned_seqs_fasta_str($type);
    Function: Gets the string of unaligned forward sequences in clustalw format
    Returns: String
    Args: -type => 'fwd' or 'rev'
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_unaligned_seqs_fastq_str

    Title: get_unaligned_seqs_fastq_str
    Usage: $my_mt_category->get_unaligned_seqs_fastq_str($type);
    Function: Gets the string of unaligned forward sequences in FASTQ format
    Returns: String
    Args: -type => 'fwd' or 'rev'
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_unaligned_quals_href

    Title: get_unaligned_quals_href
    Usage: $my_mt_category->get_unaligned_quals_href($type);
    Function: Gets the quality values of unaligned forward sequences
    Returns: Hash reference -- KEY => clustalwId, VALUE => forward quals array ref
    Args: -type => 'fwd' or 'rev'
    Throws: NA
    Comments: NA
    See Also: Match::PE::NoOverlap::getQualArrRef
    
=head2 get_con_fastq_str

    Title: get_con_fastq_str
    Usage: $my_mt_category->get_con_fastq_str($id, $type);
    Function: Gets the string of consensus sequences in FASTQ format
    Returns: String
    Args: -id => The sample id to be printed in the header
          -type => 'fwd' or 'rev'
    Throws: MyX::Generic::Undef,
            MyX::Generic::Undef::Attribute,
    Comments: NA
    See Also: _get_fwd_con_fastq_str
              _get_rev_con_fastq_str
    
=head2 _get_fwd_con_fastq_str

    Title: _get_fwd_con_fastq_str
    Usage: $my_mt_category->_get_fwd_con_fastq_str($id);
    Function: Gets the string of forward consensus sequences in FASTQ format
    Returns: String
    Args: -id => Sample ID to be printed in the header
    Throws: MyX::Generic::Undef,
            MyX::Generic::Undef::Attribute,
    Comments: NA
    See Also: NA

=head2 _get_rev_con_fastq_str

    Title: _get_rev_con_fastq_str
    Usage: $my_mt_category->_get_rev_con_fastq_str($id);
    Function: Gets the string of reverse consensus sequences in FASTQ format
    Returns: String
    Args: -id => Sample ID to be printed in the header
    Throws: MyX::Generic::Undef,
            MyX::Generic::Undef::Attribute,
    Comments: NA
    See Also: NA
    
=head2 build_con

    Title: build_con
    Usage: $my_mt_category->build_con();
    Function: Runs ConsensusBuilder methods to build a consensus
              sequence using the unaligned sequences in this MT
    Returns: (fwd_FastqConsensus, rev_FastqConsensus)
    Args: NA
    Throws: NA
    Comments: This was a method that was added in later editions.  It is
              possible that for accurate Illumina reads with few INDELs it
              is unnecessary to take the time build a MSA.
    See Also: BioUtils::ConsensusBuilder::FastqConsensus

=head2 run_command_line_muscle

    Title: run_command_line_muscle
    Usage: $my_mt_category->run_command_line_muscle($factory);
    Function: Runs muscle through the command line to build an MSA using all
              sequences in the molecule tag category
    Returns: 1 on successful completion
    Args: -fasta_dir => A directory of fasta input files
          -aln_dir => A directory to put the output files
    Throws: MyX::Generic::DoesNotExist::Dir
    Comments: Both directories must be created before running this method.  This
              method is reommended over run_command_line_clustalw because it is
              faster and more accurate.
    See Also: NA

=head2 run_command_line_clustalw

    Title: run_command_line_clustalw
    Usage: $my_mt_category->run_command_line_clustalw($factory);
    Function: Runs clustalw through the command line to build an MSA using all
              sequences in the molecule tag category
    Returns: 1 on successful completion
    Args: -fasta_dir => A directory of fasta input files
          -aln_dir => A directory to put the output files
    Throws: MyX::Generic::DoesNotExist::Dir
    Comments: NA
    See Also: NA
    
=head2 build_con_from_fasta_file

    Title: build_con_from_fasta_file
    Usage: $my_mt_category->build_con_from_fasta_file($aln_dir);
    Function: Calls the consensus building methods for forward and reverse MSAs
    Returns: ($fwd_fastq_consensus, $rev_fastq_consensus)
    Args: -aln_dir => a directory path containing the fasta MSA output files
    Throws: NA
    Comments: Because this method is named the same as in other Sample objects
              it allows any sample to build consensus sequences no matter what
              kind of Sample it is.  To build the forward and reverse consensus
              sequences seperately use _build_fwd_con_from_fasta_file and
              _build_rev_con_from_fasta_file.  Also note that this method, unlike
              the other build_con_from_fasta_file methods in the Sample objects
              returns an array of two fastqConsensus objects. The first
              fastqConsensus object is the forward consensus and the second is
              the reverse consensus.
    See Also: BioUtils::FastqConsensus,
              MTToolbox::MoleculeTagCategory::NoOverlap::run_command_line_muscle

=head2 _build_fwd_con_from_fasta_file

    Title: _build_fwd_con_from_fasta_file
    Usage: $my_mt_category->_build_fwd_con_from_fasta_file($aligned_seqs_file);
    Function: Runs ConsensusBuilder methods to build a forward consensus
              sequence using the file generate from fasta MSA
    Returns: FastqConsensus
    Args: -aligned_seqs_file => the fasta output file from a MSA
    Throws: NA
    Comments: NA
    See Also: BioUtils::FastqConsensus,
              MTToolbox::MoleculeTagCategory::NoOverlap::run_command_line_muscle

=head2 _build_rev_con_from_fasta_file

    Title: _build_rev_con_from_fasta_file
    Usage: $my_mt_category->_build_rev_con_from_fasta_file($aligned_seqs_file);
    Function: Runs ConsensusBuilder methods to build a reverse consensus
              sequence using the file generate from fasta MSA
    Returns: FastqConsensus
    Args: -aligned_seqs_file => the gde output file from a MSA
    Throws: NA
    Comments: NA
    See Also: BioUtils::FastqConsensus,
              MTToolbox::MoleculeTagCategory::NoOverlap::run_command_line_muscle

=head2 build_con_from_clustalw_file

    Title: build_con_from_clustalw_file
    Usage: $my_mt_category->build_con_from_clustalw_file($aln_dir);
    Function: Calls the consensus building methods for forward and reverse MSAs
    Returns: ($fwd_fastq_consensus, $rev_fastq_consensus)
    Args: -aln_dir => a directory path containing the clustalw gde output files
    Throws: NA
    Comments: Because this method is named the same as in other Sample objects
              it allows any sample to build consensus sequences no matter what
              kind of Sample it is.  To build the forward and reverse consensus
              sequences seperately use _build_fwd_con_from_clustalw_file and
              _build_rev_con_from_clustalw_file.  Also note that this method, unlike
              the other build_con_from_clustalw_file methods in the Sample objects
              returns an array of two fastqConsensus objects. The first
              fastqConsensus object is the forward consensus and the second is
              the reverse consensus.
    See Also: BioUtils::FastqConsensus,
              MTToolbox::MoleculeTagCategory::NoOverlap::run_command_line_clustalw

=head2 _build_fwd_con_from_clustalw_file

    Title: _build_fwd_con_from_clustalw_file
    Usage: $my_mt_category->_build_fwd_con_from_clustalw_file($aligned_seqs_file);
    Function: Runs ConsensusBuilder methods to build a forward consensus
              sequence using the file generate from commandline clustalw
    Returns: FastqConsensus
    Args: -aligned_seqs_file => the gde output file from commandline clustalw
    Throws: NA
    Comments: NA
    See Also: BioUtils::FastqConsensus,
              MTToolbox::MoleculeTagCategory::NoOverlap::run_command_line_clustalw

=head2 _build_rev_con_from_clustalw_file

    Title: _build_rev_con_from_clustalw_file
    Usage: $my_mt_category->_build_rev_con_from_clustalw_file($aligned_seqs_file);
    Function: Runs ConsensusBuilder methods to build a reverse consensus
              sequence using the file generate from commandline clustalw
    Returns: FastqConsensus
    Args: -aligned_seqs_file => the gde output file from commandline clustalw
    Throws: NA
    Comments: NA
    See Also: BioUtils::FastqConsensus,
              MTToolbox::MoleculeTagCategory::NoOverlap::run_command_line_clustalw

=head2 print_unaligned_seqs_fasta

    Title: print _unaligned_seqs_fasta
    Usage: $my_mt_category->print _unaligned_seqs_fasta($file, $type);
    Function: Print the unaligned seqeunces in clustalw format
    Returns: 1 on successful compeletion
    Args: -file => path to file for printing
          -type => 'fwd' or 'rev'
    Throws: NA
    Comments: Be careful to make sure the type parameter matches the type of
              sequence that are in the file that you are printing to.  For
              example, something like this will probably lead to headaches:
              print _unaligned_seqs_fasta("rev_seqs.fasta", 'fwd');
    See Also: NA

=head2 print_unaligned_seqs_fastq

    Title: print_unaligned_seqs_fastq
    Usage: $my_mt_category->print_unaligned_seqs_fastq($file, $type);
    Function: Print the unaligned forward seqeunces in FASTQ format
    Returns: 1 on successful compeletion
    Args: -file => path to file for printing
          -type => 'fwd' or 'rev'
    Throws: NA
    Comments: Be careful to make sure the type parameter matches the type of
              sequence that are in the file that you are printing to.  For
              example, something like this will probably lead to headaches:
              print_unaligned_seqs_fastq("rev_seqs.fastq", 'fwd');
    See Also: NA

=head2 print_con_fastq_str

    Title: print_con_fastq_str
    Usage: $my_mt_category->print_con_fastq_str($file, $id, $type);
    Function: Print the consensus seqeunces in FASTQ format
    Returns: 1 on successful compeletion
    Args: -file => path to file for printing
          -id => Sample ID to be printed in the header
          -type => 'fwd' or 'rev'
    Throws: NA
    Comments: Be careful to make sure the type parameter matches the type of
              sequence that are in the file that you are printing to.  For
              example, something like this will probably lead to headaches:
              print_con_fastq_str("rev_cons.fastq", 'fwd');
    See Also: NA

=head2 _check_type_param

    Title: _check_type_param
    Usage: _check_type_param($type);
    Function: Makes sure that $type is defined and set to either 'fwd' or 'rev'.
              If this is not the case it throws the appropriate error.
    Returns: 1 on successful compeletion
    Args: -type => 'fwd' or 'rev'
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
