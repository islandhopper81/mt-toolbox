package MTToolbox::MoleculeTagCategory;

use strict;
use warnings;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(carp croak);
use MyX::Generic 0.0.3;
use BioUtils::ConsensusBuilder::ConsensusBuilder 1.0.0 qw(build_con_from_file build_from_clustalw_file build_consensus);
use version; our $VERSION = qv('4.1.2');


{
    Readonly my $NEW_USAGE => q{ new( {tag => ,
                                       } ) };
    
    # Attributes #
    my %tag_of;
    my %seq_count_of;
    my %match_objs_href_of;
    my %fastq_con_of;
    my %is_con_defined_of;
    
    # Setters #
    sub set_tag;
    sub set_fastq_con;
    sub set_is_con_defined;
    
    # Getters #
    sub get_tag;
    sub get_seq_count;
    sub get_match_objs_href;
    sub get_unaligned_seqs_fasta_str;
    sub get_unaligned_seqs_fastq_str;
    sub get_unaligned_quals_href;
    sub get_con_fastq_str;
    sub is_con_defined;
    
    # Others #
    sub add_match;
    sub build_con; # builds a consensus sequence w/o building a MSA
    sub run_command_line_muscle;
    sub build_con_from_fasta_file; # muscle's output is fasta
    sub run_command_line_clustalw;
    sub build_con_from_clustalw_file;
    sub print_unaligned_seqs_fasta;
    sub print_unaligned_seqs_fastq;
    sub print_con_fastq_str;
    
    
    
    
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
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        # Initialize the objects attributes
        $tag_of{ident $new_obj} = $arg_href->{tag};
        $seq_count_of{ident $new_obj} = 0;
        my %matchedFastqSeqObjs = (); # A hash ref with KEY: the match ID,  VALUE: match object
        $match_objs_href_of{ident $new_obj} = \%matchedFastqSeqObjs;
        $is_con_defined_of{ident $new_obj} = 0;  # starts as undefined

    
        return $new_obj;
    }
    
    
    ###########
    # Setters #
    ###########
    sub set_tag {
        my ($self, $tag) = @_;
        $tag_of{ident $self} = $tag;
        return 1;
    }
    
    sub set_fastq_con {
        my ($self, $fastq_con) = @_;
        $fastq_con_of{ident $self} = $fastq_con;
        return 1;
    }
    
    sub set_is_con_defined {
        my ($self, $bool) = @_;
        $is_con_defined_of{ident $self} = $bool;
        return 1;
    }
    
    
    
    ###########
    # Getters #
    ###########
    sub get_tag {
        my ($self) = @_;
        return $tag_of{ident $self};
    }
    
    sub get_seq_count {
        my ($self) = @_;
        return $seq_count_of{ident $self};
    }
    
    sub get_match_objs_href {
        my ($self) = @_;
        
        return $match_objs_href_of{ident $self};
    }
    
     sub get_unaligned_seqs_fasta_str {
        my ($self) = @_;
        my $str = '';
        
        # Loop through all the matchObjs and get their clustalw fasta str
        foreach my $id (keys %{$self->get_match_objs_href()} ) {
            $str .= $self->get_match_objs_href()->{$id}->get_fasta_str();
        }
        
        return $str;
    }
    
    sub get_unaligned_seqs_fastq_str {
        my ($self) = @_;
        my $str = '';
        
        # Loop through all the matchObjs and get their fastqStr
        foreach my $id (keys %{$self->get_match_objs_href()} ) {
            $str .= $self->get_match_objs_href()->{$id}->get_fastq_str();
        }
        
        return $str;
    }
    
    sub get_unaligned_quals_href {
        my ($self) = @_;
        # I use this method when I parse the clustalw output.  As I am building the
        # consensus sequence I need to check the quality values in the event that
        # there is a tie between two bases.  If there is a tie I take the base with
        # the highest average quality score.
        
        my %qualsHash = ();
        my $temp_match_objs_href = $self->get_match_objs_href();
        foreach my $id ( keys %{$temp_match_objs_href} ) {
            #my $clustalwId = $temp_match_objs_href->{$id}->get_clustalw_id();
            $qualsHash{$id} = $temp_match_objs_href->{$id}->get_quals_aref();
        }
        
        return \%qualsHash;
    }
    
    sub get_con_fastq_str {
        my ($self, $id) = @_;
        my $temp_fastq_con = $fastq_con_of{ident $self};
        
        # Make sure the consensus sequence has been build and saved in this object
        # In the future I want to create error objects for these so I can handle them
        # in the caller method
        if ( ! defined $temp_fastq_con ) {
            carp "Consensus is not defined";
            return '';
        }
        if ( ! defined $temp_fastq_con->get_seq() ) {
            carp "Consensus sequence is not defined";
            return '';
        }
        if ( ! defined $temp_fastq_con->get_quals_str() ) {
            carp "Consensus quality values are not defined";
            return '';
        }
        
        my $str = "@" . $id . '_' . $self->get_tag() . ";";
        $str .= "barcodelabel=$id;";
        $str .= "c_score=" . $temp_fastq_con->get_c_score() . ";";
        $str .= "depth=" . $self->get_seq_count() . "\n";
        $str .= $temp_fastq_con->get_seq() . "\n";
        $str .= "+\n";
        $str .= $temp_fastq_con->get_quals_str() . "\n";
        
        return $str;
    }
    
    sub is_con_defined {
        my ($self) = @_;
        
        # return the boolean if it is defined.
        if ( my $val = $is_con_defined_of{ident $self} ) {
            return $val;
        }
        
        return 0;
    }
    
    
    ##########
    # Others #
    ##########
    sub add_match {
        my ($self, $match) = @_;
        
        $match_objs_href_of{ident $self}{$match->get_id()} = $match;
        $seq_count_of{ident $self}++;
        
        return 1;
    }
    
    sub build_con {
        my ($self) = @_;
        
        # builds a consensus sequence without building a MSA
        
        my @seqs;
        foreach my $id (keys %{$self->get_match_objs_href()} ) {
            push @seqs, $self->get_match_objs_href()->{$id}->get_fwd_amplicon();
        }
        
        # Build the fastq consensus sequence if possible
        my $fastq_con;
        eval {
            $fastq_con = build_consensus(\@seqs);
        };
        # catch errors
        if ( my $e = Exception::Class->caught(
                'BioUtils::MyX::ConsensusBuilder::TooFewSeqs') ) {
            
            carp "WARNING: " . $e->error() . ": " . $self->get_tag();
            $self->set_is_con_defined(0); # set to false
            return undef;
        }
        
        # save the fastq consensus sequence
        $self->set_fastq_con($fastq_con);
        $self->set_is_con_defined(1);  # set to TRUE
        
        # Return the fastq consensus sequence
        return $fastq_con;
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
        
        my $fastaFile = $fasta_dir . "/" . $self->get_tag() . ".fasta";
        my $outputFile = $aln_dir . "/" . $self->get_tag() . ".fasta";
        
        $self->print_unaligned_seqs_fasta($fastaFile);
        
        my $command = "muscle " .
                      "-in $fastaFile " .
                      "-out $outputFile " .
                      "-diags " .
                      "-quiet ";
        #print $command . "\n";
        system($command);
        
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
        
        my $fastaFile = $fasta_dir . "/" . $self->get_tag() . ".fasta";
        my $outputFile = $aln_dir . "/" . $self->get_tag() . ".gde";
        
        $self->print_unaligned_seqs_fasta($fastaFile);
        
        my $command = "clustalw -infile=$fastaFile " .
                      "-output=gde " .
                      "-outorder=input " .
                      "-case=upper " .
                      "-outfile=$outputFile " .
                      "-quiet " .
                      "-quicktree";
        #print $command . "\n";
        system($command);
        
        return 1;
    }
    
    sub build_con_from_fasta_file {
        my ($self, $aln_dir) = @_;
        
        # Find the gde file for this molecule tag
        my $aligned_seqs_file = $aln_dir . "/" . $self->get_tag() . ".fasta";
        if ( ! -f $aligned_seqs_file ) {
            croak "Missing aligned_seqs_file from the Muscle output";
        }
        
        # Build the fastq consensus sequence
        my $fastq_con = build_con_from_file(
                            $aligned_seqs_file,
                            "fasta",
                            $self->get_unaligned_quals_href()
                        );
        
        # Save the fastq consensus sequence
        $fastq_con_of{ident $self} = $fastq_con;
        
        # Return the fastq consensus sequence
        return $fastq_con;
    }
    
    sub build_con_from_clustalw_file {
        my ($self, $aln_dir) = @_;
        
        # Find the gde file for this molecule tag
        my $aligned_seqs_file = $aln_dir . "/" . $self->get_tag() . ".gde";
        if ( ! -f $aligned_seqs_file ) {
            croak "Missing aligned_seqs_file from the clustalw output";
        }
        
        # Build the fastq consensus sequence
        my $fastq_con = build_con_from_file(
                            $aligned_seqs_file,
                            "gde",
                            $self->get_unaligned_quals_href()
                        );
        
        # Save the fastq consensus sequence
        $fastq_con_of{ident $self} = $fastq_con;
        
        # Return the fastq consensus sequence
        return $fastq_con;
    }
    
    sub print_unaligned_seqs_fasta {
        my ($self, $file) = @_;
        
        open (OUT, ">>$file") or croak "Cannot open $file\nERROR: $!\n";
        
        print OUT $self->get_unaligned_seqs_fasta_str();
        
        close (OUT);
        
        return 1;
    }
    
    sub print_unaligned_seqs_fastq {
        my ($self, $file) = @_;
        
        open (OUT, ">>$file") or croak "Cannot open $file\nERROR: $!\n";
        
        print OUT $self->get_unaligned_seqs_fastq_str();
        
        close (OUT);
        
        return 1;
    }
    
    sub print_con_fastq_str {
        my ($self, $file, $id) = @_;
        
        # throw an error if fastq_con is undefied
        if ( ! $self->is_con_defined() ) {
            MyX::Generic::Undef::Attribute->throw(
                att_name => "fastq_con",
                error => "Undefined Fastq Consensus for " . $self->get_tag(),
            );
        }
        
        open(OUT, ">>$file") or croak "Cannot open $file\nERROR: $!\n";
        
        print OUT $self->get_con_fastq_str($id);
            
        close(OUT);
        
        return 1;
    }
}


1;
__END__

#######
# POD #
#######
=head1 MTToolbox::MoleculeTagCategory

MTToolbox::MoleculeTagCategory - A class to store and operate on a set of
sequences that have the same molecule tag

=head1 VERSION

This documentation refers to MTToolbox::MoleculeTagCategory version 4.1.2.

=head1 Included Modules

    Class::Std::Utils
    List::MoreUtils qw(any)
    Readonly
    Carp qw(carp croak)
    MyX::Generic 1.0.0
    BioUtils::ConsensusBuilder::ConsensusBuilder 1.0.0 qw(build_con_from_file build_consensus buildFromClustalwFile)
    version

=head1 Inherit

    NA

=head1 SYNOPSIS
    
    use MTToolbox::MoleculeTagCategory;
    my $my_mt_category = MTToolbox::MoleculeTagCategory->new(
                                        {tag => "ATTGATGAG_ATGG"}
                                    );
    
    # Attribute Getters
    my $tag = $my_mt_category->get_tag();
    my $seq_count = $my_mt_category->get_seq_count();
    my $match_objs_href = $my_mt_category->get_match_objs_href();
    my $fasta_str = $my_mt_category->get_unaligned_seqs_fasta_str();
    my $fastq_str = $my_mt_category->get_unaligned_seqs_fastq_str();
    my $quals_href = $my_mt_category->get_unaligned_quals_href();
    my $con_fastq_str = $my_mt_category->get_con_fastq_str();
    my $is_con_defined = $my_mt_category->is_con_defined();
    
    # Attribute Setters
    $my_mt_category->set_tag("ATTGATGAG_ATGC");
    $my_mt_category->set_fastq_con($fastq_con_obj);
    $my_mt_category->set_is_con_defined(0);
    
    # Others
    $my_mt_category->add_match($match_obj);
    $my_mt_category->build_con();
    $my_mt_category->run_command_line_muscle($fasta_dir, $aln_dir);
    $my_mt_category->run_command_line_clustalw($fasta_dir, $aln_dir);
    $my_mt_category->build_con_from_fasta_file($aln_dir);
    $my_mt_category->build_con_from_clustalw_file($aln_dir);
    $my_mt_category->print_unaligned_seqs_fasta($out_file);
    $my_mt_category->print_unaligned_seqs_fastq($out_file);
    $my_mt_category->print_con_fastq_str($out_file, $id);
    

=head1 DESCRIPTION

MoleculeTagCategory is an abstract class.  It and its repective children classes
hold a set of sequences that have the same molecule tag.  There are also
opperations that can be run on the set of sequences such as building a multiple
sequence alignment (MSA) and building a consensus sequence.

A set of sequences in a molecule tag category represent one DNA molecule that has
been PCRed and sequenced.  The goal is to take the set of sequences and collapse
it down to one sequence representing the original DNA molecule.  This process
reduces PCR and sequencing errors.

There are three main steps to this process.  First, we have to categorize the
raw reads by their molecule tags.  This operation is done in the Sample object
and its various children classes.  A Sample object can create a
MoleculeTagCategory and add sequences to it by calling add_match($match_obj).
Second, the set of sequences in a MoleculeTagCategory are used to build a
multiple sequence alignment (MSA).  Third, a consensus sequence is build from
the MSA and the quality values of the orginal raw reads.

In a previous version bioperl was used to build MSAs.  This was deamed inferior
to the commandline version of clustalw because the bioperl packages are simply
a wrapper for commandline clustalw and runs much slower.

=head1 METHODS

=over

    new
    set_tag
    set_fastq_con
    set_is_con_defined
    get_tag
    get_seq_count
    get_match_objs_href
    get_unaligned_seqs_fasta_str
    get_unaligned_seqs_fastq_str
    get_unaligned_quals_href
    get_con_fastq_str
    is_con_defined
    add_match
    build_con
    run_command_line_muscle
    build_con_from_fasta_file
    run_command_line_clustalw
    build_con_from_clustalw_file
    print_unaligned_seqs_fasta
    print_unaligned_seqs_fastq
    print_con_fastq_str
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: MTToolbox::MoleculeTagCategory->new({tag => $tag});
    Function: Creates a new MTToolbox::MoleculeTagCategory object
    Returns: MTToolbox::MoleculeTagCategory
    Args: -tag => a unique identifier representing the molecule tag(s) of a read
    Throws: MyX::Generic::Undef::Param
    Comments: NA
    See Also: NA

=head2 set_tag

    Title: set_tag
    Usage: $my_mt_category->set_tag($tag);
    Function: Sets the tag value
    Returns: 1 on successful completion
    Args: -tag => a unique identifier representing the molecule tag(s) of a read
    Throws: NA
    Comments: This NOT a raw read id.  The format of the tag depends on the types
              of reads that are input.  There are single end reads that have only
              a forward tag.  Paired end reads have a forward and reverse tag
              seperated by an underscore.
    See Also: MTToolbox::MoleculeTagCategory::NoOverlap,
              MTToolbox::MoleculeTagCategory::Overlap
              
=head2 set_fastq_con

    Title: set_fastq_con
    Usage: $my_mt_category->set_fastq_con($fastq_con);
    Function: Sets the consensus sequence attribute
    Returns: 1 on successful completion
    Args: -fastq_con => an object of type
                        BioUtils::ConsensusBuilder::FastqConsensus
    Throws: NA
    Comments: This methods should be used with caution.  It is primarily
              intended to be used internally.  It is dangerous to set the
              consensus sequence without generating it from the reads that are
              stored within the MoleculeTagCategory object.  I included this
              method for advanced users who may with to write their own
              consensus building algorithms and subsiquently set the consensus
              sequence attribute after running those aglorthims.
    See Also: BioUtils::ConsensusBuilder::FastqConsensus
    
=head2 set_is_con_defined

    Title: set_is_con_defined
    Usage: $my_mt_category->set_is_con_defined($bool);
    Function: Sets the boolean value indicating if the consensus seq is defined
    Returns: 1 on successful completion
    Args: -bool => a boolean value (1 => TRUE, 2 => FALSE) indicating if the
                   fastq consensus has been defined.
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_tag

    Title: get_tag
    Usage: $my_mt_category->get_tag();
    Function: Gets the tag value
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_seq_count

    Title: get_seq_count
    Usage: $my_mt_category->get_seq_count();
    Function: Gets the number of sequence in the molecule tag category
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_match_objs_href

    Title: get_match_objs_href
    Usage: $my_mt_category->get_match_objs_href();
    Function: Gets the hash reference for the data structure of match objects
    Returns: hash ref
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_unaligned_seqs_fasta_str

    Title: get_unaligned_seqs_fasta_str
    Usage: $my_mt_category->get_unaligned_seqs_fasta_str();
    Function: Gets the string of unaligned sequences in clustalw format
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_unaligned_seqs_fastq_str

    Title: get_unaligned_seqs_fastq_str
    Usage: $my_mt_category->get_unaligned_seqs_fastq_str();
    Function: Gets the string of unaligned sequences in FASTQ format
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_unaligned_quals_href

    Title: get_unaligned_quals_href
    Usage: $my_mt_category->get_unaligned_quals_href();
    Function: Gets the quality values of unaligned sequences
    Returns: Hash reference -- KEY => clustalwId, VALUE => quals array ref
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_con_fastq_str

    Title: get_con_fastq_str
    Usage: $my_mt_category->get_con_fastq_str($id);
    Function: Gets the string of consensus sequences in FASTQ format
    Returns: String
    Args: -id => The sample id to be printed in the header
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 is_con_defined

    Title: is_con_defined
    Usage: $my_mt_category->is_con_defined();
    Function: Gets the boolean indicating if the fastq consensus seq has been
              defined.
    Returns: bool
    Args: NA
    Throws: NA
    Comments: 1 => TRUE, 2 => FALSE
    See Also: NA

=head2 add_match

    Title: add_match
    Usage: $my_mt_category->add_match($match_obj);
    Function: Adds a match object with the same tag to this category
    Returns: 1 on successful completion
    Args: -match_obj => a match object (or one of its child classes)
                        representing a raw read
    Throws: NA
    Comments: NA
    See Also: Match
    
=head2 build_con

    Title: build_con
    Usage: $my_mt_category->build_con();
    Function: Runs ConsensusBuilder methods to build a consensus
              sequence using the unaligned sequences in this MT
    Returns: FastqConsensus
    Args: NA
    Throws: NA
    Comments: This was a method that was added in later editions.  It is
              possible that for accurate Illumina reads with few INDELs it
              is unnecessary to take the time build a MSA.
    See Also: BioUtils::ConsensusBuilder::FastqConsensus
    
=head2 run_command_line_muscle

    Title: run_command_line_muscle
    Usage: $my_mt_category->run_command_line_muscle($fasta_dir, $aln_dir);
    Function: Runs the command line version of muscle
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
    Usage: $my_mt_category->run_command_line_clustalw($fasta_dir, $aln_dir);
    Function: Runs the command line version of clustalw
    Returns: 1 on successful completion
    Args: -fasta_dir => A directory of fasta input files
          -aln_dir => A directory to put the output files
    Throws: MyX::Generic::DoesNotExist::Dir
    Comments: Both directories must be created before running this method
    See Also: NA
    
=head2 build_con_from_fasta_file

    Title: build_con_from_clustalw_file
    Usage: $my_mt_category->build_con_from_clustalw_file($aln_dir);
    Function: Runs ConsensusBuilder methods to build a consensus
              sequence using the file generate from commandline clustalw
    Returns: FastqConsensus
    Args: -aln_dir => a path to the directory with gde clustalw output files
    Throws: NA
    Comments: NA
    See Also: BioUtils::FastqConsensus,
              MTToolbox::MoleculeTagCategory::SE::runCommandlineClustalw

=head2 build_con_from_clustalw_file

    Title: build_con_from_clustalw_file
    Usage: $my_mt_category->build_con_from_clustalw_file($aln_dir);
    Function: Runs ConsensusBuilder methods to build a consensus
              sequence using the file generate from commandline clustalw
    Returns: FastqConsensus
    Args: -aln_dir => a path to the directory with gde clustalw output files
    Throws: NA
    Comments: NA
    See Also: BioUtils::FastqConsensus,
              MTToolbox::MoleculeTagCategory::SE::runCommandlineClustalw

=head2 print_unaligned_seqs_fasta

    Title: print_unaligned_seqs_fasta
    Usage: $my_mt_category->print_unaligned_seqs_fasta($file);
    Function: Print the unaligned seqeunces in clustalw format
    Returns: 1 on successful compeletion
    Args: -file => path to file for printing
    Throws: NA
    Comments: NA
    See Also: NA

=head2 print_unaligned_seqs_fastq

    Title: print_unaligned_seqs_fastq
    Usage: $my_mt_category->print_unaligned_seqs_fastq($file);
    Function: Print the unaligned seqeunces in FASTQ format
    Returns: 1 on successful compeletion
    Args: -file => path to file for printing
    Throws: NA
    Comments: NA
    See Also: NA

=head2 print_con_fastq_str

    Title: print_con_fastq_str
    Usage: $my_mt_category->print_con_fastq_str($file, $id);
    Function: Print the consensus seqeunces in FASTQ format
    Returns: 1 on successful compeletion
    Args: -file => path to file for printing
          -id => Sample id to be printed in the header
    Throws: MyX::Generic::Undef::Attribute
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



