package MTToolbox::Sample;

use strict;
use warnings;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Math::BigInt;
use Carp qw(carp croak);
use Chart::Gnuplot;
use version; our $VERSION = qv('4.1.2');

use MTToolbox::Mismatch;
use MTToolbox::MoleculeTagCategory;
use MTToolbox::SampleSummary;
use MTToolbox::MTDepthDist;
use BioUtils::FastqSeq 1.0.0;
use BioUtils::FastqIO 1.0.0;
use BioUtils::QC::FastqFilter 1.0.7;
use MyX::Generic;


{
    Readonly my $NEW_USAGE => q{ new( {id => ,
                                       name => ,
                                       barcode => ,
                                       fwd_linker => ,
                                       fwd_primer => ,
                                       rev_primer => ,
                                       fwd_mer_len => ,
                                       fwd_max_shifts => ,
                                       min_con_depth => ,
                                       diginorm_max => ,
                                       output_dir => ,
                                       } ) };
    
    Readonly my $DIGI_MAX => 1000000;

    # Attributes #
    my %id_of;
    my %name_of;
    my %barcode_of;
    my %fwd_linker_of;
    my %fwd_primer_of;
    my %rev_primer_of;
    my %fwd_mer_len_of;
    my %fwd_max_shifts_of;
    my %min_con_depth_of;
    my %diginorm_max_of;
    my %output_dir_of;
    my %mismatches_of;
    my %MT_objs_href_of;
    my %summ_info_obj_of;
    
    # Setters #
    sub set_id;
    sub set_name;
    sub set_barcode;
    sub set_fwd_linker;
    sub set_fwd_primer;
    sub set_rev_primer;
    sub set_fwd_mer_len;
    sub set_fwd_max_shifts;
    sub set_min_con_depth;
    sub set_diginorm_max;
    sub set_output_dir;
    sub _set_summ_info_obj;
    
    # Getters #
    sub get_id;
    sub get_name;
    sub get_barcode;
    sub get_fwd_linker;
    sub get_fwd_primer;
    sub get_rev_primer;
    sub get_fwd_mer_len;
    sub get_fwd_max_shifts;
    sub get_min_con_depth;
    sub get_diginorm_max;
    sub get_mismatch_seqs_aref;
    sub get_MT_objs_href;
    sub get_output_dir;
    sub get_summ_info_obj;
    
    # Others #
    sub categorize_by_MTs;
    sub build_MSAs;
    sub build_consensi;
    sub _build_con_muscle;
    sub _build_con_clustalw;
    sub _build_con_no_msa;
    sub _parse_read;
    sub print_MTs_fasta;
    sub print_MT_depth_output;
    sub print_sample_summary;
    sub print_mismatch_seqs;
    sub print_consensi;
    sub print_categorizable_reads;
    sub print_single_category_reads;
    sub qc;
    sub rm_tmp_files;
    
    
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
                               #$arg_href->{barcode},
                               $arg_href->{fwd_linker},
                               $arg_href->{fwd_primer},
                               $arg_href->{rev_primer},
                               $arg_href->{fwd_mer_len},
                               $arg_href->{fwd_max_shifts},
                               $arg_href->{min_con_depth},
                               $arg_href->{diginorm_max},
                               $arg_href->{output_dir},
            ) {
            MyX::Generic::Undef::Param->throw(
                                              error => 'Undefined parameter value',
                                              usage => $NEW_USAGE,
                                              );
        }
        
        # Bless a scalar to instantiate an object
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        # Initialize the objects attributes        
        $id_of{ident $new_obj} = $arg_href->{id};
        $name_of{ident $new_obj} = $arg_href->{name};
        $barcode_of{ident $new_obj} = $arg_href->{barcode} if (defined $arg_href->{barcode});
        $fwd_linker_of{ident $new_obj} = $arg_href->{fwd_linker};
        $fwd_primer_of{ident $new_obj} = $arg_href->{fwd_primer};
        $rev_primer_of{ident $new_obj} = $arg_href->{rev_primer};
        $fwd_mer_len_of{ident $new_obj} = $arg_href->{fwd_mer_len};
        $fwd_max_shifts_of{ident $new_obj} = $arg_href->{fwd_max_shifts};
        $min_con_depth_of{ident $new_obj} = $arg_href->{min_con_depth};
        $output_dir_of{ident $new_obj} = $arg_href->{output_dir};
        
        # there are some important operations in set_diginorm_max
        $new_obj->set_diginorm_max($arg_href->{diginorm_max});
        
        my @mismatched_fastq_seqs = ();
        $mismatches_of{ident $new_obj} = \@mismatched_fastq_seqs;
        my %MT_objs = ();  # This is a hash with KEY => molecule tag  VAL => molecule tag object
        $MT_objs_href_of{ident $new_obj} = \%MT_objs;
        
        # Set the MTToolbox::SampleSummary information object
        my $summ_info_obj = MTToolbox::SampleSummary->new();
        $summ_info_obj->set_sample_name($arg_href->{name});
        $summ_info_obj->set_sample_id($arg_href->{id});
        $new_obj->_set_summ_info_obj($summ_info_obj);
        
        
        return $new_obj;
    }
    
    
    ###########
    # Setters #
    ###########
    sub set_id {
        my ($self, $id) = @_;
        $id_of{ident $self} = $id;
        return 1;
    }
    
    sub set_name {
        my ($self, $n) = @_;
        $name_of{ident $self} = $n;
        return 1;
    }
    
    sub set_barcode {
        my ($self, $b) = @_;
        $barcode_of{ident $self} = $b;
        return 1;
    }
    
    sub set_fwd_linker {
        my ($self, $fwd_linker) = @_;
        $fwd_linker_of{ident $self} = $fwd_linker;
        return 1;
    }
    
    sub set_fwd_primer {
        my ($self, $fwd_primer) = @_;
        $fwd_primer_of{ident $self} = $fwd_primer;
        return 1;
    }
    
    sub set_rev_primer {
        my ($self, $rev_primer) = @_;
        $rev_primer_of{ident $self} = $rev_primer;
        return 1;
    }
    
    sub set_fwd_mer_len {
        my ($self, $fwd_mer_len) = @_;
        $fwd_mer_len_of{ident $self} = $fwd_mer_len;
        return 1;
    }
    
    sub set_fwd_max_shifts {
        my ($self, $fwd_max_shifts) = @_;
        $fwd_max_shifts_of{ident $self} = $fwd_max_shifts;
        return 1;
    }
    
    sub set_min_con_depth {
        my ($self, $min_con_depth) = @_;
        $min_con_depth_of{ident $self} = $min_con_depth;
        return 1;
    }
    
    sub set_diginorm_max {
        my ($self, $diginorm_max) = @_;
        if ( $diginorm_max =~ m/NA/i ) { $diginorm_max = $DIGI_MAX; }
        $diginorm_max_of{ident $self} = $diginorm_max;
        return 1;
    }
    
    sub set_output_dir {
        my ($self, $output_dir) = @_;
        $output_dir_of{ident $self} = $output_dir;
        return 1;
    }
    
    sub _set_summ_info_obj {
        my ($self, $summ_info_obj) = @_;
        $summ_info_obj_of{ident $self} = $summ_info_obj;
        return 1;
    }

    
    ###########
    # Getters #
    ###########
    sub get_id {
        my ($self) = @_;
        return $id_of{ident $self};
    }
    
    sub get_name {
        my ($self) = @_;
        return $name_of{ident $self};
    }
    
    sub get_barcode {
        my ($self) = @_;
        return $barcode_of{ident $self};
    }
    
    sub get_fwd_linker {
        my ($self) = @_;
        
        # return an empty string if undefined
        if ( ! defined $fwd_linker_of{ident $self} ) {
            return "";
        }
        
        return $fwd_linker_of{ident $self};
    }
    
    sub get_fwd_primer {
        my ($self) = @_;
        return $fwd_primer_of{ident $self};
    }
    
    sub get_rev_primer {
        my ($self) = @_;
        return $rev_primer_of{ident $self};
    }
    
    sub get_fwd_mer_len {
        my ($self) = @_;
        return $fwd_mer_len_of{ident $self};
    }
    
    sub get_fwd_max_shifts {
        my ($self) = @_;
        return $fwd_max_shifts_of{ident $self};
    }
    
    sub get_min_con_depth {
        my ($self) = @_;
        return $min_con_depth_of{ident $self};
    }
    
    sub get_diginorm_max {
        my ($self) = @_;
        return $diginorm_max_of{ident $self};
    }
    
    sub get_mismatch_seqs_aref {
        my ($self) = @_;
        return $mismatches_of{ident $self};
    }
    
    sub get_MT_objs_href {
        my ($self) = @_;
        return $MT_objs_href_of{ident $self};
    }
    
    sub get_output_dir {
        my ($self) = @_;
        return $output_dir_of{ident $self};
    }
    
    sub get_summ_info_obj {
        my ($self) = @_;
        return $summ_info_obj_of{ident $self};
    }
    
    
    
    
    ##########
    # Others #
    ##########
    sub categorize_by_MTs {
        my ($self) = @_;
        
        # Implement this in child class
    
        return -1;
    }
    
    sub build_MSAs {
        my ($self, $algo) = @_;
        
        if ( ! defined $algo ) { $algo = "Muscle"; }
        
        # Create a clustalw factory for bioPerl Clustalw
        #my @params = ('ktuple' => 2, 'matrix' => 'BLOSUM', 'quiet' => 1);
        #my $factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);
        
        # Directories for command line clustalw
        my $fasta_dir = $self->get_output_dir() . "/fasta_files";
        my $aln_dir = $self->get_output_dir() . "/aln_files";
        if ( -d $fasta_dir ) {
            `rm -rf $fasta_dir`; # A system command to erase already existing dirs
        }
        mkdir $fasta_dir;

        if ( -d $aln_dir ) {
            `rm -rf $aln_dir`; # A system command to erase already existing dirs
        }
        mkdir $aln_dir;
        
        my $temp_MT_objs_href = $self->get_MT_objs_href();
        foreach my $mt ( keys %{$temp_MT_objs_href} ) {
            if ( $temp_MT_objs_href->{$mt}->get_seq_count() >=
                $self->get_min_con_depth() ) {
                
                if ( $algo =~ m/Clustalw/ ) {
                    $temp_MT_objs_href->{$mt}->run_command_line_clustalw(
                        $fasta_dir, $aln_dir
                    );
                }
                else { # DEFAULT: Muscle
                    $temp_MT_objs_href->{$mt}->run_command_line_muscle(
                        $fasta_dir, $aln_dir
                    );
                }
            }
        }
        
        return 1;
    }
    
    sub build_consensi {
        my ($self, $method) = @_;
        
        if ( ! defined $method ) {
            # DEFAULT to clustalw if no method is given
            $self->_build_con_muscle();
            return 1;
        }
        
        if ( $method =~ /^NoMSA/i ) {
            $self->_build_con_no_msa();
        }
        elsif ( $method =~ m/^Muscle/i ) {
            $self->_build_con_muscle();
        }
        elsif ( $method =~ m/^Clustalw/i ) {
            $self->_build_con_clustalw();
        }
        else {
            # throw some error
            croak "Unrecognized con_algo ($method) in build_consensi()";
        }
        
        return 1;
    }
    
    sub _build_con_muscle {
        my ($self) = @_;
        
        my $aln_dir = $self->get_output_dir() . "/aln_files";
        
        my $temp_MT_objs_href = $MT_objs_href_of{ident $self};
        foreach my $mt ( keys %{$temp_MT_objs_href} ) {
            if ( $temp_MT_objs_href->{$mt}->get_seq_count() >=
                $self->get_min_con_depth() ) {
                $temp_MT_objs_href->{$mt}->build_con_from_fasta_file($aln_dir);
            }
        }
        
        return 1;
    }
    
    sub _build_con_clustalw {
        my ($self) = @_;
        
        my $aln_dir = $self->get_output_dir() . "/aln_files";
        
        my $temp_MT_objs_href = $MT_objs_href_of{ident $self};
        foreach my $mt ( keys %{$temp_MT_objs_href} ) {
            if ( $temp_MT_objs_href->{$mt}->get_seq_count() >=
                $self->get_min_con_depth() ) {
                $temp_MT_objs_href->{$mt}->build_con_from_clustalw_file($aln_dir);
            }
        }
        
        return 1;
    }
    
    sub _build_con_no_msa {
        my ($self) = @_;
        
        my $temp_MT_objs_href = $MT_objs_href_of{ident $self};
        foreach my $mt ( keys %{$temp_MT_objs_href} ) {
            if ( $temp_MT_objs_href->{$mt}->get_seq_count() >=
                $self->get_min_con_depth() ) {
                $temp_MT_objs_href->{$mt}->build_con();
            }
        }
        
        return 1;
    }
    
    sub _parse_read {
        my ($self) = @_;
        
        # Implement this in child class
        
        return -1;
    }
    
    sub print_MTs_fasta {
        my ($self) = @_;
        
        # Create/set the output dir
        my $out_dir = $self->get_output_dir() . "/fasta_files";
        if ( ! -d $out_dir ) { mkdir $out_dir or croak "Cannot mkdir: $out_dir\nERROR: $!\n"; }
        
        my $temp_MT_objs_href = $MT_objs_href_of{ident $self};
        foreach my $mt ( keys %{$temp_MT_objs_href} ) {
            my $out_file = $out_dir . "/" . $mt . ".fasta";
            $temp_MT_objs_href->{$mt}->print_clustalw_fasta($out_file);
        }
        
        return 1;
    }
    
    sub print_MT_depth_output {
        my ($self) = @_;
        my $dist_file = $self->get_output_dir() . "/seqs_per_mt_dist.txt";
        my $hist_file = $self->get_output_dir() . "/seqs_per_mt_hist.txt";
        my $hist_graph = $self->get_output_dir() . "/seqs_per_mt_hist.ps";
        
        my $depth_dist = MTToolbox::MTDepthDist->new({dist_file => $dist_file});
        $depth_dist->build_dist($MT_objs_href_of{ident $self}, $self->get_id());
        $depth_dist->build_hist();
        $depth_dist->print_hist($hist_file);
        $depth_dist->print_hist_graph($hist_graph,
                                      $self->get_id() . " MT Depth"
                                      );
        
        return 1;
    }
    
    sub print_sample_summary {
        my ($self) = @_;
        my $out_file = $self->get_output_dir() . "/summary.txt";
        
        open (SUM, ">$out_file") or croak "Cannot open $out_file\nERROR: $!\n";
        
        print SUM $summ_info_obj_of{ident $self}->to_string();
        
        close (SUM);
        
        return 1;
    }
    
    sub print_mismatch_seqs {
        my ($self) = @_;
        my $out_file = $self->get_output_dir() . "/mismatched.fastq";
        open (UN, ">$out_file") or croak "Cannot open $out_file\nERROR: $!\n";
        
        foreach ( @{$mismatches_of{ident $self}} ) {
            print UN $_->get_fastq_str();
        }
        
        close(UN);
        
        return 1;
    }
    
    sub print_consensi() {
        my ($self) = @_;
        
        my $out_file = $self->get_output_dir() . "/consensusSeqs.fastq";
        
        open (OUT, ">$out_file") or croak "Cannot open $out_file\nERROR: $!\n";
        
        # foreach molecule tag object that has a consensus sequence
        my $temp_MT_objs_href = $self->get_MT_objs_href();
        my $temp_id = $self->get_id();
        my $mt_obj;
        foreach my $mt ( keys %{$temp_MT_objs_href} ) {
            $mt_obj = $temp_MT_objs_href->{$mt};
            if ($mt_obj->get_seq_count() >= $self->get_min_con_depth() ) {
                # some MTs don't have a consensus sequence defined because
                # it was not able to be made (e.g. because all the reads were
                # different lengths).
                if ( $mt_obj->is_con_defined() ) {
                    print OUT $mt_obj->get_con_fastq_str($temp_id);
                }
            }
        }
        
        close(OUT);
        
        return 1;
    }
    
    sub print_categorizable_reads() {
        my ($self) = @_;
        
        my $out_file = $self->get_output_dir() . "/categorizable_reads.fastq";
        
        open (OUT, ">$out_file") or croak "Cannot open $out_file\nERROR: $!\n";
        
        my $temp_MT_objs_href = $self->get_MT_objs_href();
        foreach my $mt ( keys %{$temp_MT_objs_href} ) {
            print OUT $temp_MT_objs_href->{$mt}->get_unaligned_seqs_fastq_str();
        }
                
        close(OUT);
        
        return 1;
    }
    
    sub print_single_read_categories() {
        my ($self) = @_;
        
        # This method also sets the SRC_count value in the MTToolbox::SampleSummary obj
        
        my $out_file = $self->get_output_dir() . "/single_read_categories.fastq";
        my $count = 0;
        
        open (OUT, ">$out_file") or croak "Cannot open $out_file\nERROR: $!\n";
        
        my $temp_MT_objs_href = $self->get_MT_objs_href();
        foreach my $mt ( keys %{$temp_MT_objs_href} ) {
            if ( $temp_MT_objs_href->{$mt}->get_seq_count() == 1 ) {
                print OUT $temp_MT_objs_href->{$mt}->get_unaligned_seqs_fastq_str();
                $count++;
            }
        }
        
        close(OUT);
        
        my $sample_summ = $self->get_summ_info_obj();
        $sample_summ->set_SRC_count($count);
        
        return 1;
    }
    
    sub qc {
        my ($self, $qc_params_href) = @_;
        
        my $output_dir = $self->get_output_dir();
        
        # check if consensusSeqs.fastq file has been printed and print if needed
        if ( ! -s "$output_dir/consensusSeqs.fastq" ) {
            $self->print_consensi();
        }
        
        # filtering for consensus sequences
        my $con_filter = BioUtils::QC::FastqFilter->new(
            {
                fastq_file => "$output_dir/consensusSeqs.fastq",
                output_dir => $output_dir,
                trim_to_base => $qc_params_href->{con_trim_to_base},
                min_len => $qc_params_href->{con_min_len},
                min_avg_qual => $qc_params_href->{con_min_avg_qual},
                min_c_score => $qc_params_href->{con_min_c_score},
                allow_gaps => $qc_params_href->{con_allow_gaps},
                allow_ambig_bases => $qc_params_href->{con_allow_ambig},
                verbose => 1,
            });
        $con_filter->filter();
        
        # check if single_read_categorie.fastq file has been printed and print
        #   if needed
        if ( ! -s "$output_dir/single_read_categories.fastq" ) {
            $self->print_single_read_categories();
        }
        
        # filtering for single category reads (SRCs)
        my $SRC_filter = BioUtils::QC::FastqFilter->new(
            {
                fastq_file => "$output_dir/single_read_categories.fastq",
                output_dir => $output_dir,
                trim_to_base => $qc_params_href->{SRC_trim_to_base},
                min_len => $qc_params_href->{SRC_min_len},
                min_avg_qual => $qc_params_href->{SRC_min_avg_qual},
                allow_gaps => $qc_params_href->{SRC_allow_gaps},
                allow_ambig_bases => $qc_params_href->{SRC_allow_ambig},
                verbose => 1,
            }
        );
        $SRC_filter->filter();
        
        return 1;
    }
    
    sub rm_tmp_files {
        my ($self) = @_;
        
        my $output_dir = $self->get_output_dir();
        system("rm -rf $output_dir/aln_files");
        system("rm -rf $output_dir/fasta_files");
        
        return 1;
    }
}

1;
__END__

#######
# POD #
#######
=head1 NAME

MTToolbox::Sample - An object representing a biological sample

=head1 VERSION

This documentation refers to MTToolbox::Sample version 4.1.2.

=head1 Included Modules

    Bio::Tools::Run::Alignment::Clustalw
    Class::Std::Utils
    List::MoreUtils qw(any)
    Readonly
    Carp qw(carp croak)
    Chart::Gnuplot
    version
    MTToolbox::Mismatch
    MTToolbox::MoleculeTagCategory
    MTToolbox::SampleSummary
    MTToolbox::MTDepthDist
    BioUtils::FastqSeq 1.0.0
    BioUtils::FastqIO 1.0.0
    BioUtils::QC::FastqFilter 1.0.7
    MyX::Generic

=head1 Inherit

NA

=head1 SYNOPSIS

    use MTToolbox::Sample;
    my $my_sample = MTToolbox::Sample->new({id => $id,
                                 name => $name,
                                 barcode => $barcode,
                                 fwd_linker => $fwd_linker,
                                 fwd_primer => $fwd_primer,
                                 rev_primer => $rev_primer,
                                 fwd_mer_len => $fwd_mer_len,
                                 fwd_max_shifts => $fwd_max_shifts,
                                 min_con_depth => $min_con_depth,
                                 diginorm_max => $diginorm_max,
                                 output_dir => $output_dir,
                                 });
    
    my $id = $my_sample->get_id();
    my $name = $my_sample->get_name();
    my $barcode = $my_sample->get_barcode();
    my $fwd_linker = $my_sample->get_fwd_linker();
    my $fwd_primer = $my_sample->get_fwd_primer();
    my $rev_primer = $my_sample->get_rev_primer();
    my $fwd_mer_len = $my_sample->get_fwd_mer_len();
    my $fwd_max_shifts = $my_sample->get_fwd_max_shifts();
    my $min_con_depth = $my_sample->get_min_con_depth();
    my $diginorm_max = $my_sample->get_diginorm_max();
    my $output_dir = $my_sample->get_output_dir();
    my $summ_info_obj = $my_sample->get_summ_info_obj();
    
    $my_sample->set_id($id);
    my $sample->set_name($name);
    $my_sample->set_fwd_primer($b);
    $my_sample->set_fwd_linker($s);
    $my_sample->set_fwd_primer($f);
    $my_sample->set_rev_primer($r);
    $my_sample->set_fwd_mer_len($l);
    $my_sample->set_fwd_max_shifts($m);
    $my_sample->set_min_con_depth($c);
    $my_sample->set_diginorm_max($d);
    $my_sample->set_output_dir($o);
    $my_sample->_set_summ_info_obj($summ_onfo_obj);  # Be careful with this one!
    
    # These are the three major operations in a sample object
    $my_sample->categorize_by_MTs();  # implemented by child classes
    $my_sample->build_MSAs("Muscle");
    $my_sample->build_consensi("Muscle");
    
    # After running categorize_by_MTs() these functions can be run for summary information
    $my_sample->print_MT_depth_output();
    $my_sample->print_sample_summary();
    $my_sample->print_mismatch_seqs();
    
    # To print each sequence seperated by the molecule tag
    $my_sample->print_MTs_fasta();
    
    # To print the final consensus sequences (may be overridden by child classes)
    $my_sample->print_consensi();
    
    # To print the categorizable_reads (may be overridden by child classes)
    $my_sample->print_categorizable_reads();
    
    # To print the single read categories (may be overridden by child classes)
    $my_sample->print_single_read_categories();
    
    # To filter using BioUtils::QC::FastqFilter
    $my_sample->qc($qc_params_href);
    

=head1 DESCRIPTION

MTToolbox::Sample is an object that stores sequence information that comes from a distinct
biological sample.  During sequencing there are potential technical errors
including PCR errors and misread bases.  Using a method developed here at UNC, we
use this object to minimize these errors.  The first step to setting up a sample
is to read in the sequences and split them by their molecule tags.  After the
sequences are split, a multiple sequence alignment (MSA) can be made by the
group of sequence that represent a molecule tag.  With a MSA a consensus
sequence can be made.  This consensus sequence is an accurate representation of
the original DNA molecule that was sequenced.

=head1 METHODS

=over

    new
    set_id
    set_name
    set_barcode
    set_fwd_linker
    set_fwd_primer
    set_rev_primer
    set_fwd_mer_len
    set_fwd_max_shifts
    set_min_con_depth
    set_diginorm_max
    set_output_dir
    _set_summ_info_obj
    get_id
    get_name
    get_barcode
    get_fwd_linker
    get_fwd_primer
    get_rev_primer
    get_fwd_mer_len
    get_fwd_max_shifts
    get_min_con_depth
    get_diginorm_max
    getOuputDir
    get_mismatch_seqs_aref
    get_MT_objs_href
    get_summ_info_obj
    categorize_by_MTs
    build_MSAs
    build_consensi
    _build_con_muscle
    _build_con_clustalw
    _build_con_no_msa
    _parse_read
    print_MTs_fasta
    print_MT_depth_output
    print_sample_summary
    print_mismatch_seqs
    print_consensi
    print_categorizable_reads
    print_single_read_categories
    qc
    rm_tmp_files
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: MTToolbox::Sample->new({id => $id,
                        barcode => $barcode,
                        fwd_linker => $fwd_linker,
                        fwd_primer => $fwd_primer,
                        rev_primer => $rev_primer,
                        fwd_mer_len => $fwd_mer_len,
                        fwd_max_shifts => $fwd_max_shifts,
                        min_con_depth => $min_con_depth,
                        diginorm_max => $diginorm_max,
                        output_dir => $output_dir
                        });
    Function: Creates a new MTToolbox::Sample object
    Returns: MTToolbox::Sample
    Args: -id             => a string representing the sample id
          -barcode        => a string representing the sample barcode
          -fwd_linker     => a string representing the fwd_linker sequence
          -fwd_primer     => a string representing the forward primer sequence
          -rev_primer     => a string representing the reverse primer sequence
          -fwd_mer_len    => an int representing the length of the fwd random mer
          -fwd_max_shifts => an int representing the max fwd frameshifts
          -min_con_depth  => an int representing the min number of reads needed
                             to make a consensus sequence
          -diginorm_max   => an int representing the max number of reads to
                             consider when making a consensus. (NA uses all seqs)
          -output_dir     => a string representing the path of the output
                             directory

=head2 set_id

    Title: set_id
    Usage: $my_sample->set_id("P1");
    Function: Sets the sample ID
    Returns: 1 on successful completion
    Args: -id => a unique string representing the sample ID
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 set_name

    Title: set_name
    Usage: $my_sample->set_name("sample1");
    Function: Sets the sample name
    Returns: 1 on successful completion
    Args: -name => a string representing the sample name
    Throws: NA
    Comments: NA
    See Also: NA
 
=head2 set_barcode

    Title: set_barcode
    Usage: $my_sample->set_barcode("ATCGG");
    Function: Sets the barcode value
    Returns: 1 on successful completion
    Args: -barcode => a string representing the sample barcode sequence
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 set_fwd_linker
    
    Title: set_fwd_linker
    Usage: $my_sample->set_fwd_linker("CAGT");
    Function: Sets the forward linker value
    Returns: 1 on successful completion
    Args: -fwd_linker => a string representing the sample forward linker sequence
    Throws: NA
    Comments: NA
    See Also: NA

=head2 set_fwd_primer

    Title: set_fwd_primer
    Usage: $my_sample->set_fwd_primer("ATGTGAGATATGGC");
    Function: Sets the fwdPrimer value
    Returns: 1 on successful completion
    Args: -fwd_primer => a string representing the forward primer sequence
    Throws: NA
    Comments: Can contain abiguous bases using a regex format
              (eg. GTGCCAGC[AC]GCCGCGGTAA)
    See Also: NA
    
=head2 set_rev_primer

    Title: set_rev_primer
    Usage: $my_sample->set_rev_primer("GTAGATGTATGAG");
    Function: Sets the revPrimer value
    Returns: 1 on successful completion
    Args: -rev_primer => a string representing the reverse primer sequence
    Throws: NA
    Comments: Can contain abiguous bases using a regex format
              (eg. GTGCCAGC[AC]GCCGCGGTAA)
    See Also: NA
    
=head2 set_fwd_mer_len

    Title: set_fwd_mer_len
    Usage: $my_sample->set_fwd_mer_len(9);
    Function: Sets the length of the fwd random mer
    Returns: 1 on successful completion
    Args: -fwd_mer_len => an int representing the lenght of the fwd random mer
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 set_fwd_max_shifts

    Title: set_fwd_max_shifts
    Usage: $my_sample->set_fwd_max_shifts(5);
    Function: Sets the maximum number of frame shifts in the fwd tag
    Returns: 1 on successful completion
    Args: -fwd_max_shifts => an int representing the max number of fwd frame
                             shifts
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 set_min_con_depth

    Title: set_min_con_depth
    Usage: $my_sample->set_min_con_depth(2);
    Function: Sets min reads required to build a consensus
    Returns: 1 on successful completion
    Args: -min_con_depth => an int representing the minimum number of reads
                            in order to build a consensus sequences
    Throws: NA
    Comments: Default: 2
    See Also: NA
    
=head2 set_diginorm_max

    Title: set_diginorm_max
    Usage: $my_sample->set_diginorm_max(2);
    Function: Sets the max reads to consider for building consensus seqs
    Returns: 1 on successful completion
    Args: -diginorm_max => an int representing the number of reads to consider
                           when building a consensus sequence
    Throws: NA
    Comments: The string 'NA' signifies the program to use ALL the reads when
              building the consensus sequence.
              Default: NA
    See Also: NA
    
=head2 set_output_dir

    Title: set_output_dir
    Usage: $my_sample->set_output_dir("/home/scott/my_output/");
    Function: Sets the outputDir value
    Returns: 1 on successful completion
    Args: -output_dir => a string representing the path to a directory to print
                         the output files
    Throws: NA
    Comments: For safety always use a full path
    See Also: NA

=head2 _set_summ_info_obj

    Title: _set_summ_info_obj
    Usage: $my_sample->_set_summ_info_obj($my_summ_info_obj);
    Function: Stores a summary information object
    Returns: 1 on successful completion
    Args: -my_summ_info_obj => a MTToolbox::SampleSummary object
    Throws: NA
    Comments: Be very careful setting this.  This method is used in the child
              classes when categorizing molecule tagged sequences.  It
              automatically calculates the correct summary information for that
              sample.  The MTToolbox::SampleSummary object that it generates and
              stores can be replaced with a manually created one.  This is
              generally not recommended!
    See Also: MTToolbox::SampleSummary

=head2 get_id

    Title: get_id
    Usage: my $id = $my_sample->get_id();
    Function: Gets the sample ID
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_name

    Title: get_name
    Usage: my $id = $my_sample->get_name();
    Function: Gets the sample name
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_barcode

    Title: get_barcode
    Usage: my $barcode = $my_sample->get_barcode();
    Function: Gets the sample barcode sequence
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_fwd_linker

    Title: get_fwd_linker
    Usage: my $fwd_linker = $my_sample->get_fwd_linker();
    Function: Gets the sample forward linker sequence
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_fwd_primer

    Title: get_fwd_primer
    Usage: my $fwd_primer = $my_sample->get_fwd_primer();
    Function: Gets the sample forward primer sequence
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_rev_primer

    Title: get_rev_primer
    Usage: my $rev_primer = $my_sample->get_rev_primer();
    Function: Gets the sample reverse primer sequence
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_fwd_mer_len

    Title: get_fwd_mer_len
    Usage: my $fwd_mer_len = $my_sample->get_fwd_mer_len();
    Function: Gets the length of the fwd random mer
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_fwd_max_shifts

    Title: get_fwd_max_shifts
    Usage: my $fwd_max_shifts = $my_sample->get_fwd_max_shifts();
    Function: Gets the maximum number of allowed fwd frameshifts in the tag
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_min_con_depth

    Title: get_min_con_depth
    Usage: my $min_con_depth = $my_sample->get_min_con_depth();
    Function: Gets the minimum number of reads required to make a consensus
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_diginorm_max

    Title: get_diginorm_max
    Usage: my $diginorm_max = $my_sample->get_diginorm_max();
    Function: Gets the maximum number of reads to used for building a consensus
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_output_dir

    Title: get_output_dir
    Usage: my $ouptu_dir = $my_sample->get_output_dir();
    Function: Gets the sample output directory path
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_mismatch_seqs_aref

    Title: get_mismatch_seqs_aref
    Usage: my $mismatch_seqs = $my_sample->get_mismatch_seqs_aref();
    Function: Gets sequences that don't match the expected pattern
    Returns: Array reference
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_MT_objs_href

    Title: get_MT_objs_href
    Usage: my $mt_objs_href = $my_sample->get_MT_objs_href();
    Function: Gets the molecule tag objects
    Returns: Hash reference with KEY => tag, VALUE => MT_object
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_summ_info_obj

    Title: get_summ_info_obj
    Usage: my $summ_info_obj = $my_sample->get_summ_info_obj();
    Function: Gets the summary info object
    Returns: MTToolbox::SampleSummary object
    Args: None
    Throws: NA
    Comments: It is not recommended to use this method to get the summary
              information.  It would be better to call print_sample_summary()
    See Also: print_sample_summary()

=head2 categorize_by_MTs

    Title: categorize_by_MTs
    Usage: $my_sample->categorize_by_MTs();
    Function: Categorizes the raw sequences by their molecule tag(s)
    Returns: -1
    Args: None
    Throws: NA
    Comments: This method is not implemented in this class.  Because there are
              several patterns that can be used to parse sequences and therefore
              build several different types of match objects, this method is
              implemented in each child class of MTToolbox::Sample.  For
              example, there are sequences that are single end (SE) reads which
              are parsed differently than sequences with paired end (PE) reads.
    See Also: MTToolbox::Match.pm,
              MTToolbox::Match::SEMatch.pm,
              MTToolbox::Match::PEMatch.pm,
              MTToolbox::Match::PEMatch::Overlap.pm,
              MTToolbox::Match::PEMatch::NoOverlap.pm

=head2 build_MSAs

    Title: build_MSAs
    Usage: $my_sample->build_MSAs($algo);
    Function: Foreach molecule tag object build a multiple sequence alignment
              using each of the sequences associated with that molecule tag
    Returns: 1 on successful completion
    Args: -algo => The MSA algorithm to use.  Presently MT-Toolbox supports
                   Muscle and Clustalw.  DEFAULT: Muscle
    Throws: NA
    Comments: NA
    See Also: NA

=head2 build_consensi

    Title: build_consensi
    Usage: $my_sample->build_consensi($method);
    Function: Foreach molecule tag category object using its MSA created by
              build_MSAs(), build a consensus sequence to represent the original
              DNA molecule
    Returns: 1 on successful completion
    Args: -method => Muscle | Clustalw | NoMSA.  DEFAULT: Muscle
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 _build_con_muscle

    Title: _build_con_muscle
    Usage: $my_sample->_build_con_muscle();
    Function: Builds a consensus sequence from a muscle alignment fasta output
              file for each molecule tag category.
    Returns: 1 on successful completion
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 _build_con_clustalw

    Title: _build_con_clustalw
    Usage: $my_sample->_build_con_clustalw();
    Function: Builds a consensus sequence from a clustalw alignment gde output
              file for each molecule tag category.
    Returns: 1 on successful completion
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 _build_con_no_msa

    Title: _build_con_no_msa
    Usage: $my_sample->_build_con_no_msa();
    Function: Builds a consensus sequence from a set of sequences without doing
              a MSA for each molecule tag category
    Returns: 1 on successful completion
    Args: NA
    Throws: NA
    Comments: NA
    See Also: NA

=head2 _parse_read

    Title: _parse_read
    Usage: $my_sample->_parse_read();
    Function: Parse raw reads looking for the expected pattern
    Returns: -1
    Args: None
    Throws: NA
    Comments: This method is not implemented in this class.  Because there are
              several patterns that can be used to parse sequences and therefore
              build several different types of match objects, this method is
              implemented in each child class of MTToolbox::Sample.  For
              example, there are sequences that are single end (SE) reads which
              are parsed differently than sequences with paired end (PE) reads.
              This is also a local method meaning it should not be called from
              outside this class.
    See Also: NA

=head2 print_MTs_fasta

    Title: print_MTs_fasta
    Usage: $my_sample->print_MTs_fasta();
    Function: Prints all molecule tag objects sequences in seperate files
    Returns: 1 on successful completion
    Args: None
    Throws: NA
    Comments: This method creates a directory called fasta_files in the
              outputDir specified in the constructor or via set_output_dir.
              A fasta file is created for each molecule tag and files ared named
              by their tag sequence.  Each of the fasta files contains one or
              more sequences representing all sequences that have that particular
              molecule tag.
    See Also: NA
    
=head2 print_MT_depth_output

    Title: print_MT_depth_output
    Usage: $my_sample->print_MT_depth_output();
    Function: Prints the distribution of sequences per molecule tag
    Returns: 1 on successful completion
    Args: None
    Throws: NA
    Comments: Creates a file called 'seqs_per_mt_dist.txt' which contains the
              distribution of sequences per molecule tag.  This file can be used
              as input to gnuplot to build a graphical representation of this
              distribution.
    See Also: NA

=head2 print_sample_summary

    Title: print_sample_summary
    Usage: $my_sample->print_sample_summary();
    Function: Calls the toString method of the sample summary info object and
              prints that output
    Returns: 1 on successful completion
    Args: None
    Throws: NA
    Comments: Creates a file called 'summary.txt' which contains the string
              returned by the MTToolbox::SampleSummary object method toString()
    See Also: MTToolbox::SampleSummary.pm

=head2 print_mismatch_seqs

    Title: print_mismatch_seqs
    Usage: $my_sample->print_mismatch_seqs();
    Function: Prints the raw sequences that don't matche the expected pattern
    Returns: 1 on successful completion
    Args: None
    Throws: NA
    Comments: Creates a file called 'mismatched.fastq' with the raw sequences
              that don't match the expect pattern (defined in the child class)
              in FASTQ format.
    See Also: NA

=head2 print_consensi

    Title: print_consensi
    Usage: $my_sample->print_consensi();
    Function: Foreach molecule tag category, print its consensus sequence in
              FASTQ format
    Returns: -1
    Args: None
    Throws: NA
    Comments: This method is implemented in the child classes of
              MTToolbox::Sample
    See Also: MTToolbox::Sample::SE::Sample,
              MTToolbox::Sample::PE::Sample::NoOverlap,
              MTToolbox::Sample::PE::Sample::Overlap
    
=head2 print_categorizable_reads

    Title: print_categorizable_reads
    Usage: $my_sample->print_categorizable_reads();
    Function: Print the categorizable reads in FASTQ format
    Returns: 1 on successful completion
    Args: None
    Throws: NA
    Comments: The file that the sequences are printed to is named
              $outputDir/categorizable_reads.fastq
    See Also: NA

=head2 print_single_read_categories

    Title: print_single_read_categories
    Usage: $my_sample->print_single_read_categories();
    Function: Print the categorizable reads that have only one read in that
              category in FASTQ format
    Returns: 1 on successful completion
    Args: None
    Throws: NA
    Comments: The files that the sequences are printed to is named
              $outputDir/single_read_categories.fastq.  Also sets the single
              category reads count value in the MTToolbox::SampleSummary object.
    See Also: NA

=head2 qc

    Title: qc
    Usage: $my_sample->qc($qc_params_href);
    Function: Run quality control on the consensus seqs and single read
              category seqs.
    Returns: 1 on successful completion
    Args: qc_params_href => a hash ref with qc params
    Throws: NA
    Comments: See the BioUtils::QC::FastqFilter documentation for information
              about the parameters that should/can be included in the
              qc_params_href
    See Also: BioUtils::QC::FastqFilter
    
=head2 rm_tmp_files

    Title: rm_tmp_files
    Usage: $my_sample->rm_tmp_files();
    Function: Remove the tmp files.  Specifically the fasta_files dir and
              aln_files dir.
    Returns: 1 on successful completion
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
