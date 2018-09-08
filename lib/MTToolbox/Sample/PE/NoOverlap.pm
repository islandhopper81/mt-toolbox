package MTToolbox::Sample::PE::NoOverlap;

use strict;
use warnings;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(carp croak);
use MyX::Generic 0.0.3;
use version; our $VERSION = qv('4.1.2');

use MTToolbox::Mismatch;
use MTToolbox::Match::PE::NoOverlap;
use MTToolbox::MoleculeTagCategory::PE::NoOverlap;
use MTToolbox::SampleSummary::PE::NoOverlap;
use MTToolbox::Sample;
use base qw(MTToolbox::Sample);  # inherits from MTToolbox::Sample

{
    Readonly my $NEW_USAGE => q{ new( {id => ,
                                       barcode => ,
                                       fwd_linker => ,
                                       rev_linker => ,
                                       fwd_primer => ,
                                       rev_primer => ,
                                       fwd_mer_len => ,
                                       rev_mer_len => ,
                                       fwd_max_shifts => ,
                                       rev_max_shifts => ,
                                       min_con_depth => ,
                                       diginorm_max => ,
                                       output_dir => ,
                                       fwd_reads_file => ,
                                       rev_reads_file => ,
                                       } ) };

    # Attributes #
    my %fwd_reads_file_of;
    my %rev_reads_file_of;
    my %rev_mer_len_of;
    my %rev_max_shifts_of;
    my %rev_linker_of;
    
    # Setters #
    sub set_fwd_reads_file;
    sub set_rev_reads_file;
    sub set_rev_mer_len;
    sub set_rev_max_shifts;
    sub set_rev_linker;
    
    # Getters #
    sub get_fwd_reads_file;
    sub get_rev_reads_file;
    sub get_rev_mer_len;
    sub get_rev_max_shifts;
    sub get_rev_linker;
    
    # Others #
    sub categorize_by_MTs;
    sub print_consensi;
    sub print_categorizable_reads;
    sub print_single_read_categories;
    sub qc;
    sub _parse_read;
    sub _parse_fwd_read;
    sub _parse_rev_read;
    
    
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
                               $arg_href->{fwd_reads_file},
                               $arg_href->{rev_reads_file},
                               $arg_href->{rev_linker},
                               $arg_href->{rev_mer_len},
                               $arg_href->{rev_max_shifts},
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
        $fwd_reads_file_of{ident $new_obj} = $arg_href->{fwd_reads_file};
        $rev_reads_file_of{ident $new_obj} = $arg_href->{rev_reads_file};
        $rev_mer_len_of{ident $new_obj} = $arg_href->{rev_mer_len};
        $rev_max_shifts_of{ident $new_obj} = $arg_href->{rev_max_shifts};
        $rev_linker_of{ident $new_obj} = $arg_href->{rev_linker};
        
        # Set the MTToolbox::SampleSummary object
        my $summ_info_obj = MTToolbox::SampleSummary::PE::NoOverlap->new();
        $summ_info_obj->set_sample_name($new_obj->get_name());
        $summ_info_obj->set_sample_id($arg_href->{id});
        $new_obj->_set_summ_info_obj($summ_info_obj);
        
        return $new_obj;
    }
    
    
    ###########
    # Setters #
    ###########
    sub set_fwd_reads_file {
        my ($self, $file) = @_;
        $fwd_reads_file_of{ident $self} = $file;
        return 1;
    }
    
    sub set_rev_reads_file {
        my ($self, $file) = @_;
        $rev_reads_file_of{ident $self} = $file;
        return 1;
    }
    
    sub set_rev_mer_len {
        my ($self, $rev_mer_len) = @_;
        $rev_mer_len_of{ident $self} = $rev_mer_len;
        return 1;
    }
    
    sub set_rev_max_shifts {
        my ($self, $rev_max_shifts) = @_;
        $rev_max_shifts_of{ident $self} = $rev_max_shifts;
        return 1;
    }
    
    sub set_rev_linker {
        my ($self, $rev_linker) = @_;
        $rev_linker_of{ident $self} = $rev_linker;
        return 1;
    }
    
    
    
    
    ###########
    # Getters #
    ###########
    sub get_fwd_reads_file {
        my ($self) = @_;
        return $fwd_reads_file_of{ident $self};
    }
    
    sub get_rev_reads_file {
        my ($self) = @_;
        return $rev_reads_file_of{ident $self};
    }
    
    sub get_rev_mer_len {
        my ($self) = @_;
        return $rev_mer_len_of{ident $self};
    }
    
    sub get_rev_max_shifts {
        my ($self) = @_;
        return $rev_max_shifts_of{ident $self};
    }
    
    sub get_rev_linker {
        my ($self) = @_;
        
        # return an empty string if undefined
        if ( ! defined $rev_linker_of{ident $self} ) {
            return "";
        }
        
        return $rev_linker_of{ident $self};
    }
    
    
    
    
    ##########
    # Others #
    ##########
    sub categorize_by_MTs {
        my ($self) = @_;
        
        ### Read the fwd and rev seqs into a hash
        # KEY => seq header  VAL => seq object    
        my %seq_hash = ();
        
        # Create a FastqIO object for the forward reads
        my $fwd_seqs_in;
        eval {
            $fwd_seqs_in = BioUtils::FastqIO->new({
                                            stream_type => '<',
                                            file => $self->get_fwd_reads_file(),
                                        });
        };
        if ( my $e = MyX::Generic::Undef::Param->caught() ) {
            croak("$e->error()\n$e->trace()\n");
        }
        else { print $@; }
        
        # Reading in fwd seqs
        while ( my $seq = $fwd_seqs_in->get_next_seq() ) {
            $seq_hash{$seq->get_id()}{fwd} = $seq;
        }
        
        # Create a FastqIO object for the reverse reads
        my $rev_seqs_in;
        eval {
            $rev_seqs_in = BioUtils::FastqIO->new({
                                            stream_type => '<',
                                            file => $self->get_rev_reads_file(),
                                        });
        };
        if ( my $e = MyX::Generic::Undef::Param->caught() ) {
            croak("$e->error()\n$e->trace()\n");
        }
        else { print $@; }
        
        # Reading in rev seqs
        while ( my $seq = $rev_seqs_in->get_next_seq() ) {
            $seq_hash{$seq->get_id()}{rev} = $seq;
        }
        
        ### Create either a MTToolbox::Match::PE::NoOverlap or Mismatch object
        ### for each seq and attempt to add MTToolbox::Match::PE::NoOverlap
        ### objects to a MoleculeTag object Also, store some summary info as I
        ### go through the sequences
        my $matched_count = 0;
        my $unmatched_count = 0;
        my $total_seqs = 0;
        my $mt_count = 0;
        my $temp_mismatch_seqs_aref = $self->get_mismatch_seqs_aref();
        my $temp_MT_objs_href = $self->get_MT_objs_href();
        
        foreach my $id ( keys %seq_hash ) {
            
            # Match the fwd and rev sequences for the given id
            if ( !defined $seq_hash{$id}{fwd} ) {
                carp "No matching fwd read for $id";
                next; # skip any sequences that don't have a matching fwd read
            }
            if ( !defined $seq_hash{$id}{rev} ) {
                carp "No matching rev read for $id";
                next; # skip any sequences that don't have a matching rev read
            }
            
            # header_prefix is a variable to be the prefix to the current sequence
            # that contains the information on the sample and sequence count
            my $header_prefix = $self->get_id() . "_" . $total_seqs;
            
            ### Try creating a MTToolbox::Match::PE::NoOverlap object with the
            # paired fwd and rev read
            # After this line $match is either -1 meaning it was not a match or
            # it is a MTToolbox::Match::PE::NoOverlap object.
            my $match = $self->_parse_read( $seq_hash{$id}{fwd}, $seq_hash{$id}{rev}, $header_prefix );
            
            ### Check if the MTToolbox::Match::PE::NoOverlap object was
            # successfully created.  If not create an Mismatch object
            if ( $match == -1 ) {
                $unmatched_count++;
                push @{$temp_mismatch_seqs_aref}, MTToolbox::Mismatch->new( {
                                                                  id => $header_prefix,
                                                                  desc => $id . "_fwd",
                                                                  seq => $seq_hash{$id}{fwd}->get_seq(),
                                                                  quals_str => $seq_hash{$id}{fwd}->get_quals_str()
                                                                 });
                push @{$temp_mismatch_seqs_aref}, MTToolbox::Mismatch->new( {
                                                                  id => $header_prefix,
                                                                  desc => $id . "_rev",
                                                                  seq => $seq_hash{$id}{rev}->get_seq(),
                                                                  quals_str => $seq_hash{$id}{rev}->get_quals_str()
                                                                  });
            }
            else {
                $matched_count++;
                
                # create the tag so I can look for it in the hash or create a new hash entry
                my $tag = $match->get_fwd_tag()->get_seq() . "-" . $match->get_rev_tag()->get_seq(); 
                
                # Create a new MoleculeTag object if needed
                if ( ! defined $temp_MT_objs_href->{$tag} ) {
                    $mt_count++;
                    
                    # create the new MT object
                    my $mt_obj = MTToolbox::MoleculeTagCategory::PE::NoOverlap->new({
                                    tag => $tag
                                });
                    
                    # save the new MT object in this sample
                    $temp_MT_objs_href->{$tag} = $mt_obj;  
                }
                
                # Add the MTToolbox::Match::PE::NoOverlap to a MoleculeTag
                # object if the number or reads is less than diginorm_max
                # (e.g. digital normalization)
                if ($temp_MT_objs_href->{$tag}->get_seq_count() <
                    $self->get_diginorm_max() )
                {
                    $temp_MT_objs_href->{$tag}->add_match($match);
                }
            }
            
            # Go to the next sequence and increant the total sequence counter
            # NOTE: I am only counting the sequences that have matching fwd and rev reads here.
            $total_seqs++;
        }
        
        ### Store the summary information
        my $summ_info_obj = $self->get_summ_info_obj();
        $summ_info_obj->set_file_format('fastq');
        $summ_info_obj->set_match_count($matched_count);
        $summ_info_obj->set_mismatch_count($unmatched_count);
        $summ_info_obj->set_seq_count($total_seqs);
        $summ_info_obj->set_MT_count($mt_count);
        
        return 1;
    }
    
    sub print_consensi() {
        my ($self) = @_;
        my $fwd_out_file = $self->get_output_dir() . "/consensusSeqs_fwd.fastq";
        my $rev_out_file = $self->get_output_dir() . "/consensusSeqs_rev.fastq";
        
        open (FWD, ">$fwd_out_file") or croak "Cannot open $fwd_out_file\nERROR: $!\n";
        open (REV, ">$rev_out_file") or croak "Cannot open $rev_out_file\nERROR: $!\n";
        
        # foreach molecule tag object that has a consensus sequence
        my $temp_MT_objs_href = $self->get_MT_objs_href();
        my $temp_id = $self->get_id();
        foreach my $mt ( keys %{$temp_MT_objs_href} ) {
            if ($temp_MT_objs_href->{$mt}->get_seq_count() > 1 ) {  # must have at least 2 seqs to have a defined consensus
                print FWD $temp_MT_objs_href->{$mt}->get_con_fastq_str($temp_id, 'fwd');
                print REV $temp_MT_objs_href->{$mt}->get_con_fastq_str($temp_id, 'rev');
            }
        }
        
        close(FWD);
        close(REV);
        
        return 1;
    }
    
    sub print_categorizable_reads() {
        my ($self) = @_;
        
        my $fwd_out_file = $self->get_output_dir() . "/categorizable_reads_fwd.fastq";
        my $rev_out_file = $self->get_output_dir() . "/categorizable_reads_rev.fastq";
        
        open (FWD, ">$fwd_out_file") or croak "Cannot open $fwd_out_file\nERROR: $!\n";
        open (REV, ">$rev_out_file") or croak "Cannot open $rev_out_file\nERROR: $!\n";
        
        my $temp_MT_objs_href = $self->get_MT_objs_href();
        foreach my $mt ( keys %{$temp_MT_objs_href} ) {
            print FWD $temp_MT_objs_href->{$mt}->get_unaligned_seqs_fastq_str('fwd');
            print REV $temp_MT_objs_href->{$mt}->get_unaligned_seqs_fastq_str('rev');
        }
                
        close(FWD);
        close(REV);
        
        return 1;
    }
    
    sub print_single_read_categories() {
        my ($self) = @_;
        
        my $fwd_out_file = $self->get_output_dir() . "/single_read_categories_fwd.fastq";
        my $rev_out_file = $self->get_output_dir() . "/single_read_categories_rev.fastq";
        my $count = 0;
        
        open (FWD, ">$fwd_out_file") or croak "Cannot open $fwd_out_file\nERROR: $!\n";
        open (REV, ">$rev_out_file") or croak "Cannot open $rev_out_file\nERROR: $!\n";
        
        my $temp_MT_objs_href = $self->get_MT_objs_href();
        foreach my $mt ( keys %{$temp_MT_objs_href} ) {
            if ( $temp_MT_objs_href->{$mt}->get_seq_count() == 1 ) {
                print FWD $temp_MT_objs_href->{$mt}->get_unaligned_seqs_fastq_str('fwd');
                print REV $temp_MT_objs_href->{$mt}->get_unaligned_seqs_fastq_str('rev');
                $count++;
            }
        }
        
        close(FWD);
        close(REV);
        
        my $sample_summ = $self->get_summ_info_obj();
        $sample_summ->set_SRC_count($count);
        
        return 1;
    }
    
    sub qc {
        my ($self, $qc_params_href) = @_;
        
        my $output_dir = $self->get_output_dir();
        
        # check if consensusSeqs.fastq file has been printed and print if needed
        if ( ! -s "$output_dir/consensusSeqs_fwd.fastq" or
             ! -s "$output_dir/consensusSeqs_rev.fastq" ) {
            $self->print_consensi();
        }
        
        # filtering for consensus sequences
        my $con_filter = BioUtils::QC::FastqFilter->new(
            {
                fastq_file => "$output_dir/consensusSeqs_fwd.fastq",
                fastq_file2 => "$output_dir/consensusSeqs_rev.fastq",
                output_dir => $output_dir,
                min_len => $qc_params_href->{con_min_len},
                min_avg_qual => $qc_params_href->{con_min_avg_qual},
                min_c_score => $qc_params_href->{con_min_c_score},
                allow_gaps => $qc_params_href->{con_allow_gaps},
                allow_ambig_bases => $qc_params_href->{con_allow_ambig},
                verbose => 1,
            });
        $con_filter->filter_pairs();
        
        # check if single_read_categorie.fastq file has been printed and print
        #   if needed
        if ( ! -s "$output_dir/single_read_categories_fwd.fastq" or
             ! -s "$output_dir/single_read_categories_rev.fastq" ) {
            $self->print_single_read_categories();
        }
        
        # filtering for single category reads (SRCs)
        my $SRC_filter = BioUtils::QC::FastqFilter->new(
            {
                fastq_file => "$output_dir/single_read_categories_fwd.fastq",
                fastq_file2 => "$output_dir/single_read_categories_rev.fastq",
                output_dir => $output_dir,
                min_len => $qc_params_href->{SRC_min_len},
                min_avg_qual => $qc_params_href->{SRC_min_avg_qual},
                allow_gaps => $qc_params_href->{SRC_allow_gaps},
                allow_ambig_bases => $qc_params_href->{SRC_allow_ambig},
                verbose => 1,
            }
        );
        $SRC_filter->filter_pairs();
        
        return 1;
    }
    
    sub _parse_read($$$) {
        my ($self, $fwd_read, $rev_read, $header_prefix) = @_;
        
        my ($is_fwd_match, $fwd_tag, $fwd_linker, $fwd_primer, $fwd_amplicon, $fwd_tail) =
            $self->_parse_fwd_read($fwd_read->get_seq(), $fwd_read->get_quals_str());
        my ($is_rev_match, $rev_tag, $rev_linker, $rev_primer, $rev_amplicon, $rev_tail) =
            $self->_parse_rev_read($rev_read->get_seq(), $fwd_read->get_quals_str());
        
        if ( $is_fwd_match and $is_rev_match ) {
            # NOTE: When passing in the Id for the sequence the fwd_read->id should be the same as the rev_read->id.
            my $desc = $fwd_tag->get_seq() . "-" . $rev_tag->get_seq() . " " .
                        $fwd_read->get_header();
            my $FastqPENoOverMatch = MTToolbox::Match::PE::NoOverlap->new({
                                        id => $header_prefix,
                                        desc => $desc,
                                        fwd_tag => $fwd_tag,
                                        fwd_linker => $fwd_linker,
                                        fwd_primer => $fwd_primer,
                                        fwd_amplicon => $fwd_amplicon,
                                        fwd_tail => $fwd_tail,
                                        rev_tag => $rev_tag,
                                        rev_linker => $rev_linker,
                                        rev_primer => $rev_primer,
                                        rev_amplicon => $rev_amplicon,
                                        rev_tail => $rev_tail
                                    });
            
            return $FastqPENoOverMatch;
        }
        
        return -1; # Indicates an unmatched sequence
    }
    
    sub _parse_fwd_read($$) {
        my ($self, $fwd_read, $quals_str) = @_;
        
        # NOTE: This regex looks for the reverse primer in the event that we actually get through it.
        #       If we do get through it we want to ignore what is after it because it is most likely
        #       the reverse random mer.
        my $rev_regex = "(\\S+" . $self->get_fwd_primer() . "\\S+)" . $self->get_rev_primer();
        if ( $fwd_read =~ m/$rev_regex/i ) {
            $fwd_read = $1;
            $quals_str = substr $quals_str, 0, length $1;
        }
        
        # The regex for the fwd read
        my $fwd_min = $self->get_fwd_mer_len();
        my $fwd_max = $fwd_min + $self->get_fwd_max_shifts();
        my $regex = "(^[ATCG]{" . $fwd_min . "," . $fwd_max . "})" .
                    $self->get_fwd_linker() .
                    "(" . $self->get_fwd_primer() . ")" .
                    "(\\S+)";
        
        # Return Values
        my $is_matched = 0; # This is a boolean to help me identify in the calling method if this worked.
        my ($fwd_tag, $fwd_linker, $fwd_primer, $fwd_amplicon, $fwd_tail);
        
        if ( $fwd_read =~ m/$regex/i ) {
            my $index = 0;  # this keeps track of the end of where I am so I can index the quals array correctly
            
            my $fwd_tag_quals_str = substr $quals_str, 0, length $1;
            $fwd_tag = BioUtils::FastqSeq->new( {
                            seq => $1,
                            quals_str => $fwd_tag_quals_str
                        } );
            $index = length $1;
            
            my $temp_fwd_linker = $self->get_fwd_linker();
            my $fwd_linker_quals_str = substr $quals_str, $index, length $temp_fwd_linker;
            $fwd_linker = BioUtils::FastqSeq->new( {
                                seq => $temp_fwd_linker,
                                quals_str => $fwd_linker_quals_str
                            } );
            $index += length($temp_fwd_linker);
            
            my $fwd_primer_quals_str = substr $quals_str, $index, length $2;
            $fwd_primer = BioUtils::FastqSeq->new( {
                                seq => $2,
                                quals_str => $fwd_primer_quals_str
                            } );
            $index += length $2;
            
            # find the length of the fwd tail and amplicon
            my $tail_len = $self->get_fwd_max_shifts() -
                            ( (length $1) - $self->get_fwd_mer_len() );
            my $amp_len = (length $3) - $tail_len;
            
            # amplicon
            my $fwd_amplicon_quals_str = substr $quals_str, $index, $amp_len;
            $fwd_amplicon = BioUtils::FastqSeq->new( {
                                seq => (substr $3, 0, $amp_len),
                                quals_str => $fwd_amplicon_quals_str
                            } );
            $index += $amp_len;
            
            # tail
            if ( $tail_len > 0 ) {
                my $fwd_tail_quals_str = substr $quals_str, $index, $tail_len;
                $fwd_tail = BioUtils::FastqSeq->new( {
                                seq => (substr $3, -$tail_len),
                                quals_str => $fwd_tail_quals_str
                            } );
            }
            else {
                $fwd_tail = BioUtils::FastqSeq->new( {seq => q{}, quals_str => q{}} );
            }
            
            # This is a boolean signifying that the sequence was a match to the given regex
            $is_matched = 1;
        }
        else {
            $is_matched = 0;
        }
        
        return ($is_matched, $fwd_tag, $fwd_linker, $fwd_primer, $fwd_amplicon, $fwd_tail);
    }
    
    sub _parse_rev_read($$) {    
        my ($self, $rev_read, $quals_str) = @_;
        
        my $rev_min = $self->get_rev_mer_len();
        my $rev_max = $rev_min + $self->get_rev_max_shifts();
        my $regex = "(^[ATCG]{" . $rev_min . "," . $rev_max . "})" .
                    $self->get_rev_linker() .
                    "(" . $self->get_rev_primer() . ")" .
                    "(\\S+)";
        
        # Return Values
        my $is_matched = 0; # This is a boolean to help me identify in the calling method if this worked.
        my ($rev_tag, $rev_linker, $rev_primer, $rev_amplicon, $rev_tail);
        
        if ( $rev_read =~ m/$regex/i ) {
            my $index = 0;  # this keeps track of the end of where I am so I can index the quals array correctly
            
            my $rev_tag_quals_str = substr $quals_str, $index, length $1;
            $rev_tag = BioUtils::FastqSeq->new( {
                            seq => $1,
                            quals_str => $rev_tag_quals_str
                        } );
            $index = length $1;
            
            my $temp_rev_linker = $self->get_rev_linker();
            my $rev_linker_quals_str = substr $quals_str, $index, length $temp_rev_linker;
            $rev_linker = BioUtils::FastqSeq->new( {
                                seq => $temp_rev_linker,
                                quals_str => $rev_linker_quals_str
                            } );
            $index += length($temp_rev_linker);
            
            my $rev_primer_quals_str = substr $quals_str, $index, length $2;
            $rev_primer = BioUtils::FastqSeq->new( {
                                seq => $2,
                                quals_str => $rev_primer_quals_str
                            } );
            $index += length $2;
            
            # find the length of the rev tail and amplicon
            my $tail_len = $self->get_rev_max_shifts() -
                            ( (length $1) - $self->get_rev_mer_len() );
            my $amp_len = (length $3) - $tail_len;
            
            # amplicon
            my $rev_amplicon_quals_str = substr $quals_str, $index, $amp_len;
            $rev_amplicon = BioUtils::FastqSeq->new( {
                                seq => (substr $3, 0, $amp_len),
                                quals_str => $rev_amplicon_quals_str
                            } );
            $index += $amp_len;
            
            # tail
            if ( $tail_len > 0 ) {
                my $rev_tail_quals_str = substr $quals_str, $index, $tail_len;
                $rev_tail = BioUtils::FastqSeq->new( {
                                seq => substr($3, -$tail_len),
                                quals_str => $rev_tail_quals_str
                            } );
            }
            else {
                $rev_tail = BioUtils::FastqSeq->new( {seq => q{}, quals_str => q{}} );
            }
            
            # This is a boolean signifying that the sequence was a match to the given regex
            $is_matched = 1;
        }
        else {
            $is_matched = 0;
        }
        
        return ($is_matched, $rev_tag, $rev_linker, $rev_primer, $rev_amplicon, $rev_tail);
    }
}


1;
__END__



#######
# POD #
#######
=head1 MTToolbox::Sample::PE::NoOverlap

MTToolbox::Sample::PE::NoOverlap - A class to store and handle samples with pair
                                   end sequences that do NOT overlap

=head1 VERSION

This documentation refers to MTToolbox::Sample::PE::NoOverlap version 4.1.2.

=head1 Included Modules

    Class::Std::Utils
    List::MoreUtils qw(any)
    Readonly
    Carp qw(carp croak)
    MyX::Generic 1.0.0
    version
	MTToolbox::Mismatch
    MTToolbox::Match::PE::NoOverlap
    MTToolbox::MoleculeTagCategory::PE::NoOverlap
    MTToolbox::SampleSummary::PE::NoOverlap
    MTToolbox::Sample

=head1 Inherit

    MTToolbox::Sample

=head1 SYNOPSIS

    # NOTE: MTToolbox::Sample::PE::NoOverlap inherits from MTToolbox::Sample.
    # See MTToolbox::Sample.pm for details on inherited methods.
    
    use MTToolbox::Sample::PE::NoOverlap;
    my $my_sample = MTToolbox::Sample::PE::NoOverlap->new({
                        id => $id,
                        barcode => $barcode,
                        fwd_reads_file => $fwd_reads_file,
                        rev_reads_file => $rev_reads_file,
                        fwd_linker => $fwd_linker,
                        rev_linker => $rev_linker,
                        fwd_primer => $fwd_primer,
                        rev_primer => $rev_primer,
                        fwd_mer_len => $fwd_mer_len,
                        rev_mer_len => $rev_mer_len,
                        fwd_max_shifts => $fwd_max_shifts,
                        rev_max_shifts => $rev_max_shifts,
                        min_con_depth => $min_con_depth,
                        diginorm_max => $diginorm_max,
                        output_dir => $output_dir,
                    });
    
    my $fwd_reads_file = $my_sample->get_fwd_reads_file();
    my $rev_reads_file = $my_sample->get_rev_reads_file();
    my $rev_mer_len = $my_sample->get_rev_mer_len();
    
    $my_sample->set_fwd_reads_file($fwd_file);
    $my_sample->set_rev_reads_file($rev_file);
    $my_sample->set_rev_mer_len($rev_mer_len);
    
    # The major operations implented in MTToolbox::Sample::PE::NoOverlap are the
    # categorizing and printing methods.  Other major operations are implemented
    # in the parent class MTToolbox::Sample
    $my_sample->categorize_by_MTs();
    $my_sample->print_consensi();

=head1 DESCRIPTION

MTToolbox::Sample::PE::NoOverlap is an object that stores sequence information
that comes from a distinct biological sample with paired end (PE) reads.  This
class inherits from MTToolbox::Sample.  See the documentation in
MTToolbox::Sample.pm for important information on inherited methods.  

=head1 METHODS

=over

    new
    set_fwd_reads_file
    set_rev_reads_file
    set_rev_mer_len
    set_rev_max_shifts
    set_rev_linker
    get_fwd_reads_file
    get_rev_reads_file
    get_rev_mer_len
    get_rev_max_shifts
    get_rev_linker
    categorize_by_MTs
    print_consensi
    print_categorizable_reads
    print_single_read_categories
    qc
    _parse_read
    _parse_fwd_read
    _parse_rev_read
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: MTToolbox::Sample::PE::NoOverlap->new({
                id => $id,
                barcode => $barcode,
                fwd_reads_file => $fwd_reads_file,
                rev_reads_file => $rev_readFile,
                fwd_linker => $fwd_linker,
                rev_linker => $rev_linker,
                fwd_primer => $fwd_primer,
                rev_primer => $rev_primer,
                fwd_mer_len => $fwd_mer_len,
                rev_mer_len => $rev_mer_len,
                fwd_max_shifts => $fwd_max_shifts,
                rev_max_shifts => $rev_max_shifts,
                min_con_depth => $min_con_depth,
                diginorm_max => $diginorm_max,
                output_dir => $output_dir,
            });
    Function: Creates a new MTToolbox::Sample::PE::NoOverlap object
    Returns: MTToolbox::Sample::PE::NoOverlap
    Args: -id             => a string representing the sample id
          -barcode        => a string representing the sample barcode
          -fwd_reads_file => a string representing the path to the raw FASTQ fwd
                             reads
          -rev_reads_file => a string representing the path to the raw FASTQ rev
                             reads
          -fwd_linker     => a string representing the forward linker sequence
          -rev_linker     => a string representing the reverse linker sequence
          -fwd_primer     => a string representing the forward primer sequence
          -rev_primer     => a string representing the reverse primer sequence
          -fwd_mer_len    => an int representing the length of the fwd random mer
          -rev_mer_len    => an int representing the length of the rev random mer
          -fwd_max_shifts => an int representing the max fwd frameshifts
          -rev_max_shifts => an int representing the max rev frameshifts
          -min_con_depth  => an int representing the min number of reads needed
                             to make a consensus sequence
          -diginorm_max   => an int representing the max number of reads to
                             consider when making a consensus. (NA uses all seqs)
          -output_dir     => a string representing the path of the output
                             directory

=head2 set_fwd_reads_file

    Title: set_fwd_reads_file
    Usage: $my_sample->set_fwd_reads_file("/home/Scott/data/SE_fwd_reads.fastq");
    Function: Sets the forward reads file path
    Returns: 1 on successful completion
    Args: -fwd_reads_file => a string representing the path to the raw FASTQ fwd reads
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 set_rev_reads_file

    Title: set_rev_reads_file
    Usage: $my_sample->set_rev_reads_file("/home/Scott/data/SE_rev_reads.fastq");
    Function: Sets the reverse reads file path
    Returns: 1 on successful completion
    Args: -rev_reads_file => a string representing the path to the raw FASTQ rev reads
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 set_rev_mer_len

    Title: set_rev_mer_len
    Usage: $my_sample->set_rev_mer_len(4);
    Function: Sets the lenght of the reverse random mer
    Returns: 1 on successful completion
    Args: -rev_mer_len => an int representing the length of the rev random mer
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 set_rev_max_shifts

    Title: set_rev_max_shifts
    Usage: $my_sample->set_rev_max_shifts(5);
    Function: Sets the maximum number of frameshifts in the rev tag
    Returns: 1 on successful completion
    Args: -rev_max_shifts => an int representing the max rev frameshifts
    Throws: NA
    Comments: NA
    See Also: NA

=head2 set_rev_linker

    Title: set_rev_linker
    Usage: $my_sample->set_rev_linker("GTAC");
    Function: Sets the reverse linker sequence
    Returns: 1 on successful completion
    Args: -rev_linker => a string representing the reverse linker sequence
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_fwd_reads_file

    Title: get_fwd_reads_file
    Usage: $my_sample->get_fwd_reads_file();
    Function: Gets the forward reads file path
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_rev_reads_file

    Title: get_rev_reads_file
    Usage: $my_sample->get_rev_reads_file();
    Function: Gets the reverse reads file path
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_rev_mer_len

    Title: get_rev_mer_len
    Usage: $my_sample->get_rev_mer_len();
    Function: Gets the length of the reverse random mer
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_rev_max_shifts

    Title: get_rev_max_shifts
    Usage: $my_sample->get_rev_max_shifts();
    Function: Gets the max number of rev tag frameshifts
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_rev_linker

    Title: get_rev_linker
    Usage: $my_sample->get_rev_linker();
    Function: Gets the sample reverse linker sequence
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 categorize_by_MTs

    Title: categorize_by_MTs
    Usage: $my_sample->categorize_by_MTs();
    Function: Categorizes the raw PE sequences by their molecule tag(s)
    Returns: 1 on successful completion
    Args: None
    Throws: NA
    Comments: This method uses the local methods _parse_read, _parse_fwd_read, and
              _parse_rev_read to attempt to match each raw read to the defined
              pattern.  There is a regex in _parse_fwd_read and _parse_rev_read that
              defines that pattern.  At some point this will hopefully turn into
              a user defined parameter.  If a pattern match is successful a
              MTToolbox::Match::PE::NoOverlap object is created and added to the
              MTObjects hash inherited from MTToolbox::Sample.  If there is no
              MoleculeTagCategory with that molecule tag in the hash, a new
              MTToolbox::MoleculeTagCategory::PE::NoOverlap object is created
              and added to the hash with KEY => molecule tag,
              VALUE => MTToolbox::MoleculeTagCategory::PE::NoOverlap object
    See Also: MTToolbox::Match::PE::NoOverlap, Match::PE::_parse_read,
              MTToolbox::MoleculeTagCategory::PE::NoOverlap,
              MTToolbox::Match::PE::_parse_fwd_read,
              MTToolbox::Match::PE::_parse_rev_read

=head2 print_consensi

    Title: print_consensi
    Usage: $my_sample->print_consensi();
    Function: Prints the consensus sequences in FASTQ format
    Returns: 1 on successful completion
    Args: None
    Throws: NA
    Comments: There are two output files that are created after running this
              method: $outputDir/consensusSeqs_fwd.fastq and
              $outputDir/consensusSeqs_rev.fastq.  Because the forward and
              reverse sequences don't overlap there are seperate consensus
              sequences for them.
    See Also: NA
    
=head2 print_categorizable_reads

    Title: print_categorizable_reads
    Usage: $my_sample->print_categorizable_reads();
    Function: Print the categorizable reads in FASTQ format
    Returns: 1 on successful completion
    Args: None
    Throws: NA
    Comments: The files that the sequences are printed to are named
              $outputDir/categorizable_reads_fwd.fastq
              $outputDir/categorizable_reads_rev.fastq
    See Also: NA

=head2 print_single_read_categories

    Title: print_single_read_categories
    Usage: $my_sample->print_single_read_categories();
    Function: Print the categorizable reads that have only one read in that
              category in FASTQ format
    Returns: 1 on successful completion
    Args: None
    Throws: NA
    Comments: The files that the sequences are printed to are named
              $outputDir/single_read_categories_fwd.fastq
              $outputDir/single_read_categories_rev.fastq
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
  
=head2 _parse_read

    Title: _parse_read
    Usage: $my_sample->_parse_read($fwd_read, $rev_read, $header_prefix);
    Function: Calls the local method _parse_fwd_read
    Returns: MTToolbox::Match::PE::NoOverlap object OR -1
    Args: -fwd_read => A Bio::Seq object in FASTQ format
          -rev_read => A Bio::Seq object in FASTQ format
          -header_prefix => The sample ID and the sequence count in that sample
    Throws: NA
    Comments: I added the header prefix to the front of the sequence ID for ease
              of use when I combine several samples.  It is also required by
              OTUpipe and for building OTU tables.
    See Also: _parse_fwd_read, _parse_rev_read
    
=head2 _parse_fwd_read

    Title: _parse_fwd_read
    Usage: $my_sample->_parse_fwd_read($seq, $quals_aref);
    Function: Attempts to match the raw fwd read to a regex defining the expected
              sequence pattern for forward reads.
    Returns: Array of values: ($is_matched, $fwd_tag, $fwd_linker, $fwd_primer, $fwd_amplicon, $fwd_tail)
    Args: -seq => a string representing the sequence of the raw fwd read
          -quals_aref => an array reference with the quality values of the fwd sequence
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 _parse_rev_read

    Title: _parse_rev_read
    Usage: $my_sample->_parse_rev_read($seq, $quals_aref);
    Function: Attempts to match the raw rev read to a regex defining the expected
              sequence pattern for forward reads.
    Returns: Array of values: ($is_matched, $rev_tag, $rev_linker, $rev_primer, $rev_amplicon, $rev_tail)
    Args: -seq => a string representing the sequence of the raw rev read
          -quals_aref => an array reference with the quality values of the rev sequence
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
