package MTToolbox::Sample::PE::Overlap;

use strict;
use warnings;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(carp croak);
use MyX::Generic 0.0.3;
use version; our $VERSION = qv('4.1.2');

use MTToolbox::Mismatch;
use MTToolbox::Match::PE::Overlap;
use MTToolbox::MoleculeTagCategory::PE::Overlap;
use MTToolbox::SampleSummary::PE::Overlap;
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
                                       reads_file => ,
                                       flash_m => ,
                                       flash_M => ,
                                       flash_x => ,
                                       flash_p => ,
                                       flash_r => ,
                                       flash_f => ,
                                       flash_s => ,
                                       } ) };
    # Attributes #
    my %merged_reads_file_of;
    my %fwd_reads_file_of;
    my %rev_reads_file_of;
    my %rev_mer_len_of;
    my %rev_max_shifts_of;
    my %rev_linker_of;
    my %flash_params_of;
    
    # Setters #
    sub set_merged_reads_file;
    sub set_fwd_reads_file;
    sub set_rev_reads_file;
    sub set_rev_mer_len;
    sub set_rev_max_shifts;
    sub set_rev_linker;
    sub set_flash_params;
    
    # Getters #
    sub get_merged_reads_file;
    sub get_fwd_reads_file;
    sub get_rev_reads_file;
    sub get_rev_mer_len;
    sub get_rev_max_shits;
    sub get_rev_linker;
    sub get_flash_params;
    
    # Others #
    sub merge_with_flash;
    sub _get_flash_params_str;
    sub categorize_by_MTs;
    sub _parse_read;
    sub _parse_merged_read;
    sub _reverse_complement;
    
    
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
        $merged_reads_file_of{ident $new_obj} = $arg_href->{reads_file}
            if defined $arg_href->{reads_file};
        $fwd_reads_file_of{ident $new_obj} = $arg_href->{fwd_file}
            if defined $arg_href->{fwd_file};
        $rev_reads_file_of{ident $new_obj} = $arg_href->{rev_file}
            if defined $arg_href->{rev_file};
        $rev_mer_len_of{ident $new_obj} = $arg_href->{rev_mer_len};
        $rev_max_shifts_of{ident $new_obj} = $arg_href->{rev_max_shifts};
        $rev_linker_of{ident $new_obj} = $arg_href->{rev_linker};
        
        my %flash_params;
        $flash_params{m} = $arg_href->{flash_m} if ( defined $arg_href->{flash_m} );
        $flash_params{M} = $arg_href->{flash_M} if ( defined $arg_href->{flash_M} );
        $flash_params{x} = $arg_href->{flash_x} if ( defined $arg_href->{flash_x} );
        $flash_params{p} = $arg_href->{flash_p} if ( defined $arg_href->{flash_p} );
        $flash_params{r} = $arg_href->{flash_r} if ( defined $arg_href->{flash_r} );
        $flash_params{f} = $arg_href->{flash_f} if ( defined $arg_href->{flash_f} );
        $flash_params{s} = $arg_href->{flash_s} if ( defined $arg_href->{flash_s} );
        $flash_params_of{ident $new_obj} = \%flash_params;
        
        # Set the MTToolbox::SampleSummary object
        my $summ_info_obj = MTToolbox::SampleSummary::PE::Overlap->new();
        $summ_info_obj->set_sample_name($new_obj->get_name());
        $summ_info_obj->set_sample_id($arg_href->{id});
        $new_obj->_set_summ_info_obj($summ_info_obj);
        
        return $new_obj;
    }
    
    
    ###########
    # Setters #
    ###########
    sub set_merged_reads_file {
        my ($self, $f) = @_;
        $merged_reads_file_of{ident $self} = $f;
        return 1;
    }
    
    sub set_fwd_reads_file {
        my ($self, $f) = @_;
        $fwd_reads_file_of{ident $self} = $f;
        return 1;
    }
    
    sub set_rev_reads_file {
        my ($self, $f) = @_;
        $rev_reads_file_of{ident $self} = $f;
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
    
    sub set_flash_params {
        my ($self, $flash_params) = @_;
        $flash_params_of{ident $self} = $flash_params;
        return 1;
    }
    
    
    
    
    ###########
    # Getters #
    ###########
    sub get_merged_reads_file {
        my ($self) = @_;
        return $merged_reads_file_of{ident $self};
    }
    
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
    
    sub get_flash_params {
        my ($self) = @_;
        return $flash_params_of{ident $self};
    }
    
    
    
    
    ##########
    # Others #
    ##########
    sub merge_with_flash {
        my ($self, $flash_params_file) = @_;
        
        my $output_dir = $self->get_output_dir();
        
        # Run the merge opperation
        my $command = "flash " .
                      $self->get_fwd_reads_file() . " " .
                      $self->get_rev_reads_file() . " " .
                      $self->_get_flash_params_str($flash_params_file) .
                      " -o  merged_seqs " .
                      " -d $output_dir";
        warn "FLASH command: $command\n";
        system($command);
        
        # Save the merged file
        $self->set_merged_reads_file($output_dir . "/merged_seqs.extendedFrags.fastq");
        
        return 1;
    }
    
    sub _get_flash_params_str {
        my ($self, $file) = @_;
        
        my $param_str = '';
        
        my $flash_params_href = $flash_params_of{ident $self};
        if ( ! defined $file or ! -f $file ) {
            $param_str .= "-m " . $flash_params_href->{m} . " "
                if ( defined $flash_params_href->{m} );
            $param_str .= "-M " . $flash_params_href->{M} . " "
                if ( defined $flash_params_href->{M} );
            $param_str .= "-x " . $flash_params_href->{x} . " "
                if ( defined $flash_params_href->{x} );
            $param_str .= "-p " . $flash_params_href->{p} . " "
                if ( defined $flash_params_href->{p} );
            $param_str .= "-r " . $flash_params_href->{r} . " "
                if ( defined $flash_params_href->{r} );
            $param_str .= "-f " . $flash_params_href->{f} . " "
                if ( defined $flash_params_href->{f} );
            $param_str .= "-s " . $flash_params_href->{s} . " "
                if ( defined $flash_params_href->{s} );
        }
        else {
            open my $PAR_FH, "<", $file or croak("Cannot open file: $file\nERROR: $!\n");
            
            $param_str = '';
            foreach my $line (<$PAR_FH>) {
                
                if ( $line =~ m/(\S+)\s+(\S+)/ ) {
                    $param_str .= " -" . $1 . " " . $2;
                }
                else {
                    carp("Bad merge parameter format: $line");
                }
            }
            
            close($PAR_FH);
        }
        
        return $param_str;
    }
    
    sub categorize_by_MTs {
        my ($self) = @_;
    
        ### Create either a PE::Overlap or Mismatch object for each seq
        ### and attempt to add PE::Overlap objects to a MoleculeTag object
        ### Also, store some summary info as I go through the sequences
        my $matched_count = 0;
        my $mismatched_count = 0;
        my $merged_seq_count = 0;
        my $mt_count = 0;
        
        # Read in the sequences from the fastq file
        my $seqs_in;
        eval {
            $seqs_in = BioUtils::FastqIO->new( {
                            stream_type => '<',
                            file => $self->get_merged_reads_file()
                        } );
        };
        if ( my $e = MyX::Generic::Undef::Param->caught() ) {
            croak("$e->error()\n$e->trace()\n");
        }
        elsif ($@) {  # warn of an uncaught error
            warn $@;
        }
        
        # Start reading the sequences one at a time.
        while ( my $seq_in = $seqs_in->get_next_seq() ) {        
            
            # header_prefix is a variable to be the prefix to the current sequence
            # that contains the information on the sample and sequence count
            my $header_prefix = $self->get_id() . "_" . $merged_seq_count;
            
            ### Try creating a MTToolbox::Match::PE::Overlap object with the read
            # After this line $match is either -1 meaning it was not a match or it is a PE::Overlap.
            my $match = $self->_parse_read( $seq_in, $header_prefix );
            
            ### Check if the MTToolbox::Match::PE::Overlap object was successfully created.
            # If not create a Mismatch object
            if ( $match == -1 ) {
                $mismatched_count++;
                push @{$self->get_mismatch_seqs_aref()}, MTToolbox::Mismatch->new({
                                                                       id => $header_prefix,
                                                                       desc => $seq_in->get_header(),
                                                                       seq => $seq_in->get_seq(),
                                                                       quals_str => $seq_in->get_quals_str()
                                                                       });
            }
            else {
                $matched_count++;
                
                # create the tag so I can look for it in the hash or create a new hash entry
                my $tag = $match->get_fwd_tag()->get_seq() .
                          "-" .
                          $match->get_rev_tag()->get_seq();
                
                my $temp_MT_objs_href = $self->get_MT_objs_href();
                
                # Create a new MoleculeTag object if needed
                if ( ! defined $temp_MT_objs_href->{$tag} ) {
                    $mt_count++;
                    
                    # create the new MT object
                    my $mt_obj = MTToolbox::MoleculeTagCategory::PE::Overlap->new({tag => $tag});
                    
                    # save the new MT object in this sample
                    $temp_MT_objs_href->{$tag} = $mt_obj;
                }
                
                # Add the MTToolbox::Match::PE::Overlap object if the number of reads
                # is less than diginorm_max (e.g. digital normalization)
                if ( $temp_MT_objs_href->{$tag}->get_seq_count() <
                    $self->get_diginorm_max() )
                {
                    $temp_MT_objs_href->{$tag}->add_match($match);  
                }
            }
            
            # Go to the next sequence and increment the total sequence counter
            $merged_seq_count++;
        }
        
        ### Store the summary information
        my $summ_info_obj = $self->get_summ_info_obj();
        $summ_info_obj->set_file_format('fastq');
        $summ_info_obj->set_match_count($matched_count);
        $summ_info_obj->set_mismatch_count($mismatched_count);
        $summ_info_obj->set_MT_count($mt_count);
        $summ_info_obj->set_merged_count($merged_seq_count);
        
        # optional summary info
        if ( my $fwd_file = $self->get_fwd_reads_file() ) {
            # I assume the fwd reads file is in fastq format
            my $command = "wc -l $fwd_file " . q{| awk '{print $1/4}' };
            my $total_seqs = `$command`;
            
            chomp $total_seqs;
            $summ_info_obj->set_seq_count($total_seqs);
            $summ_info_obj->set_not_merged_count($total_seqs - $merged_seq_count);
        }
        
        return 1;
    }
    
    sub _parse_read {
        my ($self, $merged_read, $header_prefix) = @_;
        
        my (
            $is_match,
            $fwd_tag,
            $fwd_linker,
            $fwd_primer,
            $fwd_amplicon,
            $rev_primer,
            $rev_linker,
            $rev_tag,
            ) =
            $self->_parse_merged_read( $merged_read->get_seq(), $merged_read->get_quals_str() );
        
        if ( $is_match ) {
            # NOTE: When passing in the Id for the sequence the
            # fwd_read->get_header() should be the same as the
            # rev_read->get_header().
            my $desc = $fwd_tag->get_seq() . "-" . $rev_tag->get_seq() . " ".
                        $merged_read->get_header();
            my $matched_fastq_seq = MTToolbox::Match::PE::Overlap->new({
                                                            id => $header_prefix,
                                                            desc => $desc,
                                                            fwd_tag => $fwd_tag,
                                                            fwd_linker => $fwd_linker,
                                                            rev_linker => $rev_linker,
                                                            fwd_primer => $fwd_primer,
                                                            fwd_amplicon => $fwd_amplicon,
                                                            rev_primer => $rev_primer,
                                                            rev_tag => $rev_tag,
                                                            });
            return $matched_fastq_seq;
        }
        
        return -1; # Indicates an unmatched sequence
    }
    
    sub _parse_merged_read {
        my ($self, $merged_read, $quals_str) = @_;
        
        # reverse complement the reverse primer sequnece
        my $temp_rev_primer = _reverse_complement($self->get_rev_primer());
        my $temp_rev_linker = _reverse_complement($self->get_rev_linker());
        
        # The regex for the merged read
        #my $regex = "((^[ATCG]{" . $self->get_fwd_mer_len() . "})(A|TA)?)" .    # fwd_tag
        #            $self->get_fwd_linker() .                                   # fwd_linker
        #            "(" . $self->get_fwd_primer() . ")" .                       # fwd_primer
        #            "(\\S+)" .                                                  # amplicon
        #            "(" . $temp_rev_primer . ")" .                              # rev_primer
        #            $temp_rev_linker .                                          # rev_linker
        #            "((T|TA)?([ATCG]{" . $self->get_rev_mer_len() . "}))" .     # rev_tag
        #            q{$}                                                        # END
        #            ;
        
        my $fwd_min = $self->get_fwd_mer_len();
        my $fwd_max = $fwd_min + $self->get_fwd_max_shifts();
        my $rev_min = $self->get_rev_mer_len();
        my $rev_max = $rev_min + $self->get_rev_max_shifts();
        
        my $regex = "(^[ATCG]{" . $fwd_min . "," . $fwd_max . "})" .            # fwd_tag
                    $self->get_fwd_linker() .                                   # fwd_linker
                    "(" . $self->get_fwd_primer() . ")" .                       # fwd_primer
                    "(\\S+)" .                                                  # amplicon
                    "(" . $temp_rev_primer . ")" .                              # rev_primer
                    $temp_rev_linker .                                          # rev_linker
                    "([ATCG]{" . $rev_min . "," . $rev_max . "})" .             # rev_tag
                    q{$}                                                        # END
        ;
        
        # Return Values
        my $is_matched = 0; # This is a boolean to help me identify in the calling method if this worked.
        my ($fwd_tag, $fwd_linker, $fwd_primer, $fwd_amplicon, $rev_primer, $rev_linker, $rev_tag);
        
        if ( $merged_read =~ m/$regex/i ) {
            my $index = 0;  # this keeps track of the end of where I am so I can index the quals array correctly
            
            # fwd_tag
            my $fwd_tag_quals_str = substr $quals_str, $index, length $1;
            $fwd_tag = BioUtils::FastqSeq->new( {seq => $1, quals_str => $fwd_tag_quals_str} );
            $index = length $1;
            
            # fwd_linker
            my $fwd_linker_quals_str = substr $quals_str, $index, length $self->get_fwd_linker();
            $fwd_linker = BioUtils::FastqSeq->new( {seq => $self->get_fwd_linker(), quals_str => $fwd_linker_quals_str} );
            $index += length($self->get_fwd_linker());
            
            # fwd_primer
            my $fwd_primer_quals_str = substr $quals_str, $index, length $2;
            $fwd_primer = BioUtils::FastqSeq->new( {seq => $2, quals_str => $fwd_primer_quals_str} );
            $index += length $2;
            
            # amplicon
            my $fwd_amplicon_quals_str = substr $quals_str, $index, length $3;
            $fwd_amplicon = BioUtils::FastqSeq->new( {seq => $3, quals_str => $fwd_amplicon_quals_str} );
            $index += length $3;
            
            # rev_primer
            my $rev_primer_quals_str = substr $quals_str, $index, length $4;
            $rev_primer = BioUtils::FastqSeq->new( {seq => $4, quals_str => $rev_primer_quals_str} );
            $index += length $4;
            
            # rev_linker
            my $rev_linker_quals_str = substr $quals_str, $index, length $self->get_rev_linker();
            $rev_linker = BioUtils::FastqSeq->new( {seq => $temp_rev_linker, quals_str => $rev_linker_quals_str} );
            $index += length($self->get_rev_linker());
            
            # rev_tag
            my $rev_tag_quals_str = substr $quals_str, $index, length $5;
            $rev_tag = BioUtils::FastqSeq->new( {seq => $5, quals_str => $rev_tag_quals_str} );
            $index += length $5;
            
            # This is a boolean signifying that the sequence was a match to the given regex
            $is_matched = 1;
        }
        else {
            $is_matched = 0;
        }
        
        return ($is_matched, $fwd_tag, $fwd_linker, $fwd_primer, $fwd_amplicon, $rev_primer, $rev_linker, $rev_tag);
    }
    
    sub _reverse_complement {
        my ($seq) = @_;
        
        # compelement all the bases
        $seq =~ tr/ACGTacgt/TGCAtgca/;
        
        # reverse all the characters
        $seq = scalar reverse $seq;
        
        # change brackets (this is because when we reverse it they are still backward)
        $seq =~ tr/\[\]/\]\[/;
        
        return $seq;
    }
}

1;
__END__

#######
# POD #
#######
=head1 MTToolbox::Sample::PE::Overlap

MTToolbox::Sample::PE::Overlap - A class to store and handle samples with
paired-end (PE) sequences that overlap and have been merged.

=head1 VERSION

This documentation refers to MTToolbox::Sample::PE::Overlap version 4.1.2.

=head1 Included Modules


    Class::Std::Utils
    List::MoreUtils qw(any)
    Readonly
    Carp qw(carp croak)
    MyX::Generic 1.0.0
    version
	MTToolbox::Mismatch
    MTToolbox::Match::PE::Overlap
    MTToolbox::MoleculeTagCategory::PE::Overlap
    MTToolbox::SampleSummary::PE::Overlap
    MTToolbox::Sample

=head1 Inherit

    MTToolbox::Sample

=head1 SYNOPSIS

    # NOTE: MTToolbox::Sample::PE::Overlap inherits from MTToolbox::Sample.
    # See MTToolbox::Sample.pm for details on inherited methods.
    
    # Either a merged file OR a fwd and rev file must be provided.
    use MTToolbox::Sample::PE::Overlap
    my $my_sample = MTToolbox::Sample::PE::Overlap->new({
                                                id => $id,
                                                barcode => $barcode,
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
                                                merged_reads_file => $merged_reads_file,
                                                fwd_reads_file => $fwd_reads_file,
                                                rev_reads_file => $rev_reads_file,
                                                flash_m => $flash_m,
                                                flash_M => $flash_M,
                                                flash_x => $flash_x,
                                                flash_p => $flash_p,
                                                flash_r => $flash_r,
                                                flash_f => $flash_f,
                                                flash_s => $flash_s,
                                            });
    
    my $merged_reads_file = $my_sample->get_merged_reads_file();
    my $fwd_reads_file = $my_sample->get_fwd_reads_file();
    my $rev_reads_file = $my_sample->get_rev_reads_file();
    
    $my_sample->set_merged_reads_file($file);
    $my_sample->set_fwd_reads_file($file);
    $my_sample->set_rev_reads_file($file);
    
    # NOTE: Many of the major opperations are implemented in the parant class
    # (MTToolbox::Sample).  See MTToolbox::Sample.pm documentation for more
    # details.
    
    # To merge reads using FLASH
    # (NOTE: the flash exe must be in your PATH env)
    $my_sample->merge_with_flash();
    
    # The major opperation that is coded in this class is:
    $my_sample->categorize_by_MTs();

=head1 DESCRIPTION

MTToolbox::Sample::PE::Overlap is an object that stores sequence information that comes from
a distinct biological sample with paired-end (PE) reads.  This class inherits
from MTToolbox::Sample.  See the documentation in MTToolbox::Sample.pm for
important information on inherited methods.  

=head1 METHODS

=over

    new
    set_merged_reads_file
    set_fwd_reads_file
    set_rev_reads_file
    set_rev_mer_len
    set_rev_max_shifts
    set_rev_linker
    set_flash_params
    get_merged_reads_file
    get_fwd_reads_file
    get_rev_reads_file
    get_rev_mer_len
    get_rev_max_shifts
    get_rev_linker
    get_flash_params
    merge_with_flash
    _get_flash_params_str
    categorize_by_MTs
    _parse_read
    _parse_merged_read
    _reverse_complement
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: MTToolbox::Sample::PE::Overlap->new({
                                        id => $id,
                                        barcode => $barcode,
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
                                        merged_reads_file => $merged_reads_file,
                                        fwd_reads_file => $fwd_reads_file,
                                        rev_reads_file => $rev_reads_file,
                                        flash_m => $flash_m,
                                        flash_M => $flash_M,
                                        flash_x => $flash_x,
                                        flash_p => $flash_p,
                                        flash_r => $flash_r,
                                        flash_f => $flash_f,
                                        flash_s => $flash_s,
                                    });
    Function: Creates a new MTToolbox::Sample::SE::Overlap object
    Returns: MTToolbox::Sample::SE::Overlap
    Args: -id             => a string representing the sample id
          -barcode        => a string representing the sample barcode
          -fwd_linker     => a string representing the forward linker sequence
          -rev_linker     => a string representing the reverse linker sequence
          -fwd_primer     => a string representing the forward primer sequence
          -rev_primer     => a string representing the reverse primer sequence
          -fwd_mer_len    => an int representing the length of the fwd tag
          -rev_mer_len    => an int representing the length of the rev tag
          -fwd_max_shifts => an int representing the max fwd frameshifts
          -rev_max_shifts => an int representing the max rev frameshifts
          -min_con_depth  => an int representing the min number of reads needed
                             to make a consensus sequence
          -diginorm_max   => an int representing the max number of reads to
                             consider when making a consensus. (NA uses all seqs)
          -output_dir     => a string representing the path of the output
                             directory
          -merged_reads_file => a string representing the path to a fastq file
                                of merged reads
          -fwd_reads_file => a string representing the path to a fastq file
                                of fwd reads
          -rev_reads_file => a string representing the path ot a fastq file
                                of rev reads
          -flash_m => flash min overlap parameter
          -flash_M => flash max overlap parameter
          -flash_x => flash ratio parameter
          -flash_p => flash phred offset parameter
          -flash_r => flash read length parameter
          -flash_f => flash fragment length parameter
          -flash_s => flash standard deviation parameter
    Throws: NA
    Comments: Either a merged file or a fwd and rev file pair must be provided.
    See Also: MTToolbox::Sample.pm

=head2 set_merged_reads_file

    Title: set_merged_reads_file
    Usage: $my_sample->set_merged_reads_file("/home/Scott/data/PE_merged_reads.fastq");
    Function: Sets the reads file value
    Returns: 1 on successful completion
    Args: -readFile => a string representing the path to the raw FASTQ reads
    Throws: NA
    Comments: The file should have merged paired-end reads
    See Also: NA

=head2 set_fwd_reads_file

    Title: set_fwd_reads_file
    Usage: $my_sample->set_fwd_reads_file("/home/Scott/data/PE_fwd_reads.fastq");
    Function: Sets the reads file value
    Returns: 1 on successful completion
    Args: -readFile => a string representing the path to the raw FASTQ reads
    Throws: NA
    Comments: The file should have the fwd read of paired-end reads
    See Also: NA

=head2 set_rev_reads_file

    Title: set_rev_reads_file
    Usage: $my_sample->set_rev_reads_file("/home/Scott/data/PE_rev_reads.fastq");
    Function: Sets the reads file value
    Returns: 1 on successful completion
    Args: -readFile => a string representing the path to the raw FASTQ reads
    Throws: NA
    Comments: The file should have the rev read of paired-end reads
    See Also: NA
    
=head2 set_rev_mer_len

    Title: set_rev_mer_len
    Usage: $my_sample->set_rev_mer_len(9);
    Function: Sets the length of the rev random mer
    Returns: 1 on successful completion
    Args: -rev_mer_len => an int representing the lenght of the rev random mer
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 set_rev_max_shifts

    Title: set_rev_max_shifts
    Usage: $my_sample->set_rev_max_shifts(5);
    Function: Sets the max number of frameshifts in the rev tag
    Returns: 1 on successful completion
    Args: -rev_max_shifts => an int representing the max rev frameshifts
    Throws: NA
    Comments: NA
    See Also: NA
   
=head2 set_rev_linker
    
    Title: set_rev_linker
    Usage: $my_sample->set_rev_linker("CAGT");
    Function: Sets the reverse linker value
    Returns: 1 on successful completion
    Args: -rev_linker => a string representing the sample reverse linker sequence
    Throws: NA
    Comments: NA
    See Also: NA

=head2 set_flash_params
    
    Title: set_flash_params
    Usage: $my_sample->set_flash_params($params_href);
    Function: Sets the flash parameters as and href
    Returns: 1 on successful completion
    Args: -params_href => flash parameter stored in a hash
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_merged_reads_file
    
    Title: get_merged_reads_file
    Usage: my $merged_reads_file = $my_sample->get_merged_reads_file();
    Function: Gets the reads file directory path
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_fwd_reads_file
    
    Title: get_fwd_reads_file
    Usage: my $fwd_reads_file = $my_sample->get_fwd_reads_file();
    Function: Gets the reads file directory path for the fwd reads
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 get_rev_reads_file
    
    Title: get_rev_reads_file
    Usage: my $rev_reads_file = $my_sample->get_rev_reads_file();
    Function: Gets the reads file directory path for the rev reads
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_rev_mer_len

    Title: get_rev_mer_len
    Usage: my $rev_mer_len = $my_sample->get_rev_mer_len();
    Function: Gets the length of the rev random mer
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_rev_max_shifts

    Title: get_rev_max_shifts
    Usage: my $rev_max_shifts = $my_sample->get_rev_max_shifts();
    Function: Get the maximum number of rev tag frameshifts
    Returns: int
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_rev_linker

    Title: get_rev_linker
    Usage: my $rev_linker = $my_sample->get_rev_linker();
    Function: Gets the sample reverse linker sequence
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_flash_params

    Title: get_flash_params
    Usage: my $params_href = $my_sample->get_flash_params();
    Function: Gets the flash parameters in an href
    Returns: Href
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 merge_with_flash
    
    Title: merge_with_flash
    Usage: $my_sample->merge_with_flash("flash_params.txt");
    Function: Runs FLASH to merge PE reads that overlap
    Returns: 1 on successful completion
    Args: file => the FLAH parameters file as described below
    Throws: NA
    Comments: There are several output files generate by FLASH.  The main output
              file that contains the merged reads will be called
              merged_seqs.extendedFrags.fastq.  The default option to get the
              parameters that flash uses is by includeing them when the
              MTToolbox::Sample::PE::Overlap object is create.  In MTToolbox
              they are included in the config.xml file and passed down to where
              the object is created.  An alternative is to use a parameters file
              that can specifically set the FLASH parameters.  It should have
              the parameter name (as would be entered on the FLASH command)
              followed by a tab with the value.  For example (there is a tab):
              m 10
              M 250
              x .25
    See Also: NA
    
=head2 _get_flash_params_str
    
    Title: _get_flash_params_str
    Usage: my $flash_params_str = _get_flash_params_str();
    Function: Gets the parameter values for FLASH from the input parameters
              file in the form of a string
    Returns: String
    Args: None
    Throws: NA
    Comments: This is a private method and should not be used outside of
              MTToolbox::Sample::PE::Overlap.  The defualt is to use the
              flash_params attribute in the MTToolbox::Sample::PE::Overlap
              object that stores the flash parameters in hash format.  Another
              older option is to pass this function a file with the flash
              parameters.
    See Also: NA

=head2 categorize_by_MTs

    Title: categorize_by_MTs
    Usage: $my_sample->categorize_by_MTs();
    Function: Categorizes the raw PE::Overlap sequences by their molecule tag(s)
    Returns: 1 on successful completion
    Args: None
    Throws: NA
    Comments: This method uses the local methods _parse_read and
              _parse_merged_read to attempt to match each raw read to the
              defined pattern.  There is a regex in _parse_merged_read that
              defines that pattern.  At some point this will hopefully turn into
              a user defined parameter.  If a pattern match is successful a
              MTToolbox::Match::PE::Overlap object is created and added to the
              MTToolbox::MTObjects hash inherited from MTToolbox::Sample.
    See Also: MTToolbox::Match.pm,
              MTToolbox::Match::PE::Overlap.pm,
              MTToolbox::Match::PE::Overlap::_parse_read,
              MTToolbox::Match::PE::Overlap::_parse_merged_read
    
=head2 _parse_read

    Title: _parse_read
    Usage: $my_sample->_parse_read($fwd_read, $header_prefix);
    Function: Calls the local method _parse_merged_read
    Returns: MTToolbox::Match::PE::Overlap object OR -1
    Args: -fwd_read => A Bio::Seq object in FASTQ format
          -header_prefix => The sample ID and the sequence count in that sample
    Throws: NA
    Comments: I added the header prefix to the front of the sequence ID for ease
              of use when I combine several samples.  It is also required by
              OTUpipe and for building OTU tables.
    See Also: _parse_merged_read
    
=head2 _parse_merged_read

    Title: _parse_merged_read
    Usage: $my_sample->_parse_merged_read($seq, $quals_aref);
    Function: Attempts to match the raw read to a regex defining the expected
              sequence pattern.
    Returns: Array of values: (
                                $isMatched,
                                $fwdTag,
                                $fwd_linker,
                                $fwdPrimer,
                                $fwdAmplicon,
                                $rev_primer,
                                $rev_linker,
                                $rev_tag,
                                )
    Args: -seq => a string representing the sequence of the raw read
          -quals_aref => an array reference with the quality values of the
                         sequence
    Throws: NA
    Comments: NA
    See Also: NA

=head2 _reverse_complement

    Title: _reverse_complement
    Usage: _reverse_complement($seq);
    Function: Reverse complements a DNA regex (eg the fwd and rev primer regexs)
    Returns: A string (in regex form)
    Args: -seq => a string of DNA characters and regex characters
    Throws: NA
    Comments: This is a PRIVATE method.  It should not be used.
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
