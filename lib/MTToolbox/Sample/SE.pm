package MTToolbox::Sample::SE;

use strict;
use warnings;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(carp croak);
use MyX::Generic 0.0.3;
use version; our $VERSION = qv('4.1.2');

use MTToolbox::Mismatch;
use MTToolbox::MoleculeTagCategory::SE;
use MTToolbox::Match::SE;
use MTToolbox::SampleSummary::SE;
use MTToolbox::Sample;
use base qw(MTToolbox::Sample);  # inherits from MTToolbox::Sample

{
    Readonly my $NEW_USAGE => q{ new( {id => ,
                                       barcode => ,
                                       fwd_linker => ,
                                       fwd_primer => ,
                                       rev_primer => ,
                                       fwd_mer_len => ,
                                       fwd_max_shifts => ,
                                       min_con_depth => ,
                                       diginorm_max => ,
                                       output_dir => ,
                                       reads_file => ,
                                       } ) };
    # Attributes #
    my %reads_file_of;
    
    # Setters #
    sub set_reads_file;
    
    # Getters #
    sub get_reads_file;
    
    # Others #
    sub categorize_by_MTs;
    sub _parse_read;
    sub _parse_fwd_read;
    
    
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
                               $arg_href->{reads_file},
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
        $reads_file_of{ident $new_obj} = $arg_href->{reads_file};
        
        # Set the MTToolbox::SampleSummary object
        my $summ_info_obj = MTToolbox::SampleSummary::SE->new();
        $summ_info_obj->set_sample_name($new_obj->get_name());
        $summ_info_obj->set_sample_id($arg_href->{id});
        $new_obj->_set_summ_info_obj($summ_info_obj);
        
        
        return $new_obj;
    }
    
    
    ###########
    # Setters #
    ###########
    sub set_reads_file {
        my ($self, $f) = @_;
        $reads_file_of{ident $self} = $f;
        return 1;
    }
    
    
    
    
    ###########
    # Getters #
    ###########
    sub get_reads_file {
        my ($self) = @_;
        return $reads_file_of{ident $self};
    }
    
    
    
    
    ##########
    # Others #
    ##########
    sub categorize_by_MTs {
        my ($self) = @_;
    
        ### Create either a SE or FastqMismatch object for each seq
        ### and attempt to add SE objects to a MoleculeTag object
        ### Also, store some summary info as I go through the sequences
        my $matched_count = 0;
        my $mismatched_count = 0;
        my $total_seqs = 0;
        my $mt_count = 0;
        
        # Read in the sequences from the fastq file
        my $seqs_in;
        eval {
            $seqs_in = BioUtils::FastqIO->new( {
                            stream_type => '<',
                            file => $self->get_reads_file()
                        } );
        };
        if ( my $e = MyX::Generic::Undef::Param->caught() ) {
            croak("$e->error()\n$e->trace()\n");
        }
        while ( my $seq_in = $seqs_in->get_next_seq() ) {        
            
            # header_prefix is a variable to be the prefix to the current sequence
            # that contains the information on the sample and sequence count
            my $header_prefix = $self->get_id() . "_" . $total_seqs;
            
            ### Try creating a MTToolbox::Match::SE object with the read
            # After this line $match is either -1 meaning it was not a match or it is a SE.
            my $match = $self->_parse_read( $seq_in, $header_prefix );
            
            ### Check if the MTToolbox::Match::SE object was successfully created.  If not create a Mismatch object
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
                my $tag = $match->get_fwd_tag()->get_seq();
                
                my $temp_MT_objs_href = $self->get_MT_objs_href();
                
                # Create a new MoleculeTag object if needed
                if ( ! defined $temp_MT_objs_href->{$tag} ) {
                    $mt_count++;
                    
                    # create the new MT object
                    my $mt_obj = MTToolbox::MoleculeTagCategory::SE->new({ tag => $tag });
                    
                    # save the new MT object in this sample
                    $temp_MT_objs_href->{$tag} = $mt_obj;  
                }
                
                # Add the MTToolbox::Match::SE to a MoleculeTag object if the number of reads
                # is less than diginorm_max (e.g. digital normalization)
                if ($temp_MT_objs_href->{$tag}->get_seq_count() <
                    $self->get_diginorm_max() )
                {
                    $temp_MT_objs_href->{$tag}->add_match($match);
                }
            }
            
            # Go to the next sequence and increment the total sequence counter
            $total_seqs++;
        }
        
        ### Store the summary information
        my $summ_info_obj = $self->get_summ_info_obj();
        $summ_info_obj->set_file_format('fastq');
        $summ_info_obj->set_match_count($matched_count);
        $summ_info_obj->set_mismatch_count($mismatched_count);
        $summ_info_obj->set_seq_count($total_seqs);
        $summ_info_obj->set_MT_count($mt_count);
        
        return 1;
    }
    
    sub _parse_read {
        my ($self, $fwd_read, $header_prefix) = @_;
        
        my ($is_fwd_match, $fwd_tag, $fwd_linker, $fwd_primer, $fwd_amplicon, $fwd_tail) =
            $self->_parse_fwd_read($fwd_read->get_seq(), $fwd_read->get_quals_str());
        
        if ( $is_fwd_match ) {
            # NOTE: When passing in the Id for the sequence the fwd_read->get_header() should be the same as the revRead->get_header().
            my $desc = $fwd_tag->get_seq() . " " . $fwd_read->get_header();
            my $matched_fastq_seq = MTToolbox::Match::SE->new({ id => $header_prefix,
                                                     desc => $fwd_read->get_header(),
                                                     fwd_tag => $fwd_tag,
                                                     fwd_linker => $fwd_linker,
                                                     fwd_primer => $fwd_primer,
                                                     fwd_amplicon => $fwd_amplicon,
                                                     fwd_tail => $fwd_tail
                                                   });
            return $matched_fastq_seq;
        }
        
        return -1; # Indicates an unmatched sequence
    }
    
    sub _parse_fwd_read {
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
            
            my $fwd_tag_quals_str = substr $quals_str, $index, length $1;
            $fwd_tag = BioUtils::FastqSeq->new( {seq => $1, quals_str => $fwd_tag_quals_str} );
            $index = length $1;
            
            my $fwd_linker_quals_str = substr $quals_str, $index, length $self->get_fwd_linker();
            $fwd_linker = BioUtils::FastqSeq->new( {
                                seq => $self->get_fwd_linker(),
                                quals_str => $fwd_linker_quals_str
                            } );
            $index += length($self->get_fwd_linker());
            
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
}

1;
__END__

#######
# POD #
#######
=head1 MTToolbox::Sample::SE

MTToolbox::Sample::SE - A class to store and handle samples with single end
                        sequences

=head1 VERSION

This documentation refers to MTToolbox::Sample::SE version 4.1.2.

=head1 Included Modules

    Class::Std::Utils
    List::MoreUtils qw(any)
    Readonly
    Carp qw(carp croak)
    MyX::Generic 1.0.0
    version
	MTToolbox::Mismatch
    MTToolbox::MoleculeTagCategory::SE
    MTToolbox::Match::SE
    MTToolbox::Sample

=head1 Inherit

    MTToolbox::Sample

=head1 SYNOPSIS

    # NOTE: MTToolbox::Sample::SE inherits from MTToolbox::Sample.  See
    # MTToolbox::Sample.pm for details on inherited methods.
    
    use MTToolbox::Sample::SE;
    my $my_sample = MTToolbox::Sample::SE->new({id => $id,
                                     barcode => $barcode,
                                     reads_file => $reads_file,
                                     fwd_linker => $fwd_linker,
                                     fwd_primer => $fwd_primer,
                                     rev_primer => $rev_primer,
                                     min_con_depth => $min_con_depth,
                                     diginorm_max => $diginorm_max,
                                     output_dir => $output_dir,
                                     });
    
    my $reads_file = $my_sample->get_reads_file();
    
    $my_sample->set_reads_file($file);
    
    # NOTE: Many of the major opperations are implemented in the parant class
    # (MTToolbox::Sample).  See MTToolbox::Sample.pm documentation for more
    # details.
    
    # The only major opperation that is coded in this class is:
    $my_sample->categorize_by_MTs();

=head1 DESCRIPTION

MTToolbox::Sample::SE is an object that stores sequence information that comes
from a distinct biological sample with single end (SE) reads.  This class
inherits from MTToolbox::Sample.  See the documentation in MTToolbox::Sample.pm
for important information on inherited methods.  

=head1 METHODS

=over

    new
    set_reads_file
    get_reads_file
    categorize_by_MTs
    _parse_read
    _parse_fwd_read
    
=back

=head1 METHODS DESCRIPTION

=head2 new
    
    Title: new
    Usage: MTToolbox::Sample::SE->new({id => $id,
                            barcode => $barcode,
                            reads_file => $reads_file,
                            fwd_linker => $fwd_linker,
                            fwd_primer => $fwd_primer,
                            rev_primer => $rev_primer,
                            fwd_mer_len => $fwd_mer_len,
                            fwd_max_shifts => $fwd_max_shifts,
                            min_con_depth => $min_con_depth,
                            diginorm_max => $diginorm_max,
                            output_dir => $output_dir,
                            });
    Function: Creates a new MTToolbox::Sample::SE object
    Returns: MTToolbox::Sample::SE
    Args: -id             => a string representing the sample id
          -barcode        => a string representing the sample barcode
          -readFile       => a string representing the path to the raw FASTQ
                             reads
          -fwd_linker     => a string representing the forward linker sequence
          -fwd_primer     => a string representing the forward primer sequence
          -rev_primer     => a string representing the reverse primer sequence
          -fwd_mer_len    => an int representing the fwd tag length
          -fwd_max_shifts => an int representing the max fwd frameshifts
          -min_con_depth  => an int representing the min number of reads needed
                             to make a consensus sequence
          -diginorm_max   => an int representing the max number of reads to
                             consider when making a consensus. (NA uses all seqs)
          -output_dir     => a string representing the path of the output
                             directory

=head2 set_reads_file

    Title: set_reads_file
    Usage: $my_sample->set_reads_file("/home/Scott/data/SE_reads.fastq");
    Function: Sets the reads file value
    Returns: 1 on successful completion
    Args: -readFile => a string representing the path to the raw FASTQ reads
    Throws: NA
    Comments: NA
    See Also: NA
    
=head2 get_reads_file
    
    Title: get_reads_file
    Usage: my $reads_file = $my_sample->get_reads_file();
    Function: Gets the reads file directory path
    Returns: String
    Args: None
    Throws: NA
    Comments: NA
    See Also: NA

=head2 categorize_by_MTs

    Title: categorize_by_MTs
    Usage: $my_sample->categorize_by_MTs();
    Function: Categorizes the raw SE sequences by their molecule tag(s)
    Returns: 1 on successful completion
    Args: None
    Throws: NA
    Comments: This method uses the local methods _parse_read and _parse_fwd_read
              to attempt to match each raw read to the defined pattern.  There
              is a regex in _parse_fwd_read that defines that pattern.  At some
              point this will hopefully turn into a user defined parameter.  If
              a pattern match is successful a MTToolbox::Match::SE object is
              created and added to the MTObjects hash inherited from
              MTToolbox::Sample.
    See Also: MTToolbox::Match.pm,
              MTToolbox::Match::SE.pm,
              MTToolbox::Match::SE::_parse_read,
              MTToolbox::Match::SE::_parse_fwd_read
    
=head2 _parse_read

    Title: _parse_read
    Usage: $my_sample->_parse_read($fwd_read, $header_prefix);
    Function: Calls the local method _parse_fwd_read
    Returns: MTToolbox::Match::SE object OR -1
    Args: -fwd_read => A Bio::Seq object in FASTQ format
          -header_prefix => The sample ID and the sequence count in that sample
    Throws: NA
    Comments: I added the header prefix to the front of the sequence ID for ease
              of use when I combine several samples.  It is also required by
              OTUpipe and for building OTU tables.
    See Also: _parse_fwd_read
    
=head2 _parse_fwd_read

    Title: _parse_fwd_read
    Usage: $my_sample->_parse_fwd_read($seq, $quals_aref);
    Function: Attempts to match the raw read to a regex defining the expected
              sequence pattern.
    Returns: Array of values: ($isMatched, $fwdTag, $fwd_linker, $fwdPrimer, $fwdAmplicon, $fwdTail)
    Args: -seq => a string representing the sequence of the raw read
          -quals_aref => an array reference with the quality values of the sequence
    Throws: NA
    Comments: The pattern also looks for the reverse primer sequence in the event
              that the SE read is longer than the amplicon
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
