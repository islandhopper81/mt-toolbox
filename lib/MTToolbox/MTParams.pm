package MTToolbox::MTParams;

use warnings;
use strict;

use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(carp croak);
use XML::Simple;
use Scalar::Util qw(looks_like_number);
use Exception::Class;
use MyX::Generic 0.0.3;
use MTToolbox::MyX::MTParams 4.1.2;
use version; our $VERSION = qv('4.1.2');

{
	# Usage statement
	Readonly my $NEW_USAGE => q{ new({[xml_file => ],
                                      [href => ],
                                      } ) };

	# Required global tags in xml configuration file
	Readonly::Hash my %required_global_tags => map { $_ => 1 } qw(
												output_root
												fwd_primer
												rev_primer
												seq_type
												sample
												parallelize
											);
	
	# Optional global tags in xml configuration file
	Readonly::Hash my %optional_global_tags => map { $_ => 1 } qw(
                                                split_samples_root_dir
												bin_sample_exec
												merge_params_file
												fwd_mer_len
												rev_mer_len
												fwd_max_shifts
												rev_max_shifts
												fwd_linker
												rev_linker
                                                min_con_depth
                                                diginorm_max
                                                flash_params
                                                qc_params
												keep_tmp_files
												con_algo
											);

	# Required for each sample in the configuration file
	Readonly::Hash my %required_sample_tags => map { $_ => 1 } qw(
													   sample_id
													);
	
	# Optional for each sample in the configuration file
	Readonly::Hash my %optional_sample_tags => map { $_ => 1 } qw(
													   barcode
													   seq_type
													   fwd_primer
													   rev_primer
													   fwd_mer_len
													   rev_mer_len
													   fwd_max_shifts
													   rev_max_shifts
													   fwd_linker
													   rev_linker
                                                       min_con_depth
                                                       diginorm_max
													   output_dir
													   fwd_file
													   rev_file
													   merged_file
													   flash_params
													   qc_params
													   con_algo
													);
	
	# Optional flash params tags in xml configuration file
	Readonly::Hash my %optional_flash_params_tags => map { $_ => 1 } qw(m
																		M
																		x
																		p
																		r
																		f
																		s
																	);

    # Optional QC params tags in xml configuration file
    Readonly::Hash my %optional_qc_params_tags => map { $_ => 1 } qw(
		con_trim_to_base
        con_min_len
		con_min_avg_qual
		con_allow_gaps
		con_allow_ambig
		con_min_c_score
		SRC_trim_to_base
		SRC_min_len
		SRC_min_avg_qual
		SRC_allow_gaps
		SRC_allow_ambig
    );
	
	# Some DEFAULT parameters
	Readonly::Scalar my $FWD_MER_LEN => 8;
	Readonly::Scalar my $REV_MER_LEN => 5;
	Readonly::Scalar my $FWD_MAX_SHIFTS => 5;
	Readonly::Scalar my $REV_MAX_SHIFTS => 5;
	Readonly::Scalar my $PARALLELIZE => "Y";  # Yes, use multiple processors
	Readonly::Scalar my $KEEP_TMP_FILES => "N";  # No, do NOT keep
	Readonly::Scalar my $MIN_CON_DEPTH => 2;
	Readonly::Scalar my $DIGINORM_MAX => "NA";  # NA == all reads are kept
	Readonly::Scalar my $CON_ALGO => "NoMSA";
	# NOTE: Some parameters (fwd_primer, rev_primer, fwd_linker, rev_linker, and
	#       seq_type) MUST be specifically defined by the user.  If they are
	#       not set then errors are thrown.
	# NOTE: FLASH parameters defaults are set by FLASH itself
	# NOTE: QC parameters defaults are set by perl fitler scripts

    # Required Global Attributes #
	my %params_href_of;

    # Subroutines #
	sub get_params_href;
	sub set_params_href;
    sub check_params;
    sub check_sample_entry;
    sub print_xml_file;
	
	sub _get_required_global_tags;
	sub _get_optional_global_tags;
	sub _get_required_sample_tags;
	sub _get_optional_sample_tags;
	sub _get_optional_flash_params_tags;
    sub _get_optional_qc_params_tags;
	
    sub _is_required_tag_present;
    sub _is_required_tag_defined;
    sub _is_optional_tag_present;
    sub _is_optional_tag_defined;
	
	sub _check_samples;
    sub _check_bin_sample_exec;
    sub _check_split_samples_root_dir;
    sub _check_output_root;
    sub _check_fwd_linker;
	sub _check_rev_linker;
    sub _check_fwd_primer;
    sub _check_rev_primer;
    sub _check_fwd_mer_len;
    sub _check_rev_mer_len;
	sub _check_fwd_max_shifts;
	sub _check_rev_max_shifts;
    sub _check_seq_type;
    sub _check_parallelize;
	sub _check_keep_tmp_files;
	sub _check_min_con_depth;
	sub _check_dignorm_max;
	sub _check_con_algo;
	sub _check_merge_params_file;
    sub _check_fastq_file;
	
	sub _check_min_overlap;
	sub _check_max_overlap;
	sub _check_ratio;
	sub _check_phred_offset;
	sub _check_read_len;
	sub _check_frag_len;
	sub _check_sd;
	
	sub _check_trim_to_base;
	sub _check_min_len;
	sub _check_min_avg_qual;
	sub _check_allow_gaps;
	sub _check_allow_ambig;
	sub _check_min_c_score;
	
    sub _translate_seq_type;
    
    
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
		
		# attributes
        if ( defined $arg_href ) {
            $new_obj->set_params_href($arg_href);
        }
        
        return $new_obj;
    }
    
    ###############
	# Subroutines #
	###############
    sub set_params_href {
        my ($self, $arg_href) = @_;
        
        # if an href is passed in use this as the href.
        if ( defined $arg_href->{href} ) {
            $params_href_of{ident $self} = $arg_href->{href};
        }
        elsif ( defined $arg_href->{xml_file} ) {
            $params_href_of{ident $self} =
                XML::Simple->new()->XMLin($arg_href->{xml_file});
        }
        else {
            my %href = ();
            $params_href_of{ident $self} = \%href;
        }
        
        return 1;
    }
	
	sub get_params_href {
		my ($self) = @_;
		return $params_href_of{ident $self};
	}
    
    sub check_params {
		my ($self, $data_href) = @_;
		
		# Check for:
		#   - Unknown tag values
		#   - required parameters
		#   - logic between seq_type and file values
		#   - directories exist
		#   - files exist
        
        # NOTE:
        # Errors throw in the methods called from this method are automatically
        # re-thrown.  they can be caught and handled in whom ever calls
        # check_params().
		
        # If no params href is provided the one stored in this object is used
        if ( ! defined $data_href ) {
            $data_href = $self->get_params_href();
        }
		
		# Check for required tags and their values in the TOP level of the XML file
		# This is not checking for required tags and values in the samples tag
		foreach my $tag ( keys %required_global_tags ) {
			_is_required_tag_defined($tag, $data_href);
		}
		
		# Check for tags that are unknown
		foreach my $tag ( keys %{$data_href} ) {
			if ( ! defined $required_global_tags{$tag} and
				 ! defined $optional_global_tags{$tag} ) {
					MTToolbox::MyX::MTParams::UnrecognizedTag->throw(
						error => "Unrecognized parameter",
						tag_name => $tag,
                );
			}
		}
		
		# Check for unknown FLASH paramters
		if ( defined $data_href->{flash_params} ) {
			foreach my $tag ( keys %{$data_href->{flash_params}} ) {
				if ( ! defined $optional_flash_params_tags{$tag} ) {
					MTToolbox::MyX::MTParams::UnrecognizedTag->throw(
                        error => "Unrecognized FLASH parameter",
                        tag_name => $tag,
                    );
				}
			}
		}
        
        # Check for unknown QC parameters
        if ( defined $data_href->{qc_params} ) {
			foreach my $tag ( keys %{$data_href->{qc_params}} ) {
				if ( ! defined $optional_qc_params_tags{$tag} ) {
					MTToolbox::MyX::MTParams::UnrecognizedTag->throw(
                        error => "Unrecognized QC parameter",
                        tag_name => $tag,
                    );
				}
			}
		}
		
		# Check the global variables
		_check_bin_sample_exec($data_href);
		_check_split_samples_root_dir($data_href);
		_check_output_root($data_href);
		_check_fwd_primer($data_href);
		_check_rev_primer($data_href);
		_check_seq_type($data_href);
		_check_parallelize($data_href);
		
		# Check optional global varaible
		_check_fwd_linker($data_href);
		_check_rev_linker($data_href);
		_check_merge_params_file($data_href);
		_check_fwd_mer_len($data_href);
		_check_rev_mer_len($data_href);
		_check_fwd_max_shifts($data_href);
		_check_rev_max_shifts($data_href);
        _check_min_con_depth($data_href);
        _check_diginorm_max($data_href);
		_check_con_algo($data_href);
		_check_keep_tmp_files($data_href);
		
		# Check optional FLASH parameters
		_check_min_overlap($data_href->{flash_params});
		_check_max_overlap($data_href->{flash_params});
		_check_ratio($data_href->{flash_params});
		_check_phred_offset($data_href->{flash_params});
		_check_read_len($data_href->{flash_params});
		_check_frag_len($data_href->{flash_params});
		_check_sd($data_href->{flash_params});
        
        # Check QC parameters
		_check_trim_to_base($data_href->{qc_params});
		_check_min_len($data_href->{qc_params});
		_check_min_avg_qual($data_href->{qc_params});
		_check_allow_gaps($data_href->{qc_params});
		_check_allow_ambig($data_href->{qc_params});
		_check_min_c_score($data_href->{qc_params});
		
		# Check each sample for correct tags
		_check_samples($data_href);
		
		return $data_href;
	}
	
	sub _check_samples {
		my ($data_href) = @_;
		
		if ( ref($data_href->{sample}) eq "ARRAY" ) {
			
			if ( @{$data_href->{sample}} == 0 ) {
				# there are no sample entries
				MTToolbox::MyX::MTParams::MissingTagValue->throw(
					error => "Missing tag value in params href",
					tag_name => 'sample',
				);
			}
			
			foreach my $sample_entry_href ( @{ $data_href->{sample} } ) {
				check_sample_entry($sample_entry_href,
								   $data_href->{seq_type},
								   $data_href->{split_samples_root_dir},
				);
			}
		}
		else {
			# there is only one sample meaning it is stored in the config
			# file NOT as an array ref but as a hash entry
			check_sample_entry($data_href->{sample},
							   $data_href->{seq_type},
							   $data_href->{split_samples_root_dir},
			);
			
			# I change the $data_href->{sample} object to point to an array ref
			# because that is what remainin parts of the program rely on
			my @arr = ($data_href->{sample});
			$data_href->{sample} = \@arr;
		}
		
		return 1;
	}
    
    sub check_sample_entry {
		my ($sample_href, $seq_type_g, $split_samples_root_dir_g) = @_;
		
		#print Dumper($sample_href);

		# Make sure all required tags are present
		foreach my $tag ( keys %required_sample_tags ) {
            _is_required_tag_defined($tag, $sample_href);
		}

		# Check for tags that are unknown
		foreach my $tag ( keys %{$sample_href} ) {
			if ( ! defined $required_sample_tags{$tag} and
				 ! defined $optional_sample_tags{$tag}
				) {
				MTToolbox::MyX::MTParams::UnrecognizedTag->throw(
                    error => "Unrecognized Tag",
                    tag_name => $tag,
                );
			}
		}
		
		# Check for unknown FLASH paramters
		if ( defined $sample_href->{flash_params} ) {
			foreach my $tag ( keys %{$sample_href->{flash_params}} ) {
				if ( ! defined $optional_flash_params_tags{$tag} ) {
					MTToolbox::MyX::MTParams::UnrecognizedTag->throw(
                        error => "Unrecognized FLASH Tag",
                        tag_name => $tag,
                    );
				}
			}
		}
        
        # Check local sample variables
        _check_fwd_primer($sample_href) if defined $sample_href->{fwd_primer};
		_check_rev_primer($sample_href) if defined $sample_href->{rev_primer};
		_check_seq_type($sample_href) if defined $sample_href->{seq_type};
		_check_fwd_linker($sample_href) if defined $sample_href->{fwd_linker};
		_check_rev_linker($sample_href) if defined $sample_href->{rev_linker};
		_check_merge_params_file($sample_href) if defined $sample_href->{merge_params_file};
		_check_fwd_mer_len($sample_href) if defined $sample_href->{fwd_mer_len};
		_check_rev_mer_len($sample_href) if defined $sample_href->{rev_mer_len};
		_check_fwd_max_shifts($sample_href) if defined $sample_href->{fwd_max_shifts};
		_check_rev_max_shifts($sample_href) if defined $sample_href->{rev_max_shifts};
        _check_min_con_depth($sample_href) if defined $sample_href->{min_con_depth};
        _check_diginorm_max($sample_href) if defined $sample_href->{diginorm_max};
		_check_con_algo($sample_href) if defined $sample_href->{con_algo};
		
		# Check local sample FLASH parameters
        my $flash_href = $sample_href->{flash_params};
		_check_min_overlap($flash_href) if defined $flash_href->{m};
		_check_max_overlap($flash_href) if defined $flash_href->{M};
		_check_ratio($flash_href) if defined $flash_href->{x};
		_check_phred_offset($flash_href) if defined $flash_href->{p};
		_check_read_len($flash_href) if defined $flash_href->{r};
		_check_frag_len($flash_href) if defined $flash_href->{f};
		_check_sd($flash_href) if defined $flash_href->{s};
		
		# check local sample QC parameters
		my $qc_href = $sample_href->{qc_params};
		_check_trim_to_base($qc_href) if defined $qc_href->{con_trim_to_base};
		_check_min_len($qc_href) if defined $qc_href->{con_min_len};
		_check_min_avg_qual($qc_href) if defined $qc_href->{con_min_avg_qual};
		_check_allow_gaps($qc_href) if defined $qc_href->{con_allow_gaps};
		_check_allow_ambig($qc_href) if defined $qc_href->{con_allow_ambig};
		_check_min_c_score($qc_href) if defined $qc_href->{con_min_c_score};
		_check_trim_to_base($qc_href) if defined $qc_href->{SRC_trim_to_base};
		_check_min_len($qc_href) if defined $qc_href->{SRC_min_len};
		_check_min_avg_qual($qc_href) if defined $qc_href->{SRC_min_avg_qual};
		_check_allow_gaps($qc_href) if defined $qc_href->{SRC_allow_gaps};
		_check_allow_ambig($qc_href) if defined $qc_href->{SRC_allow_ambig};
		_check_min_c_score($qc_href) if defined $qc_href->{SRC_min_c_score};

		# Make sure seq_type and files agree
		{
			my $_seq_type = undef;  # a local version of seq_type
			
			# Set the local _seq_type by looking first in the sample and then in
			# the global href
			if ( _is_optional_tag_defined('seq_type', $sample_href) ) {
				_check_seq_type($sample_href);
				$_seq_type = $sample_href->{seq_type};
			}
			else {  # else use the global seq_type
				# I don't have to check_seq_type here because it was already done
				# on a global level in _parse_params_file
				$_seq_type = $seq_type_g;
			}
            
            # translate the _seq_type
            $_seq_type = _translate_seq_type($_seq_type);
			
			# SE
			if ( $_seq_type == 1 ) {
                _is_required_tag_defined('fwd_file', $sample_href);
			}
			
			# PE w/ Overlap
			if ( $_seq_type == 2 ) {
				# Do they provide a merged file
				if ( _is_optional_tag_defined('merged_file', $sample_href) ) {
					# Do nothing.  But this means that fwd and rev files are
                    # not required
				}
				# or fwd and rev files
				else {
                    _is_required_tag_defined('fwd_file', $sample_href);
                    _is_required_tag_defined('rev_file', $sample_href);
				}
			}
			
			# PE w/o Overlap
			if ( $_seq_type == 3 ) {
                _is_required_tag_defined('fwd_file', $sample_href);
                _is_required_tag_defined('rev_file', $sample_href);
			}    
		}
		
		# Check the given sequence files
		{
			# If they are defined we have already confirmed that we need it
            # based on seq_type
			
			# fwd
			if ( my $_fwd_file = $sample_href->{fwd_file} ) {
                _check_fastq_file("$split_samples_root_dir_g/$_fwd_file");
			}
			
			# rev
			if ( my $_rev_file = $sample_href->{rev_file} ) {
                _check_fastq_file("$split_samples_root_dir_g/$_rev_file");
			}
			
			# merged
			if ( my $_merged_file = $sample_href->{merged_file} ) {
                _check_fastq_file("$split_samples_root_dir_g/$_merged_file");
			}
		}
		
		return 1; # true
	}
    
    sub print_xml_file {
        my ($self, $out_file) = @_;
        
        open my $OUT, ">", $out_file or
            croak("Cannot open file: $out_file");
    
        my $xml = XMLout($self->get_params_href());
        print $OUT $xml;
        
        close($OUT);
        
        return 1;
    }
	
	sub _get_required_global_tags {
		return \%required_global_tags;
	}
	
	sub _get_optional_global_tags {
		return \%optional_global_tags;
	}
	
	sub _get_required_sample_tags {
		return \%required_sample_tags;
	}
	
	sub _get_optional_sample_tags {
		return \%optional_sample_tags;
	}
	
	sub _get_optional_flash_params_tags {
		return \%optional_flash_params_tags;
	}
    
    sub _get_optional_qc_params_tags {
        return \%optional_qc_params_tags;
    }
	
    sub _is_required_tag_present {
		my ($tag, $href) = @_;
		
		# Is the XML tag in the params href
		if ( ! defined $href->{$tag} ) {
            MTToolbox::MyX::MTParams::MissingRequiredTag->throw(
                error => "Missing tag in params href",
                tag_name => $tag,
            );
		}
		
		return 1;
	}
	
	sub _is_required_tag_defined {
		my ($tag, $href) = @_;
        		
		# A tag must be present to be defined
        _is_required_tag_present($tag, $href);
		
		
		# XML::simple returns an empty hash ref when there is no value in between a tag.
		# I check for the situtation where the tag is present but no value has been entered
		if ( ref($href->{$tag}) eq "HASH" and keys %{$href->{$tag}} == 0 ) {
            MTToolbox::MyX::MTParams::MissingTagValue->throw(
                error => "Missing tag value in params href",
                tag_name => $tag,
            );
		}
        
        # it is just white space
        if ( $href->{$tag} =~ m/^\s+$/ or $href->{$tag} eq '' ) {
            MTToolbox::MyX::MTParams::MissingTagValue->throw(
                error => "Missing tag value in params href",
                tag_name => $tag,
            );
		}
		
		return 1;
	}
	
	sub _is_optional_tag_present {
		my ($tag, $href) = @_;
		
		# Is the XML tag actually in the XML file
		if ( ! defined $href->{$tag} ) {
			return 0;
		}
		
		return 1;
	}
	
	sub _is_optional_tag_defined {
		my ($tag, $href) = @_;
		
		# A tag must be present to be defined
		if ( ! _is_optional_tag_present($tag, $href) ) {
			return 0;  # return false if not present
		}
		
		# XML::simple returns an empty hash ref when there is no value inbetween a tag.
		# I check for the situtation where the tag is present but no value has been entered
		if ( ref($href->{$tag}) eq "HASH" and keys %{$href->{$tag}} == 0 ) {
			return 0;  # return false if not defined
		}
		
		return 1;  # true
	}
    
    sub _check_bin_sample_exec {
		my ($data_href) = @_;
		
		if ( ! defined $data_href->{bin_sample_exec} ) {
			$data_href->{bin_sample_exec} = "bin_sample.pl";
			return 1;
			
			# NOTE: I don't need to do the check below if the user wants to use
			# the bin_sample.pl executable that comes installed with MTToolbox.
		}
        
        # I added this because sometimes bin_sample.pl is not executable as
        # tested in the if statement below because it is not in the current
        # directory.
        if ( $data_href->{bin_sample_exec} eq "bin_sample.pl" ) {
            return 1;
        }
		
		my $bin_sample_exec = $data_href->{bin_sample_exec};
		if ( ! -x $bin_sample_exec ) {
			MTToolbox::MyX::MTParams::MissingExec->throw(
                error => "Can't find executable file",
                exec_name => $bin_sample_exec,
            );
		}
		
		return 1; # true
	}
	
    sub _check_split_samples_root_dir {
		my ($data_href) = @_;
		
		my $split_samples_root_dir = $data_href->{split_samples_root_dir};
		
		if ( $split_samples_root_dir ne '' and
			! -d $split_samples_root_dir) {
            MyX::Generic::DoesNotExist::Dir->throw(
                error => "Directory does not exist",
                dir_name => $split_samples_root_dir,
            );
		}
		
		return 1; # true
	}
	
    sub _check_output_root {
		my ($data_href) = @_;
		
		my $output_root = $data_href->{output_root};
		
		if ( ! -d $output_root ) {
			mkdir $output_root;
			carp("Making output_root directory: $output_root");
		}
		
		return 1; # true
	}
	
    sub _check_fwd_linker {
		my ($data_href) = @_;
		# NOTE: the forward linker is an optional global parameter.
		# If it is not defined globally it should be defined in EACH sample entry
		
		my $fwd_linker = $data_href->{fwd_linker};
		
		# skip these tests if the forward linker is not defined
		if ( ! defined $fwd_linker or
			 ref($fwd_linker) eq "HASH" and ! keys %{$fwd_linker} ) {
			return "";
		}
		
		# Warn if the length is a little long
		if ( length $fwd_linker > 10 ) {
			MTToolbox::MyX::MTParams::LongLinker->throw(
                error => "FWD linker is longer than 10 bp",
                orientation => 'fwd',
                length => length $fwd_linker,
            );
		}
		
		# Make all characters in linker upper case
		$fwd_linker = uc $fwd_linker;
		
		# Check for bad characters in forward linker sequence
		my %good_chars = map { $_ => 1 } qw(A T C G);
		foreach my $char ( split //, $fwd_linker ) {
			if ( ! defined $good_chars{$char} ) {
                MTToolbox::MyX::MTParams::UnrecognizedChar->throw(
                    error => "Unrecognized Char in fwd linker",
                    char => $char,
                );
			}
		}
		
		# because I change the characters to upper case I am re-saving the
		# forward linker sequence just incase it changed to upper case
		$data_href->{fwd_linker} = $fwd_linker;
		
		return 1; # true
	}
	
	sub _check_rev_linker {
		my ($data_href) = @_;
		# NOTE: the reverse linker is an optional global parameter.
		# If it is not defined globally it should be defined in each sample entry
		
		my $rev_linker = $data_href->{rev_linker};
		
		# skip these tests if the reverse linker is not defined
		if ( ! defined $rev_linker or
			 ref($rev_linker) eq "HASH" and ! keys %{$rev_linker} ) {
			return "";
		}
		
		# Warn if the lenght is a little long
		if ( length $rev_linker > 10 ) {
			MTToolbox::MyX::MTParams::LongLinker->throw(
                error => "REV linker is longer than 10 bp",
                orientation => 'rev',
                length => length $rev_linker,
            );
		}
		
		# Make all characters in reverse linker upper case
		$rev_linker = uc $rev_linker;
		
		# Check for bad characters in reverse linker sequence
		my %good_chars = map { $_ => 1 } qw(A T C G);
		foreach my $char ( split //, $rev_linker ) {
			if ( ! defined $good_chars{$char} ) {
				MTToolbox::MyX::MTParams::UnrecognizedChar->throw(
                    error => "Unrecognized Char in rev linker",
                    char => $char,
                );
			}
		}
		
		# because I change the characters to upper case I am re-saving the
		# reverse linker sequence just incase it changed to upper case
		$data_href->{rev_linker} = $rev_linker;
		
		return 1; # true
	}
	
    sub _check_fwd_primer {
		my ($data_href) = @_;
		
		return 1; # true
	}
	
	sub _check_rev_primer {
		my ($data_href) = @_;
		
		return 1; # true
	}
	
    sub _check_fwd_mer_len {
		my ($data_href) = @_;
		
		# An optional global value
		
		my $_fwd_mer_len = $data_href->{fwd_mer_len};
		
		# set to DEFAULT if not defined
		if ( ! defined $_fwd_mer_len or
			 ( ref($_fwd_mer_len) eq "HASH" and ! keys %{$_fwd_mer_len} ) ) {
			$data_href->{fwd_mer_len} = $FWD_MER_LEN;
			return 1;
		}
		
		if ( ! looks_like_number($_fwd_mer_len) ) {
			MyX::Generic::Digit::MustBeDigit->throw(
                error => "fwd mer len must be a digit",
                value => $_fwd_mer_len,
            )
		}
		
		if ( $_fwd_mer_len < 0 ) {
			MyX::Generic::Digit::TooSmall->throw(
                error => "fwd mer len must be > 0",
                value => $_fwd_mer_len,
                MIN => 0,
            );
		}
		
		return 1; # true
	}

	sub _check_rev_mer_len {
		my ($data_href) = @_;
		
		# An optional global value
		
		my $_rev_mer_len = $data_href->{rev_mer_len};
		
		# set to DEFAULT if not defined
		if ( ! defined $_rev_mer_len or
			 (ref($_rev_mer_len) eq "HASH" and ! keys %{$_rev_mer_len}) ) {
			$data_href->{rev_mer_len} = $REV_MER_LEN;
			return 1;
		}
		
		if ( ! looks_like_number($_rev_mer_len) ) {
			MyX::Generic::Digit::MustBeDigit->throw(
                error => "rev mer len must be a digit",
                value => $_rev_mer_len,
            )
		}
		
		if ( $_rev_mer_len < 0 ) {
			MyX::Generic::Digit::TooSmall->throw(
                error => "rev mer len must be > 0",
                value => $_rev_mer_len,
                MIN => 0,
            );
		}
		
		return 1; # true
	}
	
	sub _check_fwd_max_shifts {
		my ($data_href) = @_;
		
		# An optional global value
		
		my $_fwd_max_shifts = $data_href->{fwd_max_shifts};
		
		# set to DEFAULT if not defined
		if ( ! defined $_fwd_max_shifts or
			 (ref($_fwd_max_shifts) eq "HASH" and ! keys %{$_fwd_max_shifts}) ) {
			$data_href->{fwd_max_shifts} = $FWD_MAX_SHIFTS;
			return 1;
		}
		
		if ( ! looks_like_number($_fwd_max_shifts) ) {
			MyX::Generic::Digit::MustBeDigit->throw(
                error => "fwd max shifts must be a digit",
                value => $_fwd_max_shifts,
            )
		}
		
		if ( $_fwd_max_shifts < 0 ) {
			MyX::Generic::Digit::TooSmall->throw(
                error => "fwd max shifts must be > 0",
                value => $_fwd_max_shifts,
                MIN => 0,
            );
		}
		
		return 1; # true
	}
	
	sub _check_rev_max_shifts {
		my ($data_href) = @_;
		
		# An optional global value
		
		my $_rev_max_shifts = $data_href->{rev_max_shifts};
		
		# set to DEFAULT if not defined
		if ( ! defined $_rev_max_shifts or
			 (ref($_rev_max_shifts) eq "HASH" and ! keys %{$_rev_max_shifts}) ) {
			$data_href->{rev_max_shifts} = $REV_MAX_SHIFTS;
			return 1;
		}
		
		if ( ! looks_like_number($_rev_max_shifts) ) {
			MyX::Generic::Digit::MustBeDigit->throw(
                error => "rev max shifts must be a digit",
                value => $_rev_max_shifts,
            )
		}
		
		if ( $_rev_max_shifts < 0 ) {
			MyX::Generic::Digit::TooSmall->throw(
                error => "rev max shifts must be > 0",
                value => $_rev_max_shifts,
                MIN => 0,
            );
		}
		
		return 1; # true
	}
    
	sub _check_seq_type {
		my ($data_href) = @_;
		
		# NOTE: No DEFUALT is set here because this should be carefully
		# considered by the user.
		
		my $seq_type = $data_href->{seq_type};
        
        # translate if needed
        if ( $seq_type =~ m/SE/ ){ $seq_type = 1;}
        elsif ( $seq_type =~ m/PE w\/ Overlap/ ) {$seq_type = 2;}
        elsif ($seq_type =~ m/PE w\/o Overlap/ ) { $seq_type = 3;}
		
		if ( ! looks_like_number($seq_type) ) {
            my $message = "seq_type must be either 1, 2, or 3\n" .
                          "1 => Single end reads\n" .
                          "2 => Paired end overlapping reads\n" .
                          "3 => Paired end non-overlapping reads\n";
            
			MTToolbox::MyX::MTParams::BadSeqType->throw(
                error => $message,
                value => $seq_type,
            );
		}
		
		if ( $seq_type < 1 or $seq_type > 3 ) {
			my $message = "seq_type must be either 1, 2, or 3\n" .
                          "1 => Single end reads\n" .
                          "2 => Paired end overlapping reads\n" .
                          "3 => Paired end non-overlapping reads\n";
            
			MTToolbox::MyX::MTParams::BadSeqType->throw(
                error => $message,
                value => $seq_type,
            );
		}
		
		return 1; # true
	}

	sub _check_parallelize{
		my ($data_href) = @_;
		
		my $parallelize = $data_href->{parallelize};
		
		# right now this is only a boolean.  That might change as I get more
		# sofisticated.
		
		# set to DEFAULT if not defined
		if ( ! defined $parallelize or
			 (ref($parallelize) eq "HASH" and ! keys %{$parallelize}) ) {
			$data_href->{parallelize} = $PARALLELIZE;
			return 1;
		}
		
		# Make sure parallelize is upper case
		$parallelize = uc $parallelize;
		
		# Check for bad characters
		my %good_yes_values = map { $_ => 1 } qw(Y YES y yes);
		my %good_no_values = map { $_ => 1 } qw(N NO n no);
		if (     ! defined $good_yes_values{$parallelize}
			 and ! defined $good_no_values {$parallelize}
             and $parallelize ne 0
             and $parallelize ne 1
			) {
            MyX::Generic::BadValue->throw(
                error => "Unrecognized value in parallelize",
                value => $parallelize,
            );
		}
		
		if ( defined $good_yes_values{$parallelize} or
             $parallelize eq 1 ) {
			if ( ! `which bsub` ) {
                MTToolbox::MyX::MTParams::MissingExec->throw(
                    error => "Cannot find bsub system command",
                    exec_name => 'bsub',
                );
			}
			#$data_href->{parallelize} = 1;  # resets parallelize to a boolean value
			return 1;
		}
		
		#$data_href->{parallelize} = 0;  # resets parallelize to a boolean value
		return 1;  # true
	}
	
	sub _check_keep_tmp_files {
		my ($data_href) = @_;
		
		my $keep_tmp_files = $data_href->{keep_tmp_files};
		
		# right now this is only a boolean.  That might change as I get more
		# sofisticated.
		
		# set to DEFAULT if not defined
		if ( ! defined $keep_tmp_files or
			 (ref($keep_tmp_files) eq "HASH" and ! keys %{$keep_tmp_files}) ) {
			$data_href->{keep_tmp_files} = $KEEP_TMP_FILES;
			return 1;
		}
		
		# Make sure keep_tmp_files is upper case
		$keep_tmp_files = uc $keep_tmp_files;
		
		# Check for bad characters
		my %good_yes_values = map { $_ => 1 } qw(Y YES y yes Yes);
		my %good_no_values = map { $_ => 1 } qw(N NO n no No);
		if (     ! defined $good_yes_values{$keep_tmp_files}
			 and ! defined $good_no_values {$keep_tmp_files}
             and $keep_tmp_files ne 0
             and $keep_tmp_files ne 1
			) {
            MyX::Generic::BadValue->throw(
                error => "Unrecognized value in keep_tmp_files",
                value => $keep_tmp_files,
            );
		}
		
		if ( defined $good_yes_values{$keep_tmp_files} or
             $keep_tmp_files eq 1 ) {
			# resets keep_tmp_files to a boolean value
			#$data_href->{keep_tmp_files} = 1;  
			return 1;
		}
		
		# resets keep_tmp_files to a boolean value
		#$data_href->{keep_tmp_files} = 0;  
		return 1;  # true
	}
    
    sub _check_min_con_depth {
        my ($data_href) = @_;
        
        my $_min_con_depth = $data_href->{min_con_depth};
		
		# set to DEFAULT if not defined
		if ( ! defined $_min_con_depth or
			 (ref($_min_con_depth) eq "HASH" and ! keys %{$_min_con_depth}) ) {
			$data_href->{min_con_depth} = $MIN_CON_DEPTH;
			return 1;
		}
        
        if ( ! looks_like_number($_min_con_depth) ) {
            MyX::Generic::Digit::MustBeDigit->throw(
                error => "min_con_depth must be an int > 1\n",
                value => $_min_con_depth,
            );
		}
		
		if ( $_min_con_depth < 2 ) {
			my $message = "min_con_depth must be > 1\n";
            MyX::Generic::Digit::TooSmall->throw(
                error => "min_con_depth must be > 1\n",
                value => $_min_con_depth,
                MIN => 2,
            );
		}
		
		return 1; # true
    }
    
    sub _check_diginorm_max {
        my ($data_href) = @_;
        
        my $_diginorm_max = $data_href->{diginorm_max};
		
		# set to DEFAULT if not defined
		if ( ! defined $_diginorm_max or
			 (ref($_diginorm_max) eq "HASH" and ! keys %{$_diginorm_max}) ) {
			$data_href->{diginorm_max} = $DIGINORM_MAX;
			return 1;
		}
        
		# NA means that NO digital normalization should be done
        if ( $_diginorm_max =~ m/NA/i ) {
            return 1; # true
        }
        
        if ( ! looks_like_number($_diginorm_max) ) {
            MyX::Generic::Digit::MustBeDigit->throw(
                error => "diginorm_max must be an int > 1 or NA\n",
                value => $_diginorm_max,
            );
		}
		
		if ( $_diginorm_max < 2 ) {
            MyX::Generic::Digit::TooSmall->throw(
                error => "diginorm_max must be an int > 1 or NA\n",
                value => $_diginorm_max,
                MIN => 2,
            );
		}
		
		return 1; # true
    }
	
	sub _check_con_algo {
		my ($data_href) = @_;
		
		my $_con_algo = $data_href->{con_algo};
		
		# set to DEFAULT if not defined
		if ( ! defined $_con_algo or
			 (ref($_con_algo) eq "HASH" and ! keys %{$_con_algo}) ) {
			$data_href->{con_algo} = $CON_ALGO;
			return 1;
		}
		
		# Check for bad characters
		my %good_values = ('Muscle' => 1,
						   'muscle' => 1,
						   'Clustalw' => 1,
						   'ClustalW' => 1,
						   'clustalw' => 1,
						   'NoMSA' => 1,
						   'nomsa' => 1,
						  );
		
		if ( ! defined $good_values{$_con_algo} ) {
            MyX::Generic::BadValue->throw(
                error => 'Unrecognized value in con_algo.  ' .
						 'Options: Muscle | Clustalw | NoMSA',
                value => $_con_algo,
            );
		}
		
		return 1;
	}
	
	sub _check_min_overlap {
		my ($data_href) = @_;
		
		if ( my $m = $data_href->{m} ) {
			if ( ! looks_like_number($m) ) {
                MyX::Generic::Digit::MustBeDigit->throw(
                    error => "FLASH param \'m\' is not a digit",
                    value => $m,
                );
			}
			if ( $m < 0 ) {
				MyX::Generic::Digit::TooSmall->throw(
                    error => "FLASH param \'m\' must be > 0",
                    value => $m,
                    MIN => 0,
                );
			}
		}
		
		return 1;
	}
	
	sub _check_max_overlap {
		my ($data_href) = @_;
		
		if ( my $M = $data_href->{M} ) {
			if ( ! looks_like_number($M) ) {
				MyX::Generic::Digit::MustBeDigit->throw(
                    error => "FLASH param \'M\' is not a digit",
                    value => $M,
                );
			}
			if ( $M < 0 ) {
				MyX::Generic::Digit::TooSmall->throw(
                    error => "FLASH param \'M\' must be > 0",
                    value => $M,
                    MIN => 0,
                );
			}
		}
		
		return 1;
	}
	
	sub _check_ratio {
		my ($data_href) = @_;
		
		if ( my $x = $data_href->{x} ) {
			if ( ! looks_like_number($x) ) {
				MyX::Generic::Digit::MustBeDigit->throw(
                    error => "FLASH param \'x\' is not a digit",
                    value => $x,
                );
			}
			if ( $x < 0 or $x > 1 ) {
				MyX::Generic::Digit::OOB->throw(
                    error => "FLASH param \'x\' must be > 0 and < 1",
                    value => $x,
                    MIN => 0,
                    MAX => 1,
                );
			}
		}
		
		return 1;
	}
	
	sub _check_phred_offset {
		my ($data_href) = @_;

		if ( my $p = $data_href->{p} ) {
			if ( ! looks_like_number($p) ) {
				MyX::Generic::Digit::MustBeDigit->throw(
                    error => "FLASH param \'p\' is not a digit",
                    value => $p,
                );
			}
			if ( $p < 0 ) {
				MyX::Generic::Digit::TooSmall->throw(
                    error => "FLASH param \'p\' must be > 0",
                    value => $p,
                    MIN => 0,
                );
			}
		}
		
		return 1;
	}
	
	sub _check_read_len {
		my ($data_href) = @_;

		if ( my $r = $data_href->{r} ) {
			if ( ! looks_like_number($r) ) {
				MyX::Generic::Digit::MustBeDigit->throw(
                    error => "FLASH param \'r\' is not a digit",
                    value => $r,
                );
			}
			if ( $r < 0 ) {
				MyX::Generic::Digit::TooSmall->throw(
                    error => "FLASH param \'r\' must be > 0",
                    value => $r,
                    MIN => 0,
                );
			}
		}
		
		return 1;
	}
	
	sub _check_frag_len {
		my ($data_href) = @_;

		if ( my $f = $data_href->{f} ) {
			if ( ! looks_like_number($f) ) {
				MyX::Generic::Digit::MustBeDigit->throw(
                    error => "FLASH param \'f\' is not a digit",
                    value => $f,
                );
			}
			if ( $f < 0 ) {
				MyX::Generic::Digit::TooSmall->throw(
                    error => "FLASH param \'f\' must be > 0",
                    value => $f,
                    MIN => 0,
                );
			}
		}
		
		return 1;
	}
	
	sub _check_sd {
		my ($data_href) = @_;

		if ( my $s = $data_href->{s} ) {
			if ( ! looks_like_number($s) ) {
				MyX::Generic::Digit::MustBeDigit->throw(
                    error => "FLASH param \'s\' is not a digit",
                    value => $s,
                );
			}
			if ( $s < 0 ) {
				MyX::Generic::Digit::TooSmall->throw(
                    error => "FLASH param \'s\' must be > 0",
                    value => $s,
                    MIN => 0,
                );
			}
		}
		
		return 1;
	}
	
	sub _check_trim_to_base {
		my ($data_href) = @_;
		
		# check consensus min length
		if ( my $last_base = $data_href->{con_trim_to_base} ) {
            if ( ! defined $last_base ) {
                $data_href->{con_trim_to_base} = 'NA';
            }
			elsif ( $last_base =~ m/NA/i ) {;} # do nothing if min_len is NA
			elsif ( ! looks_like_number($last_base) ) {
				MyX::Generic::Digit::MustBeDigit->throw(
                    error => "Con filter param trim_to_base is not a digit",
                    value => $last_base,
                );
			}
			elsif ( $last_base < 0 ) {
				MyX::Generic::Digit::TooSmall->throw(
                    error => "Con filter param trim_to_base must be > 0",
                    value => $last_base,
                    MIN => 0,
                );
            }
		}
		
		# check single read category (SRC) min len
		if ( my $last_base = $data_href->{SRC_trim_to_base} ) {
            if ( ! defined $last_base ) {
                $data_href->{SRC_trim_to_base} = 'NA';
            }
			elsif ( $last_base =~ m/NA/i ) {;} # do nothing if min_len is NA
			elsif ( ! looks_like_number($last_base) ) {
				MyX::Generic::Digit::MustBeDigit->throw(
                    error => "SRC filter param trim_to_base is not a digit",
                    value => $last_base,
                );
			}
			elsif ( $last_base < 0 ) {
				MyX::Generic::Digit::TooSmall->throw(
                    error => "SRC filter param trim_to_base must be > 0",
                    value => $last_base,
                    MIN => 0,
                );
			}
		}
		
		return 1;
	}
	
	sub _check_min_len {
		my ($data_href) = @_;
		
		# check consensus min length
		if ( my $min_len = $data_href->{con_min_len} ) {
            if ( ! defined $min_len ) {
                $data_href->{con_min_len} = 'NA';
            }
			elsif ( $min_len =~ m/NA/i ) {;} # do nothing if min_len is NA
			elsif ( ! looks_like_number($min_len) ) {
				MyX::Generic::Digit::MustBeDigit->throw(
                    error => "Con filter param min_len is not a digit",
                    value => $min_len,
                );
			}
			elsif ( $min_len < 0 ) {
				MyX::Generic::Digit::TooSmall->throw(
                    error => "Con filter param min_len must be > 0",
                    value => $min_len,
                    MIN => 0,
                );
            }
		}
		
		# check single read category (SRC) min len
		if ( my $min_len = $data_href->{SRC_min_len} ) {
            if ( ! defined $min_len ) {
                $data_href->{SRC_min_len} = 'NA';
            }
			elsif ( $min_len =~ m/NA/i ) {;} # do nothing if min_len is NA
			elsif ( ! looks_like_number($min_len) ) {
				MyX::Generic::Digit::MustBeDigit->throw(
                    error => "SRC filter param min_len is not a digit",
                    value => $min_len,
                );
			}
			elsif ( $min_len < 0 ) {
				MyX::Generic::Digit::TooSmall->throw(
                    error => "SRC filter param min_len must be > 0",
                    value => $min_len,
                    MIN => 0,
                );
			}
		}
		
		return 1;
	}
	
	sub _check_min_avg_qual {
		my ($data_href) = @_;
		
		# check consensus min avg qual
		if ( my $min_avg_qual = $data_href->{con_min_avg_qual} ) {
			if ( ! looks_like_number($min_avg_qual) ) {
				MyX::Generic::Digit::MustBeDigit->throw(
                    error => "Con filter param min_avg_qual is not a digit",
                    value => $min_avg_qual,
                );
			}
			if ( $min_avg_qual < 0 ) {
				MyX::Generic::Digit::TooSmall->throw(
                    error => "Con filter param min_avg_qual must be > 0",
                    value => $min_avg_qual,
                    MIN => 0,
                );
			}
		}
		
		# check single read category (SRC) min avg qual
		if ( my $min_avg_qual = $data_href->{SRC_min_avg_qual} ) {
			if ( ! looks_like_number($min_avg_qual) ) {
				MyX::Generic::Digit::MustBeDigit->throw(
                    error => "SRC filter param min_avg_qual is not a digit",
                    value => $min_avg_qual,
                );
			}
			if ( $min_avg_qual < 0 ) {
				MyX::Generic::Digit::TooSmall->throw(
                    error => "SRC filter param min_avg_qual must be > 0",
                    value => $min_avg_qual,
                    MIN => 0,
                );
			}
		}
		
		return 1;
	}
	
	sub _check_allow_gaps {
		my ($data_href) = @_;
		
		my %vals = map { $_ => 1 } qw(Y y yes N n No);
		
		# check consensus allow gaps
		if ( my $allow_gaps = $data_href->{con_allow_gaps} ) {
			if ( ! defined $vals{$allow_gaps} ) {
				MyX::Generic::BadValue->throw(
                    error => "Con allow_gaps must be Y or N",
                    value => $allow_gaps,
                );
			}
		}
		
		# check single read category (SRC) allow gaps
		if ( my $allow_gaps = $data_href->{SRC_allow_gaps} ) {
			if ( ! defined $vals{$allow_gaps} ) {
				MyX::Generic::BadValue->throw(
                    error => "SRC allow_gaps must be Y or N",
                    value => $allow_gaps,
                );
			}
		}
		
		return 1;
	}
	
	sub _check_allow_ambig {
		my ($data_href) = @_;
		
		my %vals = map { $_ => 1 } qw(Y y YES Yes yes N n NO No no);
		
		# check consensus allow ambiguous bases
		if ( my $allow_ambig = $data_href->{con_allow_ambig} ) {
			if ( ! defined $vals{$allow_ambig} ) {
				MyX::Generic::BadValue->throw(
                    error => "Con allow_ambig must be Y or N",
                    value => $allow_ambig,
                );
			}
		}
		
		# check single read category (SRC) allow ambiguous bases
		if ( my $allow_ambig = $data_href->{SRC_allow_ambig} ) {
			if ( ! defined $vals{$allow_ambig} ) {
				MyX::Generic::BadValue->throw(
                    error => "SRC allow_ambig must be Y or N",
                    value => $allow_ambig,
                );
			}
		}
		
		return 1;
	}
	
	sub _check_min_c_score {
		my ($data_href) = @_;
		
		# check consensus min c-score
		if ( my $min_c_score = $data_href->{con_min_c_score} ) {
			if ( ! looks_like_number($min_c_score) ) {
				MyX::Generic::Digit::MustBeDigit->throw(
                    error => "Con filter param min_c_score is not a digit",
                    value => $min_c_score,
                );
			}
			if ( $min_c_score < 0 ) {
				MyX::Generic::Digit::TooSmall->throw(
                    error => "Con filter param min_c_score must be > 0",
                    value => $min_c_score,
                    MIN => 0,
                );
			}
		}
		
		# there is no check for the single read category (SRC) min c-score
		
		return 1;
	}
	
	sub _check_merge_params_file {
		my ($data_href) = @_;
		
		my $file = $data_href->{merge_params_file};
		my $seq_type = $data_href->{seq_type};
		
		# Make sure we even need a merge params file
		if ( $seq_type != 2 ) {
			return 1;
		}
		
		# Skip these tests if the paramter is not used.
		# The merging paramters can also be passed in via the flash_params href
		if ( ! defined $file ) {
			return 1;
		}
		
		# Make sure the file exists
		if ( ! -e $file ) {
			MyX::Generic::DoesNotExist::File->throw(
                error => "Merge params file does not exist",
                file_name => $file,
            );
		}
		
		# Make sure the file has something in it
		if ( ! -s $file ) {
			MyX::Generic::File::Empty->throw(
                error => "Merge params file is empty",
                file_name => $file,
            );
		}
		
		return 1; # true
	}

	sub _check_fastq_file {
		my ($file) = @_;
		
		# Make sure the file exists
		if ( ! -e $file ) {
			MyX::Generic::DoesNotExist::File->throw(
                error => "Fastq file does not exist",
                file_name => $file,
            );
		}
		
		# Make sure the file has something in it
		if ( ! -s $file ) {
			MyX::Generic::File::Empty->throw(
                error => "Fastq file is empty",
                file_name => $file,
            );
		}
		
		# Make sure the file ends in an appropriate fastq file extension
		my %file_extensions = map { $_ => 1} qw(fastq fq);
		my ($ext) = $file =~ /\.([^.]+)$/;
		if ( ! defined $file_extensions{$ext} ) {
            MyX::Generic::File::BadExtension->throw(
                error => "Unrecognized FASTQ file extension.  Try fastq or fq",
                file_name => $file,
                ext => $ext,
            );
		}
		
		return 1; # true
	}
    
    sub _translate_seq_type {
        my ($seq_type) = @_;
        
        my $seq_type_code = $seq_type;
        if ( $seq_type =~ m/SE/ ){ $seq_type_code = 1;}
        elsif ( $seq_type =~ m/PE w\/ Overlap/ ) {$seq_type_code = 2;}
        elsif ($seq_type =~ m/PE w\/o Overlap/ ) { $seq_type_code = 3;}
        
        return $seq_type_code;
    }
}



1; # Magic true value required at end of module
__END__

=head1 NAME

MTToolbox::MTParams - A module for holding and checking the paramters used in
					  MTToolbox


=head1 VERSION

This document describes MTToolbox::MTParams version 4.1.2


=head1 SYNOPSIS

    use MTToolbox::MTParams;
	
    # Build the object -- there are three ways
	my $mt_params = MTToolbox::MTParams->new({href => $my_href});
    my $mt_params = MTToolbox::MTParams->new({xml_file => $my_xml_file});
    my $mt_params = MTToolbox::MTParams->new();  #initialize params individually
    
    # 2 ways to set the params href other than through the constructor
    $mt_params->set_params_href({href => $my_href});
    $mt_params->set_params_href({xml_file = $my_xml_file});
    
    # Get the params href
    my $params_href = $mt_params->get_params_href();
    
    # Check the paramters
    eval{ $mt_params->check_params };
    if ( my $e = Exception::Class->caught('[An error]') ) {
        # Handle the error
    }
    
    # print the params to an xml file
    $mt_params->print_xml_file($xml_output_file);
	

  
=head1 DESRCIPTION

MTToolbox::MTParams holds all the parameters required for MTToolbox in a hash
reference that can easily be printed and parse using XML format.
MTToolbox::MTParams also checks the parameters to ensure the user inputs are
valid.  Error are thrown in the even of invalid parameters.


=head1 CONFIGURATION AND ENVIRONMENT

    None


=head1 DEPENDENCIES

    Class::Std::Utils
    List::MoreUtils qw(any)
    Readonly
    Carp qw(carp croak)
    XML::Simple
    Scalar::Util qw(looks_like_number)
    Exception::Class
    MyX::Generic 1.0.0
    MTToolbox::MyX::MTParams 4.1.2
    version our $VERSION = qv('4.1.2')


=head1 INCOMPATIBILITIES

    None reported.

=head1 METHODS

	get_params_href
	set_params_href
	check_params
	check_sample_entry
	print_xml_file
	
	_get_required_global_tags
	_get_optional_global_tags
	_get_required_sample_tags
	_get_optional_sample_tags
	_get_optional_flash_params_tags
    _get_optional_qc_params_tags
	
    _is_required_tag_present
    _is_required_tag_defined
    _is_optional_tag_present
    _is_optional_tag_defined
	
	_check_samples
    _check_bin_sample_exec
    _check_split_samples_root_dir
    _check_output_root
    _check_fwd_linker
	_check_rev_linker
    _check_fwd_primer
    _check_rev_primer
    _check_fwd_mer_len
    _check_rev_mer_len
    _check_seq_type
    _check_parallelize
	_check_keep_tmp_files
	_check_min_con_depth
	_check_dignorm_max
	_check_merge_params_file
    _check_fastq_file
	
	_check_min_overlap
	_check_max_overlap
	_check_ratio
	_check_phred_offset
	_check_read_len
	_check_frag_len
	_check_sd
    
	_check_trim_to_base
    _check_min_len
	_check_min_avg_qual
	_check_allow_gaps
	_check_allow_ambig
	_check_min_c_score

    _translate_seq_type    


=head1 METHODS DESRCIPTION

=head2 new

	Title: new
	Usage: my $mt_params = MTToolbox::MTParams->new({
								xml_file => "params_file.xml"
							});
	Function: Create a new MTToolbox::MTParams object
	Returns: Reference to a MTToolbox::MTParams object
	Args: -xml_file => a XML formated file with all parameters
          -href => a hash ref with all parameters
	Throws: croak 'Constructor called on existing object instead of class'
	Comments: NA
	See Also: NA

=head2 get_params_href

	Title: get_params_href
	Usage: my $params_href = $mt_params->get_params_href()
	Function: Returns the hash ref of parameters
	Returns: Hash Ref
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA

=head2 set_params_href

	Title: set_params_file
	Usage: $toolbox->set_params_file({xml_file => "params_file.xml"})
	Function: Sets the parameters file
	Returns: 1 on completion
	Args: -xml_file => a XML formated file with all parameters
          -href => a hash ref with all parameters
	Throws: NA
	Comments: If no parameters are passed an empty href is saved
	See Also: NA
    
=head2 check_params

	Title: check_params
	Usage: $toolbox->check_params()
	Function: Checks all parameters in the saved parameters href.  Checks
              for unknown tag values, required parameters, logic amoung
              paramter options (e.g. seq_type and file values), directories,
              files, and executables.
	Returns: Hash ref
	Args: NA
	Throws: MTToolbox::MyX::MTParams::UnrecognizedTag
            MTToolbox::MyX::MTParams::MissingRequiredTag
            MTToolbox::MyX::MTParams::MissingTagValue
            MTToolbox::MyX::MTParams::MissingExec
            MTToolbox::MyX::MTParams::LongLinker
            MTToolbox::MyX::MTParams::UrecognizedChar
            MTToolbox::MyX::MTParams::BadSeqType
            MyX::Generic::Undef::Param
            MyX::Generic::Undef::Attribute
            MyX::Generic::DoesNotExist::File
            MyX::Generic::DoesNotExist::Dir
            MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
            MyX::Generic::Digit::TooBig
            MyX::Generic::Digit::OOB
            MyX::Generic::File::Empty
            MyX::Generic::File::BadExtension
            MyX::Generic::BadValue
	Comments: There are many calls to internal methods in check_params.  Those
              calls are not encapsulated in eval statements.  However, any
              exceptions thrown in the internal methods are automatically
              passed to the caller of check_params.  It is recommended that
              those execptions be caught and handled by the caller of
              check_params.
	See Also: NA

=head2 check_sample_entry

	Title: check_sample_entry
	Usage: check_sample_entry($global_href, $sample_href)
	Function: Checks the merge_params_file parameter
	Returns: 1 on completion
	Args: -global_href => a hash ref of global parameter tags and thier values
		  -sample_href => a hash ref of sample parameter tags and their values
	Throws: -Croaks if required tags are missing
			-Croaks if there are unknown tags
			-Croaks with missing tag errors if give seq_types and file types do
			 not match
			-Croaks if given file types are not good fastq files
	Comments: Some of the Croak errors come from internal method calls
			  (ie _check_fastq_file)
	See Also: NA
    
=head2 print_xml_file

	Title: print_xml_file
	Usage: $toolbox->print_xml_file($output_file)
	Function: Prints the params href as an XML file
	Returns: 1 on completion
	Args: -output_file => file path to print the xml output to
	Throws: croak("Cannot open file: $out_file");
	Comments: NA
	See Also: NA

=head2 _get_required_global_tags

	Title: _get_required_global_tags
	Usage: _get_required_global_tags()
	Function: Returns the required global tags
	Returns: hash ref
	Args: NA
	Throws: NA
	Comments: The required global tags are automatically set as Readonly objects
			  when a MTToolbox is created.
	See Also: NA

=head2 _get_optional_global_tags

	Title: _get_optional_global_tags
	Usage: _get_optional_global_tags()
	Function: Returns the optional global tags
	Returns: hash ref
	Args: NA
	Throws: NA
	Comments: The optional global tags are automatically set as Readonly objects
			  when a MTToolbox is created.
	See Also: NA

=head2 _get_required_sample_tags

	Title: _get_required_sample_tags
	Usage: _get_required_sample_tags()
	Function: Returns the required sample tags
	Returns: hash ref
	Args: NA
	Throws: NA
	Comments: The required sample tags are automatically set as Readonly objects
			  when a MTToolbox is created.
	See Also: NA

=head2 _get_optional_sample_tags

	Title: _get_optional_sample_tags
	Usage: _get_optional_sample_tags()
	Function: Returns the optional sample tags
	Returns: hash ref
	Args: NA
	Throws: NA
	Comments: The optional sample tags are automatically set as Readonly objects
			  when a MTToolbox is created.
	See Also: NA

=head2 _get_optional_flash_params_tags

	Title: _get_optional_flash_params_tags
	Usage: _get_optional_flash_params_tags()
	Function: Returns the optional flash_params tags
	Returns: hash ref
	Args: NA
	Throws: NA
	Comments: The optional flash_params tags are automatically set as Readonly
			  objects when a MTToolbox is created.
	See Also: NA
    
=head2 _get_optional_qc_params_tags

	Title: _get_optional_qc_params_tags
	Usage: _get_optional_qc_params_tags()
	Function: Returns the optional quality control tags
	Returns: hash ref
	Args: NA
	Throws: NA
	Comments: The optional qc tags are automatically set as Readonly
			  objects when a MTToolbox is created.
	See Also: NA

=head2 _is_required_tag_present

	Title: _is_required_tag_present
	Usage: _is_required_tag_present($tag, $href)
	Function: Checks the given href for a required tag
	Returns: 1 if present
	Args: -tag => the tag to look for in the href
		  -href => a hash reference of tags and thier values
	Throws: MTToolbox::MyX::MTParams::MissingRequiredTag
	Comments: NA
	See Also: NA

=head2 _is_required_tag_defined

	Title: _is_required_tag_defined
	Usage: _is_required_tag_defined($tag, $href)
	Function: Checks the given href for a required tag
	Returns: 1 if present
	Args: -tag => the tag to look for in the href
		  -href => a hash reference of tags and thier values
	Throws: MTToolbox::MyX::MTParams::MissingRequiredTag
            MTToolbox::MyX::MTParams::MissingTagValue
	Comments: This calls _is_required_tag_present.  A tag has to be present to
			  be defined.
	See Also: NA

=head2 _is_optional_tag_present

	Title: _is_optional_tag_present
	Usage: _is_optional_tag_present($tag, $href)
	Function: Checks the given href for a optional tag
	Returns: boolean
	Args: -tag => the tag to look for in the href
		  -href => a hash reference of tags and thier values
	Throws: NA
	Comments: NA
	See Also: NA

=head2 _is_optional_tag_defined

	Title: _is_optional_tag_defined
	Usage: _is_optional_tag_defined($tag, $href)
	Function: Checks the given href for a optional tag
	Returns: boolean
	Args: -tag => the tag to look for in the href
		  -href => a hash reference of tags and thier values
	Throws: NA
	Comments: This calls _is_optional_tag_present.  A tag has to be present to
			  be defined.
	See Also: NA

=head2 _check_samples

	Title: _check_samples
	Usage: _check_samples($href)
	Function: Checks the sample tag parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MTToolbox::MyX::MTParams::MissingTagValue
	Comments: XML::Simple will store multiple samples as an array of hash refs.
			  However, if there is only one sample it stores the sample as a
			  single hash ref.  _check_samples() handles both cases and converts
			  a single hash ref to an array with a single hash ref element to be
			  compatable with the first case (i.e. an array of hash ref
			  samples).  If no samples are provided a MissingTagValue is thrown.
	See Also: NA

=head2 _check_bin_sample_exec

	Title: _check_bin_sample_exec
	Usage: _check_bin_sample_exec($href)
	Function: Checks the bin_sample_exec parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MTToolbox::MyX::MTParams::MissingExec
	Comments: NA
	See Also: NA
	
=head2 _check_split_samples_root_dir

	Title: _check_split_samples_root_dir
	Usage: _check_split_samples_root_dir($href)
	Function: Checks the split_samples_root_dir parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MyX::Generic::DoesNotExist::Dir
	Comments: NA
	See Also: NA

=head2 _check_output_root

	Title: _check_output_root
	Usage: _check_output_root($href)
	Function: Checks the output_root parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: carp("Making output_root directory: $output_root");
	Comments: NA
	See Also: NA

=head2 _check_fwd_linker

	Title: _check_fwd_linker
	Usage: _check_fwd_linker($href)
	Function: Checks the forward linker parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MTToolbox::MyX::MTParams::LongLinker
			MTToolbox::MyX::MTParams::UnrecognizedChar
	Comments: Returns an empty string if fwd linker is not globally defined.
			  This may be the case if a fwd linker is completly omitted from the
			  expectedvsequence pattern.  Also, the if the fwd linker has lower
			  case characeters they are converted to upper case.
	See Also: NA

=head2 _check_rev_linker

	Title: _check_rev_linker
	Usage: _check_rev_linker($href)
	Function: Checks the reverse linker parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MTToolbox::MyX::MTParams::LongLinker
			MTToolbox::MyX::MTParams::UnrecognizedChar
	Comments: Returns an empty string if rev linker is not globally defined.
			  This may be the case if a rev linker is completly omitted from the
			  expected sequence pattern.  Also, the if the rev linker has lower
			  case characeters they are converted to upper case.
	See Also: NA

=head2 _check_fwd_primer

	Title: _check_fwd_primer
	Usage: _check_fwd_primer($href)
	Function: Unimplemented
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: NA
	Comments: NA
	See Also: NA

=head2 _check_rev_primer

	Title: _check_rev_primer
	Usage: _check_rev_primer($href)
	Function: Unimplemented
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: NA
	Comments: NA
	See Also: NA

=head2 _check_fwd_len

	Title: _check_fwd_mer_len
	Usage: _check_fwd_mer_len($href)
	Function: Checks the fwd_mer_len parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
	Comments: DEFAULT = 8
	See Also: NA

=head2 _check_rev_mer_len

	Title: _check_rev_mer_len
	Usage: _check_rev_mer_len($href)
	Function: Checks the rev_mer_len parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
	Comments: DEFAULT = 5
	See Also: NA
	
=head2 _check_fwd_max_shift

	Title: _check_fwd_max_shift
	Usage: _check_fwd_max_shift($href)
	Function: Checks the fwd_max_shift parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: -MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
	Comments: DEFAULT = 5
	See Also: NA
	
=head2 _check_rev_max_shift

	Title: _check_rev_max_shift
	Usage: _check_rev_max_shift($href)
	Function: Checks the rev_max_shift parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
	Comments: DEFAULT = 5
	See Also: NA

=head2 _check_seq_type

	Title: _check_seq_type
	Usage: _check_seq_type($href)
	Function: Checks the seq_type parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MTToolbox::MyX::MTParams::BadSeqType
	Comments: 1 => Single end reads
			  2 => Paired end overlapping reads
			  3 => Paired end non-overlapping reads
	See Also: NA

=head2 _check_parallelize

	Title: _check_parallelize
	Usage: _check_parallelize($href)
	Function: Checks the parallelize parameter AND resets it to a boolean value
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MyX::Generic::BadValue
            MTToolbox::MyX::MTParams::MissingExec
	Comments: The parallelize parameter is reset in this method to be a boolean
	See Also: NA

=head2 _check_keep_tmp_files

	Title: _check_keep_tmp_files
	Usage: _check_keep_tmp_files($href)
	Function: Checks the keep_tmp_files parameter AND resets it to a boolean
			  value
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MyX::Generic::BadValue
	Comments: The keep_tmp_files parameter is reset in this method to be a
			  boolean
	See Also: NA

=head2 _check_min_con_depth

	Title: _check_min_con_depth
	Usage: _check_min_con_depth($href)
	Function: Checks the min_con_depth parameter value
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MyX::Generic::Digit::MustBeDigit
			MyX::Generic::Digit::TooSmall
	Comments: The min legal value is 2
	See Also: NA

=head2 _check_dignorm_max

	Title: _check_dignorm_max
	Usage: _check_dignorm_max($href)
	Function: Checks the diginorm_max parameter value
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MyX::Generic::Digit::MustBeDigit
			MyX::Generic::Digit::TooSmall
	Comments: The min legal value is 2
	See Also: NA

=head2 _check_merge_params_file

	Title: _check_merge_params_file
	Usage: _check_merge_params_file($href)
	Function: Checks the merge_params_file parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MyX::Generic::DoesNotExist::File
            MyX::Generic::File::Empty
	Comments: The merge_params_file was original used to pass the merging
			  parameters to flash.  There is now a hash object in the data_href
			  called flash_params which contains the the parameters for merging.
			  However, the merge_params_file argument can still be used but is
			  optional.  It is recommended for running MTToolbox to use the GUI
			  which create the hash object for flash_params.  When using the
			  commandline you can use the merge_params_file argument.
	See Also: NA

=head2 _check_fastq_file

	Title: _check_fastq_file
	Usage: _check_fastq_file($file)
	Function: Checks that file is a valid fastq file
	Returns: 1 on completion
	Args: -file => a file path to the fastq file to check
	Throws: -Croaks if file is empty
			-Croaks if there is an unrecognized extension (not fastq or fq)
	Comments: NA
	See Also: NA

=head2 _check_min_overlap

	Title: _check_min_overlap
	Usage: _check_min_overlap($href)
	Function: Checks the min_overlap parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
	Comments: A FLASH parameter
	See Also: NA

=head2 _check_max_overlap

	Title: _check_max_overlap
	Usage: _check_max_overlap($href)
	Function: Checks the max_overlap parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
	Comments: A FLASH parameter
	See Also: NA

=head2 _check_ratio

	Title: _check_ratio
	Usage: _check_ratio($href)
	Function: Checks the ratio parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::OOB
	Comments: A FLASH parameter
	See Also: NA

=head2 _check_phred_offset

	Title: _check_phred_offset
	Usage: _check_phred_offset($href)
	Function: Checks the phred_offset parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
	Comments: A FLASH parameter
	See Also: NA

=head2 _check_read_len

	Title: _check_read_len
	Usage: _check_read_len($href)
	Function: Checks the read_len parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
	Comments: A FLASH parameter
	See Also: NA

=head2 _check_frag_len

	Title: _check_frag_len
	Usage: _check_frag_len($href)
	Function: Checks the frag_len parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
	Comments: A FLASH parameter
	See Also: NA

=head2 _check_sd

	Title: _check_sd
	Usage: _check_sd($href)
	Function: Checks the sd parameter
	Returns: 1 on completion
	Args: -href => a hash reference of tags and thier values
	Throws: MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
	Comments: A FLASH parameter
	See Also: NA

=head2 _check_trim_to_base

	Title: _check_trim_to_base
	Usage: _check_trim_to_base($href)
	Function: Checks the trim to base position filter parameter
	Returns: 1 on completion
	Args: -href => a hash reference of filter parameters
	Throws: MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
	Comments: A Filtering parameter
	See Also: NA
    
=head2 _check_min_len

	Title: _check_min_len
	Usage: _check_min_len($href)
	Function: Checks the min length filter parameter
	Returns: 1 on completion
	Args: -href => a hash reference of filter parameters
	Throws: MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
	Comments: A Filtering parameter
	See Also: NA
    
=head2 _check_min_avg_qual

	Title: _check_min_avg_qual
	Usage: _check_min_avg_qual($href)
	Function: Checks the min average quality value filter parameter
	Returns: 1 on completion
	Args: -href => a hash reference of filter parameters
	Throws: MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
	Comments: A Filtering parameter
	See Also: NA
    
=head2 _check_allow_gaps

	Title: _check_allow_gaps
	Usage: _check_allow_gaps($href)
	Function: Checks the allow gaps filter parameter
	Returns: 1 on completion
	Args: -href => a hash reference of filter parameters
	Throws: MyX::Generic::BadValue
	Comments: A Filtering parameter
	See Also: NA
    
=head2 _check_allow_ambig

	Title: _check_allow_ambig
	Usage: _check_allow_ambig($href)
	Function: Checks the allow ambiguous bases filter parameter
	Returns: 1 on completion
	Args: -href => a hash reference of filter parameters
	Throws: MyX::Generic::BadValue
	Comments: A Filtering parameter
	See Also: NA
	
=head2 _check_min_c_score

	Title: _check_min_c_score
	Usage: _check_min_c_score($href)
	Function: Checks the min c-score value filter parameter
	Returns: 1 on completion
	Args: -href => a hash reference of filter parameters
	Throws: MyX::Generic::Digit::MustBeDigit
            MyX::Generic::Digit::TooSmall
	Comments: A Filtering parameter
	See Also: NA
    
=head2 _translate_seq_type

	Title: _translate_seq_type
	Usage: _translate_seq_type($seq_type)
	Function: Translates string versions of seq_type to number code
	Returns: Int
	Args: -seq_type => either the string or number code
	Throws: NA
	Comments: If the number code is passed in as the parameters it is not
              changed and returned again.  The code is as follows:
              1 => SE
              2 => PE w/ Overlap
              3 => PE w/o Overlap
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

