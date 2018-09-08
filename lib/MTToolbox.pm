package MTToolbox;

use warnings;
use strict;

use version; our $VERSION = qv('4.1.2');
use Class::Std::Utils;
use List::MoreUtils qw(any);
use Readonly;
use Carp qw(carp croak);
use File::Temp qw(tempfile);
use XML::Simple;
use Data::Dumper; # For debuging XML config files
use Statistics::Descriptive;
use Scalar::Util qw(looks_like_number);
use Benchmark;
use Text::Wrap;      # <hc /> for writing readme files

use BioUtils::FastqIO 1.0.7;
use BioUtils::FastaIO 1.0.7;
use BioUtils::SeqSet::Diagnostics 1.0.7;
use BioUtils::QC::FastqFilter 1.0.7;
use MyX::Generic 0.0.3;
use MTToolbox::MTDepthDist 4.1.2;
use MTToolbox::RunSummary 4.1.2;
use MTToolbox::MTParams 4.1.2;
use MTToolbox::MTEnv 4.1.2;


{
	# Usage statement
	Readonly my $NEW_USAGE => q{ new( {
										params_file => ,
										} ) };
	
    # Attributes #
	my %params_file_of;
    
    # Subroutines #
	sub run;
	sub get_params_file;
	sub set_params_file;
	
    sub _get_command;
    sub _get_serial_command;
    sub _get_parallel_command;
	
	sub _stall;
	sub _cat_outputs;
	sub _consolidate_summaries;
	sub _gaps_and_ambig_bases_analysis;
	sub _print_mt_depth_hist;
	
	sub _clean_up;
	sub _to_bool;
    
    
    ###############
    # Constructor #
    ###############
    sub new {
        my ($class, $arg_href) = @_;
        
        # Croak if calling 'new' on an already blessed reference
        croak 'Constructor called on existing object instead of class'
            if ref $class;
		
		# Make sure the required parameters are defined
        if ( any {!defined $_} $arg_href->{params_file},
            ) {
            MyX::Generic::Undef::Param->throw(
                                              error => 'Undefined parameter value',
                                              usage => $NEW_USAGE,
                                              );
        }
        
        # Bless a scalar to instantiate an object
        my $new_obj = bless \do{my $anon_scalar}, $class;
		
		# attributes
		$params_file_of{ident $new_obj} = $arg_href->{params_file};
        
        return $new_obj;
    }
    
    ###############
	# Subroutines #
	###############
	sub run {
		my ($self) = @_;
		
		# start the Benchmark
		my $t0 = Benchmark->new();
		
		# Save the params file for future analysis
		my $params_file = $params_file_of{ident $self};
		
		# Parameters file checks
		warn "### Check Parameters\n\n";
		{
			# Make sure the params files is defined
			if ( ! defined $params_file or $params_file eq "" ) {
				croak("Please provide a parameters file\n$NEW_USAGE");
			}

			# Make sure the params file exists
			if ( ! -e $params_file ) {
				croak("Parameters file, $params_file, doesn't exist\n");
			}
			
			# Make sure the params file has something in it
			if ( ! -s $params_file ) {
				croak("Parameters file, $params_file, is empty\n");
			}
		}
		
		# Parse the params file and check its inputs
		my $params_obj = MTToolbox::MTParams->new({xml_file => $params_file});
		$params_obj->check_params();
		my $data_href = $params_obj->get_params_href();
		
		# check the environment
		my $mt_env = MTToolbox::MTEnv->new();
		$mt_env->check_env();
		
		# save output_root for repeated use
		my $output_root = $data_href->{output_root};
		
		# Run the commands
		foreach my $sample_entry_href ( @{ $data_href->{sample} } ) {
			my $command = _get_command($data_href, $sample_entry_href);
			
			warn "### Running Command --\n$command\n\n";
			system("$command");
		}
		
		# Wait for commands to finish if using a parallel system
		_stall($data_href->{parallelize}, $output_root);
		
		
		# Concatenate individual sample outputs
		warn "### Cat outputs\n\n";
		_cat_outputs($output_root);
		
		
		# Consolidate individual sample summary files
		warn "### Consolidate Summaries\n\n";
		_consolidate_summaries($output_root);
		
		
		# Analyze gaps and ambiguous bases
		# I broke this when I moved the QC into the sample folders
		#warn "### Analyze Gaps and Ambiguous Bases\n\n";
		#_gaps_and_ambig_bases_analysis($output_root);
		
		
		# Print seqs_per_mt_hist files (i.e. a txt file and ps file)
		warn "### Print MT Depth Histogram\n\n";
		_print_mt_depth_hist($output_root);
		
		
		# Some clean up stuff like deleting temp files
		warn "### Clean up\n\n";
		_clean_up($output_root,
				  $data_href->{keep_tmp_files},
				  $data_href->{parallelize}
				  );
		

        _file_master($output_root);

		# FINISHED
		warn "### FINISHED\n\n";
		
		# finishing benchmark
		my $t1 = Benchmark->new();
		my $td = timediff($t1, $t0);
		warn "Time: ", timestr($td), "\n";
		
		return 1;
		
	}
	
	sub get_params_file {
		my ($self) = @_;
		return $params_file_of{ident $self};
	}
	
	sub set_params_file {
		my ($self, $file) = @_;
		$params_file_of{ident $self} = $file;
		return 1;
	}

    sub _get_command {
		my ($global_href, $sample_entry_href) = @_;
		
		my $command = undef;

		if ( _to_bool( $global_href->{parallelize} ) ) {
			$command = _get_parallel_command($global_href, $sample_entry_href);
		}
		else {
			$command = _get_serial_command($global_href, $sample_entry_href);
		}
		
		return $command;
	}
	
    sub _get_serial_command {
		my ($global_href, $sample_entry_href) = @_;
		
		# Some local variables to help build the command
		my $_seq_type = undef;
		
		# Add the executable
		my $command = "$global_href->{bin_sample_exec} ";
		
		# Add the unique sample ID
		$command .= "--id $sample_entry_href->{sample_id} ";
		
		# Add the seq_type
		if ( defined $sample_entry_href->{seq_type} ) {
			$command .= "--seq_type $sample_entry_href->{seq_type} ";
			$_seq_type = $sample_entry_href->{seq_type}; # a local variable for later use
		}
		else {
			$command .= "--seq_type $global_href->{seq_type} ";
			$_seq_type = $global_href->{seq_type};  # a local variable for later use
		}
		
		# Add the barcode -- There will always be a barcode
		$command .= "--barcode $sample_entry_href->{barcode} ";
		
		# Add the merged file -- There may NOT always be a merged file
		my $split_samples_root_dir = $global_href->{split_samples_root_dir};
		if ( $_seq_type == 2 and defined $sample_entry_href->{merged_file} ) {
			$command .= "--merged_file $split_samples_root_dir/$sample_entry_href->{merged_file} ";
		}
		
		# Add forward file -- There may NOT always be a fwd file
		if (  $_seq_type == 1 or
			( $_seq_type == 2 and ! defined $sample_entry_href->{merged_file}) or 
			  $_seq_type == 3
			) {
			$command .= "--fwd_file $split_samples_root_dir/$sample_entry_href->{fwd_file} ";
		}
		
		# Add reverse file -- There may NOT always be a rev file based on the seq_type
		if ( ($_seq_type == 2 and ! defined $sample_entry_href->{merged_file}) or
			  $_seq_type == 3
			) {
			$command .= "--rev_file $split_samples_root_dir/$sample_entry_href->{rev_file} ";
		}
		
		# Add the forward linker
		if ( defined $sample_entry_href->{fwd_linker} ) {
			$command .= "--fwd_linker $sample_entry_href->{fwd_linker} ";
		}
		else {
			$command .= "--fwd_linker $global_href->{fwd_linker} ";
		}
		
		# Add the reverse linker
		if ( defined $sample_entry_href->{rev_linker} ) {
			$command .= "--rev_linker $sample_entry_href->{rev_linker} ";
		}
		else {
			$command .= "--rev_linker $global_href->{rev_linker} ";
		}
		
		# Add the forward primer
		if ( defined $sample_entry_href->{fwd_primer} ) {
			$command .= "--fwd_primer \"$sample_entry_href->{fwd_primer}\" ";
		}
		else {
			$command .= "--fwd_primer \"$global_href->{fwd_primer}\" ";
		}
		
		# Add the reverse primer
		if ( defined $sample_entry_href->{rev_primer} ) {
			$command .= "--rev_primer \"$sample_entry_href->{rev_primer}\" ";
		}
		else {
			$command .= "--rev_primer \"$global_href->{rev_primer}\" ";
		}
		
		# Add the forward random mer length
		if ( defined $sample_entry_href->{fwd_mer_len} ) {
			$command .= "--fwd_mer_len $sample_entry_href->{fwd_mer_len} ";
		}
		else {
			$command .= "--fwd_mer_len $global_href->{fwd_mer_len} ";
		}
		
		# Add the reverse random mer length
		if ( defined $sample_entry_href->{rev_mer_len} ) {
			$command .= "--rev_mer_len $sample_entry_href->{rev_mer_len} ";
		}
		elsif( defined $global_href->{rev_mer_len} ) {
			$command .= "--rev_mer_len $global_href->{rev_mer_len} ";
		}
		
		# Add the forward max shifts
		if ( defined $sample_entry_href->{fwd_max_shifts} ) {
			$command .= "--fwd_max_shifts $sample_entry_href->{fwd_max_shifts} ";
		}
		elsif ( defined $global_href->{fwd_max_shifts} ) {
			$command .= "--fwd_max_shifts $global_href->{fwd_max_shifts} ";
		}
		
		# Add the reverse max shifts
		if ( defined $sample_entry_href->{rev_max_shifts} ) {
			$command .= "--rev_max_shifts $sample_entry_href->{rev_max_shifts} ";
		}
		elsif ( defined $global_href->{rev_max_shifts} ) {
			$command .= "--rev_max_shifts $global_href->{rev_max_shifts} ";
		}
		
		# Add the min con depth
		if ( defined $sample_entry_href->{min_con_depth} ) {
			$command .= "--min_con_depth $sample_entry_href->{min_con_depth} ";
		}
		elsif ( defined $global_href->{min_con_depth} ) {
			$command .= "--min_con_depth $global_href->{min_con_depth} ";
		}
		
		# Add the diginorm max
		if ( defined $sample_entry_href->{diginorm_max} ) {
			$command .= "--diginorm_max $sample_entry_href->{diginorm_max} ";
		}
		elsif ( defined $global_href->{diginorm_max} ) {
			$command .= "--diginorm_max $global_href->{diginorm_max} ";
		}
		
		# Add the consensus algorithm
		if ( defined $sample_entry_href->{con_algo} ) {
			$command .= "--con_algo $sample_entry_href->{con_algo} ";
		}
		elsif ( defined $global_href->{con_algo} ) {
			$command .= "--con_algo $global_href->{con_algo} ";
		}
		
		# Add the output dir
		my $output_root = $global_href->{output_root};
		if ( defined $sample_entry_href->{output_dir} ) {
			$command .= "--output_dir $output_root/$sample_entry_href->{output_dir} ";
		}
		else {  # DEFAULT to id_barcode (ie P15_ATGTGA)
			$command .= "--output_dir $output_root/$sample_entry_href->{sample_id}_$sample_entry_href->{barcode} ";
		}
		
		# Add the keep_tmp_files value
		if ( defined $global_href->{keep_tmp_files} ) {
			$command .= "--keep_tmp_files $global_href->{keep_tmp_files} ";
		}
		
		# Add the merge parameters file
		# (only needed if the fwd and rev reads need to be merged)
		if ( my $merge_params_file = $global_href->{merge_params_file} ) {
			$command .= "--merge_params_file $merge_params_file ";
		}
		
		# Add the local flash parameters
		my %flash_params_added = ();
		if ( defined $sample_entry_href->{flash_params} ) {
			foreach my $tag ( keys %{$sample_entry_href->{'flash_params'}} ) {
				
				# a special case to deal with Getopt::Long being case-insensitive
				if ( $tag =~ m/M/ ) {
					$command .= "--flash_Max " .
								$sample_entry_href->{'flash_params'}->{$tag} .
								" ";
								next;
				}
				
				# all normal cases (every parameter but M)
				$command .= "--flash_" . $tag .
							" " .
							$sample_entry_href->{'flash_params'}->{$tag} .
							" ";
				$flash_params_added{$tag} = 1;
			}
		}
		
		# Add the global flash parameters
		if ( defined $global_href->{flash_params} ) {
			foreach my $tag ( keys %{$global_href->{'flash_params'}} ) {
				
				# a special case to deal with Getopt::Long being case-insensitive
				if ( $tag =~ m/M/ ) {
					$command .= "--flash_Max " .
								$global_href->{'flash_params'}->{$tag} .
								" ";
								next;
				}
				
				# all normal cases (every paramter but M)
				if ( ! defined $flash_params_added{$tag} ) {
					$command .= "--flash_" . $tag .
							" " .
							$global_href->{'flash_params'}->{$tag} .
							" ";
				}
			}
		}
		
		# Add the local QC parameters
		my %qc_params_added = ();
		if ( defined $sample_entry_href->{qc_params} ) {
			foreach my $tag ( keys %{$sample_entry_href->{'qc_params'}} ) {
				$command .= "--qc_" . $tag .
							" " .
							$sample_entry_href->{'qc_params'}->{$tag} .
							" ";
				$qc_params_added{$tag} = 1;
			}
		}
		
		# Add the global QC parameters
		if ( defined $global_href->{qc_params} ) {
			foreach my $tag ( keys %{$global_href->{'qc_params'}} ) {
				if ( ! defined $qc_params_added{$tag} ) {
					$command .= "--qc_" . $tag .
							" " .
							$global_href->{'qc_params'}->{$tag} .
							" ";
				}
			}
		}
		
		return $command;
	}
	
    sub _get_parallel_command {
		my ($global_href, $sample_entry_href) = @_;
		
		my $output_root = $global_href->{output_root};
		my $_output_dir = undef;  # a local variable for later use
		
		# Figure out the output directory and store in a local variable for later use
		if ( defined $sample_entry_href->{output_dir} ) {
			$_output_dir .= "$output_root/$sample_entry_href->{output_dir}";
		}
		else {  # DEFAULT to barcode name
			$_output_dir .= "$output_root/$sample_entry_href->{sample_id}_$sample_entry_href->{barcode}";
		}
		
		# Add the bsub portion of the command to the beginning of the serial portion of the command
		my $command = ("bsub -q week -n 1 " .
					   "-J MTBIN_$sample_entry_href->{sample_id} " .
					   "-o $_output_dir/log.out " .
					   "-e $_output_dir/log.err " .
					   _get_serial_command($global_href, $sample_entry_href)
					   );
		
		return $command;
	}
	
	sub _stall {
		my ($parallelize, $output_root) = @_;
		
		if ( _to_bool($parallelize) ) {
			my $command = 'bjobs > ' .  "$output_root/ACTIVE_JOBS" . '
						  while [ -s ' . "$output_root/ACTIVE_JOBS" . ' ]
						  do
							  sleep 10
							  bjobs -J MTBIN\* > ' . "$output_root/ACTIVE_JOBS" . '
						  done
						  rm ' . "$output_root/ACTIVE_JOBS";
			system($command);
		}
		
		return 1;
	}
	
	sub _cat_outputs {
		my ($output_root) = @_;
		# File to look for
		my @files_to_cat = (
			"consensusSeqs_HQ.fastq",
			"consensusSeqs_fwd_HQ.fastq",
			"consensusSeqs_rev_HQ.fastq",
			"consensusSeqs_LQ.fastq",
			"consensusSeqs_fwd_LQ.fastq",
			"consensusSeqs_rev_LQ.fastq",
			"categorizable_reads.fastq",
			"categorizable_reads_fwd.fastq",
			"categorizable_reads_rev.fastq",
			"single_read_categories_HQ.fastq",
			"single_read_categories_fwd_HQ.fastq",
			"single_read_categories_rev_HQ.fastq",
			"single_read_categories_LQ.fastq",
			"single_read_categories_fwd_LQ.fastq",
			"single_read_categories_rev_LQ.fastq",
			"seqs_per_mt_dist.txt",
		);
		
		foreach my $file ( @files_to_cat ) {
			my $com = "cat " . $output_root . "/*/" . $file . " > " .
				$output_root . "/all_" . $file;
			warn $com;
			`$com`;
			
			# remove the file all_* file if it is empty
			if ( -z "$output_root/all_" . $file ) {
				warn "Removing empty file: all_" . $file;
				system("rm $output_root/all_" . $file);
			}
		}
		
		# cat the diagnostics tables
		my @diag_files = ("consensusSeqs_diagnostics.txt",
						  "consensusSeqs_fwd_diagnostics.txt",
						  "single_read_categories_diagnostics.txt",
						  "single_read_categories_fwd_diagnostics.txt"
						  );
		foreach my $file ( @diag_files ) {
			my @con_diag_files = `find $output_root -name "$file"`;
			my $out = "$output_root/all_" . $file;
			my $first = 1;

			foreach my $found ( @con_diag_files ) {
				chomp $found;
				
				if ( $first ) {
					`cat $found > $out`;
					$first = 0;

				}
				else {
					`grep -v \"#\" $found >> $out`;
				}
			}
		}

		
		# a system command for percent merged file for all samples
		# NOTE: I commented this out because it is in the summary file already
		#       - Mar 1, 2013
		#system("grep \"Percent merged\" $output_root/P*/summary.txt | " .
		#	   "sed 's/:Percent merged: / /g' | " .
		#	   "sort -n -r -k 2 | " .
		#	   q{sed 's/[ |\/]/\t/g' | } .
		#	   q{awk '{OFS="\t"}{print $(NF-2), $NF}' } . 
		#	   "> $output_root/all_merged_percents.txt"
		#	   );
		#
		## remove the above file if it is empty
		#if ( -z "$output_root/all_merged_percents.txt" ) {
		#	#carp("Removing empty file: all_merged_percents.txt");
		#	system("rm $output_root/all_merged_percents.txt");
		#}
		
		return 1;
	}
		
	sub _consolidate_summaries {
		my ($output_root) = @_;
		
		my $run_summary = MTToolbox::RunSummary->new();
		
		$run_summary->load_sample_summaries($output_root);
		
		$run_summary->print_summary_files($output_root);
		
		return 1;
	}
	
	sub _gaps_and_ambig_bases_analysis() {
		my ($output_root) = @_;
		# -s == file exists and has non-zero size
		if ( -s "$output_root/all_consensusSeqs_HQ.fastq" ) { 
			my $diag_obj = BioUtils::SeqSet::Diagnostics->new();
			$diag_obj->run_diagnosis("$output_root/all_consensusSeqs_HQ.fastq");
			$diag_obj->graph_ambig_base_pos_dist("$output_root/con_ambig_base_pos_dist");
			$diag_obj->graph_ambig_base_comp_dist("$output_root/con_ambig_base_comp_dist");
			$diag_obj->graph_ambig_bases_per_seq("$output_root/con_ambig_bases_per_seq");
		}
		if ( -s "$output_root/all_consensusSeqs_fwd_HQ.fastq" ) {
			my $diag_obj_fwd = BioUtils::SeqSet::Diagnostics->new();
			$diag_obj_fwd->run_diagnosis("$output_root/all_consensusSeqs_fwd_HQ.fastq");
			$diag_obj_fwd->graph_ambig_base_pos_dist("$output_root/con_ambig_base_pos_dist_fwd");
			$diag_obj_fwd->graph_ambig_base_comp_dist("$output_root/con_ambig_base_comp_dist_fwd");
			$diag_obj_fwd->graph_ambig_bases_per_seq("$output_root/con_ambig_bases_per_seq_fwd");
		}
		if ( -s "$output_root/all_consensusSeqs_rev_HQ.fastq" ) {
			my $diag_obj_rev = BioUtils::SeqSet::Diagnostics->new();
			$diag_obj_rev->run_diagnosis("$output_root/all_consensusSeqs_rev_HQ.fastq");
			$diag_obj_rev->graph_ambig_base_pos_dist("$output_root/con_ambig_base_pos_dist_rev");
			$diag_obj_rev->graph_ambig_base_comp_dist("$output_root/con_ambig_base_comp_dist_rev");
			$diag_obj_rev->graph_ambig_bases_per_seq("$output_root/con_ambig_bases_per_seq_rev");
		}
		
		return 1;
	}
	
	sub _print_mt_depth_hist {
		my ($output_root) = @_;
		
		# prints both hist file and hist graph in ps format
		
		my $dist_file = "$output_root/all_seqs_per_mt_dist.txt";
		my $hist_file = "$output_root/all_seqs_per_mt_hist.txt";
		my $hist_graph = "$output_root/all_seqs_per_mt_hist.ps";
		
		my $mt_depth_dist = MTToolbox::MTDepthDist->new({dist_file => $dist_file});
		
		# NOTE: the dist is already built when _cat_outputs is called.
		
		$mt_depth_dist->build_hist();
		$mt_depth_dist->print_hist($hist_file);
		$mt_depth_dist->print_hist_graph($hist_graph, "All Seqs MT Depth");
		
		return 1;
	}

    # <hc> - arrange files in a clear way and write readme files
    sub _file_master {
        # moves the files output by base MT-Toolbox to their respective directories
        ( my $output_root ) = @_;

        # move summary files
        ( -d "$output_root/1_summary_info" ) or mkdir("$output_root/1_summary_info") or carp("Could not make directory: $output_root/1_summary_info");
        system("mv $output_root/all_seqs_per_* $output_root/all_summary* $output_root/1_summary_info") == 0 or carp("Could not move summary files into $output_root/1_summary_info");

        # move conSeq files
        ( -d "$output_root/2_consensus_seqs" ) or mkdir("$output_root/2_consensus_seqs") or carp("Could not make directory: $output_root/2_consensus_seqs");
        system("mv $output_root/all_c* $output_root/2_consensus_seqs") == 0 or carp("Could not move consensus sequence files to: $output_root/2_consensus_seqs");
    
        # move single reads
        ( -d "$output_root/3_single_read_categories" ) or mkdir("$output_root/3_single_read_categories") or carp("Could not make directory: $output_root/3_single_read_categories");
        system("mv $output_root/all_single_read_categories_* $output_root/3_single_read_categories") == 0 or carp("Could not move single read categories files to: $output_root/3_single_read_caregories");


        #
        ## Begin writing READMEs
        #
    
        # hash to hold all the files with readme information.
        my %readme_hash = (    
            # summary files
            '1_summary_info' => "directory of summary files",
            'all_seqs_per_mt_dist.txt' => "text counts of number of sequences per molecule tag",
            'all_seqs_per_mt_hist.ps' => "histogram of number of MTs with each number of sequences",
            'all_seqs_per_mt_hist.txt' => "text representation of all_seqs_per_mt_hist.ps",
            'all_summary.ps' => "graph of sample summaries",
            'all_summary.txt' => "text representation of all_summary.ps",

            # conSeq files
            '2_consensus_seqs' => "directory of consensus sequence FASTA, FASTQ, and diagnostic files",
            'all_categorizable_reads.fastq' => "precursor to the consensusSeqs files",
            'all_consensus_and_singles_HQ.fasta' => "FASTA of all high-quality sequences (both consensus and singles)",
            'all_consensusSeqs_diagnostics.txt' => "table of information about the consensus sequences that were formed (quality, etc)",
            'all_consensusSeqs_HQ.fasta' => "FASTA of all high-quality consensus sequences",
            'all_consensusSeqs_HQ.fastq' => "FASTQ of all high-quality consensus sequences",
            'all_consensusSeqs_LQ.fastq', => "FASTQ of all low-quality consensus sequences",

            # Singles files
            '3_single_read_categories' => "directory of FASTA, FASTQ, and diagnostic files for reads that were not incorporated in consensus sequences",
            'all_single_read_categories_diagnostics.txt' => "table of information about the single reads (quality, etc)",
            'all_single_read_categories_HQ.fasta' => "FASTA of all high-quality single reads",
            'all_single_read_categories_HQ.fastq' => "FASTQ of all high-quality single reads",

            # sample files
            'P\d+_[GCAT]+' => "directory containing information about the sample identified with the specified barcode",

            # clean up files
            'log.out' => "stdout from removing temporary files",
            'log.err' => "stderr from removing temporary files",
        );

        _write_readmes($output_root, \%readme_hash);

           # write a special readme for the root 
        open(my $README, ">>", "$output_root/README");
        print $README "\n";
        print $README "Main Logs: (may be different if you ran MT-MTDriver manually)\n";
        print $README "*" x 60, "\n\n";
 
        my %root_hash = (
            'MT_MTDriver.err' => "STDERR from the main pipeline of MT-MT-Toolbox",
            'MT_MTDriver.log' => "STDOUT from the main pipeline of MT-MT-Toolbox",
        );

        # set the number of characters for the first column
        my $first_field_len = 50;

        local $Text::Wrap::columns = 110 - $first_field_len;
        foreach my $file ( sort keys %root_hash ) {
       
            printf($README "%-${first_field_len}s", $file);
            my $desc = Text::Wrap::wrap("","", $root_hash{$file});
            my @lines = split("\n", $desc);

            # write the first line of description
            printf($README "%-${Text::Wrap::columns}s\n", shift(@lines));
            # write any subsequent lines
            foreach ( @lines ) {
                printf($README "%-${first_field_len}s%-${Text::Wrap::columns}s\n", "", $_);
            }
            print $README "\n";

        }   # End writing of root readme
    }


    sub _write_readmes {

        my ( $output_root, $readme_hash ) = @_;
    
        my %to_write;
        opendir(my $DIR, $output_root) or die("Couldn't open directory $output_root\n");

        my @files = readdir($DIR) or die("Couldn't read from directory\n");

        # check if there were any files
        if (@files) {

            # add the info from each file to the hash (enter recursive mode if directory)
            foreach my $file ( @files) {

                next if $file =~ m/^\./;        # skip hidden files or folders
                foreach my $key ( keys %{$readme_hash} ) {
                    $to_write{$file} = $readme_hash->{$key} if ( $file =~ m/$key/ );
                }
                if (-d "$output_root/$file") {
                    _write_readmes("$output_root/$file", $readme_hash);
                }
            }
            # write readme file using my pretty printing method --still could use some work to ensure error free
            my $README;
            if ( -f "$output_root/README" ) {
                if ( keys %to_write ) {

                    open($README, ">>", "$output_root/README");
                    print $README "\n";
                    print $README "Additional associated files:\n";
                    print $README "*" x 60, "\n\n";
                }
                else {  # return if there are no files
                    return 
                }
            }
            else {
                open($README, ">", "$output_root/README");

                print $README "This is the README for $output_root\n";
                print $README scalar localtime(), "\n";
                print $README "*" x 60, "\n\n";

                # output message if no files
                if ( ! keys %to_write ) {
                    print $README "This directory was not used in your MT-MT-Toolbox run.";
                    close $README;
                    return;
                }
            }

            # set column widths for pretty printing
            # get the max file field may want to consider an option to turncate filename in readme if over a certain number of characters.
            my $max = 0;

            $to_write{"README"} = "this file";
            foreach ( keys %to_write ) { $max = length($_) if length($_) > $max; }

            # set the number of characters for the first column
            my $first_field_len = $max + 4;
            $first_field_len = 75 if $first_field_len > 75;     # weak attempt preventing long (really long) filenames from messing up the format, needs to be better

            local $Text::Wrap::columns = 110 - $first_field_len;
            foreach my $file ( sort keys %to_write ) {
                printf($README "%-${first_field_len}s", $file);

                my $desc = Text::Wrap::wrap("","", $to_write{$file});
                my @lines = split("\n", $desc);

                # write the first line of description
                printf($README "%-${Text::Wrap::columns}s\n", shift(@lines));
                # write any subsequent lines
                foreach ( @lines ) {
                    printf($README "%-${first_field_len}s%-${Text::Wrap::columns}s\n", "", $_);
                }
                print $README "\n";

            }   # End pretty print

            close($README)
        }   # end if(files)
    }


	sub _clean_up {
		my ($output_root, $keep_tmp_files, $parallelize) = @_;
		
		if ( ! _to_bool($keep_tmp_files) ) {
			if ( _to_bool($parallelize) ) {
				opendir( my $DIR, $output_root ) or
					carp("cannot open $output_root for deleting");
				
				while ( my $file = readdir($DIR) ) {
					if ( $file =~ m/(P\d+)_/ ) {
						my $command = ("bsub -q week -n 1 " .
						   "-J DEL_$1 " .
						   "-o $output_root/log.out " .
						   "-e $output_root/log.err " .
						   "rm -rf $output_root/$file"
						   );
						print "command: $command\n";
						system($command);
					}
				}
			}
			else {
				system("rm -rf $output_root/P*");
			}
			
			# these are single files that can be delited without parallelism
			# I don't need to delete these anymore because only the HQ and LQ
			# 	files are concatonated outside the sample directories
			#system("find $output_root -name all_single_read_categories.fastq -delete");
			#system("find $output_root -name all_consensusSeqs.fastq -delete");
			#system("find $output_root -name \'*fwd.fastq\' -delete");
			#system("find $output_root -name \'*rev.fastq\' -delete");
		}
		
		return 1;
	}
	
	sub _to_bool {
		my ($parallelize) = @_;
		
		if ( $parallelize eq 1 or $parallelize eq 0 ) {
			return $parallelize;
		}
		
		my %good_yes_values = map { $_ => 1 } qw(Y YES Yes y yes);
		if ( defined $good_yes_values{$parallelize} ) {
			return 1;
		}
		
		# else -- meaning no parallel
		return 0;
	}
}



1; # Magic true value required at end of module
__END__

=head1 NAME

MTToolbox - A module for running many types of molecule tagged sequence analysis


=head1 VERSION

This document describes MTToolbox version 4.1.2


=head1 SYNOPSIS

    use MTToolbox;
	
	my $mt_toolbox = MTToolbox->new({params_file => "params_file.xml"});
	
	# Getting and setting the parameters file
	$mt_toolbox->set_params_file("params_file.xml");
	my $params_file = $mt_toolbox->get_parmas_file();
	
	# Run the analysis
	$mt_toolbox->run();

  
=head1 DESRCIPTION

MTToolbox is a module for running many types of molecule tagged sequence
analysis.  Limitting the number of sequencing errors is crutial for precise
sequence comparisons.  Sequencing errors are often introduced during the PCR
steps and/or the sequencing step.  In an attempt to limit the number of errors
we have developed a molecular biology technique in which individual DNA
molecules are tagged with a unique barcode-like sequences.  We call this
sequence a molecule tag (MT).  After each DNA molecule has been tagged it can be
PCRed and sequenced.  After sequencing we recieve raw reads that contain a
molecule tag at the beginning.  MTToolbox.pl sorts those sequences based on their
molecule tags and builds a consensus sequence to represent the original DNA
molecule that was tagged.

=head1 PARAMETERS FILE

The only parameter required is a XML formated parameters file.  The parameters
file should look something like this:

	<?xml version="1.0" encoding="UTF-8" ?>
	<MTbin_run>
		<bin_sample_exec>/Users/Scott/Projects_test/MTToolbox/trunk/bin/bin_sample.pl</bin_sample_exec>
		<split_samples_root_dir>/Users/Scott/temp/temp_split_by_sample</split_samples_root_dir>
		<output_root>/Users/Scott/temp/temp_output_20120808</output_root>
		<parallelize>N</parallelize>
		<fwd_linker>CAGT</fwd_linker>
		<rev_linker>CAGT</rev_linker>
		<fwd_primer>GTGCCAGC[AC]GCCGCGGTAA</fwd_primer>
		<rev_primer>GGACTAC[ACT][ACG]GGGT[AT]TCTAAT</rev_primer>
		<fwd_primer_len>9</fwd_primer_len>
		<rev_primer_len>4</rev_primer_len>
		<seq_type>1</seq_type>
		<flash_params>
			<m>10</m>
			<M>70</M>
			<x>.20</x>
			<p>33</p>
			<r>250</r>
			<f>310</f>
			<s>20</s>
		</flash_params>
		<sample>
			<barcode>ATCACG</barcode>
			<fwd_file>A01_ATCACG_L001_R1_001.fastq</fwd_file>
			<output_dir></output_dir>
		</sample>
		<sample>
			<barcode>GATCAG</barcode>
			<fwd_file>A02_GATCAG_L001_R1_001.fastq</fwd_file>
			<output_dir></output_dir>
		</sample>
	</MTbin_run>

=head1 OUTPUTS
    
    The main output file is called all_consensusSeqs.fastq and contains the
    consensus sequences for each molecule tag bin in each sample.  This file
    can be found amoung other intermediate files and directories in the
    output_root directory specified in the parameters file.  The first value in
    the header of the sequence is the sample ID as defined in the configuration
    XML file.  The second value in the header is the molecule tag sequence for
    that bin.
    
    A secondary output file is called summary.txt.  This file simply contains
    information about each sample with some summary statistics for each sample
    and a combination of all samples.  A file called seqs_per_mt_dist.txt
    contains the distribution of the number of reads found in each molecule tag
    bin.  A file called categorizable_reads.fastq is a fastq file of all the
    reads that match the expected pattern and can therefore be classified in
    catigories by their molecule tag.  A file called single_read_categories.fastq
    is a fastq formated file that contains all the reads that are categorizable
    but only have one read that fits in that category.
    
    There are also several intermediate files and directories that are saved
    to facilitate downstream analysis.  Each sample has a directory that
    contains information about the sequences originating from that sample.
    There are several output files in these directories.  For a complete
    description of these files see the documentation for bin_sample.pl.


=head1 CONFIGURATION AND ENVIRONMENT

- gnuplot must be installed and locatable by Chart::Graph::Gnuplot qw(gnuplot)
- flash must be installed and executable on the commandline to merge paired-end
  sequences with overlapping ends


=head1 DEPENDENCIES

	version our $VERSION = qv('4.1.2')
	Class::Std::Utils
	List::MoreUtils qw(any)
	Readonly
	Carp qw(carp croak)
	File::Temp qw(tempfile)
	XML::Simple
	Data::Dumper # For debuging XML config files
	Statistics::Descriptive
	Scalar::Util qw(looks_like_number)

	BioUtils::FastqIO 1.0.0
	BioUtils::FastaIO 1.0.0
	BioUtils::SeqSet::Diagnostics 1.0.0
	BioUtils::QC::FastqFilter 1.0.0
	BioUtils::QC::ContaminantFilter 1.0.0
	MyX::Generic 1.0.0
	MTToolbox::MTDepthDist 4.1.2
	MTToolbox::RunSummary 4.1.2
	MTToolbox::MTParams 4.1.2


=head1 INCOMPATIBILITIES

	None reported.

=head1 METHODS

	run
	get_params_file
	set_params_file
    
    _get_command
    _get_serial_command
    _get_parallel_command
	
	_stall
	_cat_outputs
	_consolidate_summaries
	_gaps_and_ambig_bases_analysis
	_print_mt_depth_hist
	
	_clean_up
	__to_bool


=head1 METHODS DESRCIPTION

=head2 new

	Title: new
	Usage: my $toolbox = MTToolbox->new({params_file => "params_file.xml"});
	Function: Create a new MTToolbox object
	Returns: Reference to a MTToolbox object
	Args: -params_file => a XML formated file with all parameters
	Throws: MyX::Generic::Undef::Param
	Comments: NA
	See Also: NA

=head2 run

	Title: run
	Usage: $toolbox->run()
	Function: Runs a full analysis on molecule tagged sequences
	Returns: 1 on completion
	Args: NA
	Throws: Croaks on parameter file problems
	Comments: NA
	See Also: NA

=head2 get_params_file

	Title: get_params_file
	Usage: $toolbox->get_params_file()
	Function: Returns the file name of the parameters file
	Returns: String
	Args: NA
	Throws: NA
	Comments: NA
	See Also: NA

=head2 set_params_file

	Title: set_params_file
	Usage: $toolbox->set_params_file("params.xml")
	Function: Sets the parameters file
	Returns: 1 on completion
	Args: -file => a XML formated file with all parameters
	Throws: NA
	Comments: NA
	See Also: NA

=head2 _get_command

	Title: _get_command
	Usage: _get_command($global_href, $sample_href)
	Function: Gets the bin_sample.pl command for the given sample
	Returns: Sting
	Args: -global_href => a hash ref of global parameter tags and thier values
		  -sample_href => a hash ref of sample parameter tags and their values
	Throws: NA
	Comments: When this is called all the parameters should have been checked
	See Also: NA

=head2 _get_serial_command

	Title: _get_serial_command
	Usage: _get_serial_command($global_href, $sample_href)
	Function: Gets the bin_sample.pl serial command for the given sample
	Returns: Sting
	Args: -global_href => a hash ref of global parameter tags and thier values
		  -sample_href => a hash ref of sample parameter tags and their values
	Throws: NA
	Comments: When this is called all the parameters should have been checked
	See Also: NA

=head2 _get_parallel_command

	Title: _get_parallel_command
	Usage: _get_parallel_command($global_href, $sample_href)
	Function: Gets the bin_sample.pl parallel command for the given sample
	Returns: Sting
	Args: -global_href => a hash ref of global parameter tags and thier values
		  -sample_href => a hash ref of sample parameter tags and their values
	Throws: NA
	Comments: When this is called all the parameters should have been checked
	See Also: NA

=head2 _stall

	Title: _stall
	Usage: _stall($parallelize, $output_root)
	Function: Stalls a parallel command waiting for all the samples to finish
	Returns: 1 on completion
	Args: -parallelize => a boolean indicating if a run is parallel or serial
		  -output_root => where to put the ACTIVE_JOBS file
	Throws: NA
	Comments: Will only work on a LSF system with a bjobs command
	See Also: NA

=head2 _cat_outputs

	Title: _cat_outputs
	Usage: _cat_outputs($output_root)
	Function: concatenates some output files from individual sample output
			  directories into one larger output file in the output_root
			  directory
	Returns: 1 on completion
	Args: -output_root => the output root directory
	Throws: Carp if a file is not found
	Comments: Generally the carps from this method can be ignored because they
			  correspond to files that are not associated with this run type.
	See Also: NA

=head2 _consolidate_summaries

	Title: _consolidate_summaries
	Usage: _consolidate_summaries($output_root)
	Function: combines all the summary information for each of the samples into
			  a large summary file
	Returns: 1 on completion
	Args: -output_root => the output root directory
	Throws: NA
	Comments: The output file is called "all_summary.txt"
	See Also: NA

=head2 _gaps_and_ambig_bases_analysis

	Title: _gaps_and_ambig_bases_analysis
	Usage: _gaps_and_ambig_bases_analysis($output_root)
	Function: runs an analysis on gaps and ambiguous bases in the consensus seqs
	Returns: 1 on completion
	Args: -output_root => the output root directory
	Throws: NA
	Comments: NA
	See Also: NA
	
=head2 _print_mt_depth_hist

	Title: _print_mt_depth_hist
	Usage: _print_mt_depth_hist($output_root)
	Function: print seqs per mt histogram file and graph
	Returns: 1 on completion
	Args: -output_root => the output root directory
	Throws: NA
	Comments: Print both txt file and postscript file into output_root
	See Also: NA
	
=head2 _clean_up

	Title: _clean_up
	Usage: _clean_up($output_root, $keep_tmp_files, $parallelize)
	Function: Removes unneed and temp files
	Returns: 1 on completion
	Args: -output_root => the output root directory
		  -keep_tmp_files => boolean to indicate if temp files should be kept
		  -parallelize => boolean to indicate if jobs should be parallized
	Throws: NA
	Comments: NA
	See Also: NA
	
=head2 _file master

        Title: _file_master
        Usage: _file_master($output_root)
        Function: Moves files to a more organized location
        Returns: NA
        Args: -output_root => the base output path
        Throws: NA
        Comments: NA
        See Also: NA

=head2 _write_readmes

        Title: _write_readmes
        Usage: _write_readmes($output_root, $readme_hash)
        Function: Writes a README file in each directory
        Returns: NA
        Args: -output_root => the base output path
          -readme_hash => a hash of files to write READMEs for
        Throws: NA
        Comments: NA
        See Also: NA

=head2 __to_bool

	Title: __to_bool
	Usage: __to_bool($parallelize)
	Function: Translates No, N, Yes, or Y into boolean
	Returns: boolean
	Args: -parallelize => variable to translate to boolean value
	Throws: NA
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
