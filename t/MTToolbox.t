
use strict;
use warnings;

use MTToolbox;
use Test::More tests => 39;
use Test::Exception;
use Test::Warn;
use File::Temp qw(tempfile tempdir);
use File::Basename;

BEGIN { use_ok( 'MTToolbox'); }

# Create a temporary params file
my ($fh, $filename) = tempfile();

# Create a SampleSummary object to test
my $mt_toolbox;
lives_ok( sub { $mt_toolbox = MTToolbox->new( {params_file => $filename} ) },
        "calling new -- expected to live");
isa_ok($mt_toolbox, 'MTToolbox', "is MTToolbox");

# test run
{
    # TODO
    ;
}

# test get_params_file
{
    is( $mt_toolbox->get_params_file(), $filename, "get_params_file");
}

# test set_param_file
{
    my ($fh2, $filename2) = tempfile();
    lives_ok( sub{ $mt_toolbox->set_params_file($filename2) }, "setting new params file");
    is( $mt_toolbox->get_params_file(), $filename2, "getting new params file");
}

# test _get_command
{
    my $global_href = ();
    my $sample_href = ();
    
    # all globally assigned when possible
    $global_href->{bin_sample_exec} = "my_exec";
    $global_href->{seq_type} = 1;
    $global_href->{split_samples_root_dir} = "samples_dir";
    $global_href->{output_root} = "output_root";
    $global_href->{fwd_linker} = "ATCG";
	$global_href->{rev_linker} = "GCTA";
    $global_href->{fwd_primer} = "GGGG";
    $global_href->{rev_primer} = "AAAA";
    $global_href->{fwd_mer_len} = 9;
    $global_href->{parallelize} = 0;  # This is normally set as "N" and changed in _check_parallelize
	$global_href->{min_con_depth} = 2;
	$global_href->{diginorm_max} = "NA";
	$global_href->{con_algo} = 'MSA';
    $sample_href->{sample_id} = "P1";
    $sample_href->{barcode} = "AAAAA";
    $sample_href->{fwd_file} = "fwd_file";
    $sample_href->{output_dir} = "output_dir";
    
    # testing _get_command for a serial command
    is( MTToolbox::_get_command($global_href, $sample_href),
       "my_exec " .
       "--id P1 " .
       "--seq_type 1 " .
       "--barcode AAAAA " .
       "--fwd_file samples_dir/fwd_file " .
       "--fwd_linker ATCG " .
	   "--rev_linker GCTA " .
       "--fwd_primer \"GGGG\" " .
       "--rev_primer \"AAAA\" " .
       "--fwd_mer_len 9 " .
	   "--min_con_depth 2 " .
	   "--diginorm_max NA " .
	   "--con_algo MSA " .
       "--output_dir output_root/output_dir ",
       "_get_command() -- a serial one" );
    
    # testing _get_command for a parallel command
    $global_href->{parallelize} = 1;  # This is normally set as "Y" and changed in _check_parallelize
    my $expected_command = "bsub -q week -n 1 " .
                           "-J MTBIN_P1 " .
                           "-o output_root/output_dir/log.out " .
                           "-e output_root/output_dir/log.err " .
                           MTToolbox::_get_serial_command($global_href, $sample_href);
    is( MTToolbox::_get_command($global_href, $sample_href),
        $expected_command,
        "_get_command() -- a parallel one" );
}

# test _get_serial_command
{
    my $global_href = ();
    my $sample_href = ();
    
    my $command_index = 0;
    
    # all globally assigned when possible
    $global_href->{bin_sample_exec} = "my_exec";
    $global_href->{seq_type} = 1;
    $global_href->{split_samples_root_dir} = "samples_dir";
    $global_href->{output_root} = "output_root";
    $global_href->{fwd_linker} = "ATCG";
	$global_href->{rev_linker} = "GCTA";
    $global_href->{fwd_primer} = "GGGG";
    $global_href->{rev_primer} = "AAAA";
    $global_href->{fwd_mer_len} = 9;
	$global_href->{min_con_depth} = 2;
	$global_href->{diginorm_max} = "NA";
	$global_href->{con_algo} = 'MSA';
    $sample_href->{sample_id} = "P1";
    $sample_href->{barcode} = "AAAAA";
    $sample_href->{fwd_file} = "fwd_file";
    $sample_href->{output_dir} = "output_dir";
    is( MTToolbox::_get_serial_command($global_href, $sample_href),
       "my_exec " .
       "--id P1 " .
       "--seq_type 1 " .
       "--barcode AAAAA " .
       "--fwd_file samples_dir/fwd_file " .
       "--fwd_linker ATCG " .
	   "--rev_linker GCTA " .
       "--fwd_primer \"GGGG\" " .
       "--rev_primer \"AAAA\" " .
       "--fwd_mer_len 9 " .
	   "--min_con_depth 2 " .
	   "--diginorm_max NA " .
	   "--con_algo MSA " .
       "--output_dir output_root/output_dir ",
       "_get_serial_command - $command_index" );
    $command_index++;
    
    # with an overridden seq_type -- this is actually a wrong command but
    # effectivly test overridding seq_type in the sample entry
    $sample_href->{seq_type} = 2;
    $sample_href->{rev_file} = "rev_file";
    is( MTToolbox::_get_serial_command($global_href, $sample_href),
       "my_exec " .
       "--id P1 " .
       "--seq_type 2 " .
       "--barcode AAAAA " .
       "--fwd_file samples_dir/fwd_file " .
       "--rev_file samples_dir/rev_file " .
       "--fwd_linker ATCG " .
	   "--rev_linker GCTA " .
       "--fwd_primer \"GGGG\" " .
       "--rev_primer \"AAAA\" " .
       "--fwd_mer_len 9 " .
	   "--min_con_depth 2 " .
	   "--diginorm_max NA " .
	   "--con_algo MSA " .
       "--output_dir output_root/output_dir ",
       "_get_serial_command - $command_index" );
    $command_index++;
    
    # overrridden seq_type (2) with a merged file - no fwd or rev files appear
    $sample_href->{merged_file} = "merged_file";
    is( MTToolbox::_get_serial_command($global_href, $sample_href),
       "my_exec " .
       "--id P1 " .
       "--seq_type 2 " .
       "--barcode AAAAA " .
       "--merged_file samples_dir/merged_file " .
       "--fwd_linker ATCG " .
	   "--rev_linker GCTA " .
       "--fwd_primer \"GGGG\" " .
       "--rev_primer \"AAAA\" " .
       "--fwd_mer_len 9 " .
	   "--min_con_depth 2 " .
	   "--diginorm_max NA " .
	   "--con_algo MSA " .
       "--output_dir output_root/output_dir ",
       "_get_serial_command - $command_index" );
    $command_index++;
	
	# undo the sample_href override and default back to global
	delete $sample_href->{merged_file};
    delete $sample_href->{seq_type};
    
    # override the fwd_linker when seq_type == 1 as in the first few tests
    $sample_href->{fwd_linker} = "GCTA";
    is( MTToolbox::_get_serial_command($global_href, $sample_href),
       "my_exec " .
       "--id P1 " .
       "--seq_type 1 " .
       "--barcode AAAAA " .
       "--fwd_file samples_dir/fwd_file " .
       "--fwd_linker GCTA " .
	   "--rev_linker GCTA " .
       "--fwd_primer \"GGGG\" " .
       "--rev_primer \"AAAA\" " .
       "--fwd_mer_len 9 " .
	   "--min_con_depth 2 " .
	   "--diginorm_max NA " .
	   "--con_algo MSA " . 
       "--output_dir output_root/output_dir ",
       "_get_serial_command - $command_index" );
    $command_index++;
	
	# undo the sample_href override and default back to global
	delete $sample_href->{merged_file};
    delete $sample_href->{seq_type};
	delete $sample_href->{fwd_linker};
	
	# override the rev_linker when seq_type == 1 as in the first few tests
    $sample_href->{rev_linker} = "ATCG";
    is( MTToolbox::_get_serial_command($global_href, $sample_href),
       "my_exec " .
       "--id P1 " .
       "--seq_type 1 " .
       "--barcode AAAAA " .
       "--fwd_file samples_dir/fwd_file " .
       "--fwd_linker ATCG " .
	   "--rev_linker ATCG " .
       "--fwd_primer \"GGGG\" " .
       "--rev_primer \"AAAA\" " .
       "--fwd_mer_len 9 " .
	   "--min_con_depth 2 " .
	   "--diginorm_max NA " .
	   "--con_algo MSA " .
       "--output_dir output_root/output_dir ",
       "_get_serial_command - $command_index" );
    $command_index++;
    
	# undo the sample_href override and default back to global
	delete $sample_href->{rev_linker};
	
    # override the fwd primer
    $sample_href->{fwd_primer} = "TTTT";
    is( MTToolbox::_get_serial_command($global_href, $sample_href),
       "my_exec " .
       "--id P1 " .
       "--seq_type 1 " .
       "--barcode AAAAA " .
       "--fwd_file samples_dir/fwd_file " .
       "--fwd_linker ATCG " .
	   "--rev_linker GCTA " .
       "--fwd_primer \"TTTT\" " .
       "--rev_primer \"AAAA\" " .
       "--fwd_mer_len 9 " .
	   "--min_con_depth 2 " .
	   "--diginorm_max NA " .
	   "--con_algo MSA " .
       "--output_dir output_root/output_dir ",
       "_get_serial_command - $command_index" );
    $command_index++;
	
	# undo the sample_href override and default back to global
    delete $sample_href->{fwd_primer};
	
    # override the rev primer
    $sample_href->{rev_primer} = "TTTT";
    is( MTToolbox::_get_serial_command($global_href, $sample_href),
       "my_exec " .
       "--id P1 " .
       "--seq_type 1 " .
       "--barcode AAAAA " .
       "--fwd_file samples_dir/fwd_file " .
       "--fwd_linker ATCG " .
	   "--rev_linker GCTA " .
       "--fwd_primer \"GGGG\" " .
       "--rev_primer \"TTTT\" " .
       "--fwd_mer_len 9 " .
	   "--min_con_depth 2 " .
	   "--diginorm_max NA " .
	   "--con_algo MSA " .
       "--output_dir output_root/output_dir ",
       "_get_serial_command - $command_index" );
    $command_index++;
	
	# undo the sample_href override and default back to global
    delete $sample_href->{rev_primer};
	
    # override the fwd_mer_len
    $sample_href->{fwd_mer_len} = 5;
    is( MTToolbox::_get_serial_command($global_href, $sample_href),
       "my_exec " .
       "--id P1 " .
       "--seq_type 1 " .
       "--barcode AAAAA " .
       "--fwd_file samples_dir/fwd_file " .
       "--fwd_linker ATCG " .
	   "--rev_linker GCTA " .
       "--fwd_primer \"GGGG\" " .
       "--rev_primer \"AAAA\" " .
       "--fwd_mer_len 5 " .
	   "--min_con_depth 2 " .
	   "--diginorm_max NA " .
	   "--con_algo MSA " .
       "--output_dir output_root/output_dir ",
       "_get_serial_command - $command_index" );
    $command_index++;
	
	# undo the sample_href override and default back to global
    delete $sample_href->{fwd_mer_len};
	
	# override the min_con_depth
    $sample_href->{min_con_depth} = 5;
    is( MTToolbox::_get_serial_command($global_href, $sample_href),
       "my_exec " .
       "--id P1 " .
       "--seq_type 1 " .
       "--barcode AAAAA " .
       "--fwd_file samples_dir/fwd_file " .
       "--fwd_linker ATCG " .
	   "--rev_linker GCTA " .
       "--fwd_primer \"GGGG\" " .
       "--rev_primer \"AAAA\" " .
       "--fwd_mer_len 9 " .
	   "--min_con_depth 5 " .
	   "--diginorm_max NA " .
	   "--con_algo MSA " .
       "--output_dir output_root/output_dir ",
       "_get_serial_command - $command_index" );
    $command_index++;
	
	# undo the sample_href override and default back to global
    delete $sample_href->{min_con_depth};
	
	# override the diginorm_max
    $sample_href->{diginorm_max} = 5;
    is( MTToolbox::_get_serial_command($global_href, $sample_href),
       "my_exec " .
       "--id P1 " .
       "--seq_type 1 " .
       "--barcode AAAAA " .
       "--fwd_file samples_dir/fwd_file " .
       "--fwd_linker ATCG " .
	   "--rev_linker GCTA " .
       "--fwd_primer \"GGGG\" " .
       "--rev_primer \"AAAA\" " .
       "--fwd_mer_len 9 " .
	   "--min_con_depth 2 " .
	   "--diginorm_max 5 " .
	   "--con_algo MSA " .
       "--output_dir output_root/output_dir ",
       "_get_serial_command - $command_index" );
    $command_index++;
	
	# undo the sample_href override and default back to global
    delete $sample_href->{diginorm_max};
	
	# override the con_algo
    $sample_href->{con_algo} = 'NoMSA';
    is( MTToolbox::_get_serial_command($global_href, $sample_href),
       "my_exec " .
       "--id P1 " .
       "--seq_type 1 " .
       "--barcode AAAAA " .
       "--fwd_file samples_dir/fwd_file " .
       "--fwd_linker ATCG " .
	   "--rev_linker GCTA " .
       "--fwd_primer \"GGGG\" " .
       "--rev_primer \"AAAA\" " .
       "--fwd_mer_len 9 " .
	   "--min_con_depth 2 " .
	   "--diginorm_max NA " .
	   "--con_algo NoMSA " .
       "--output_dir output_root/output_dir ",
       "_get_serial_command - $command_index" );
    $command_index++;
	
	# undo the sample_href override and default back to global
    delete $sample_href->{con_algo};
	
	delete $sample_href->{output_dir}; #must go last (after override tests)
	
    # test automatic output dir based on sample_id and barcode
    is( MTToolbox::_get_serial_command($global_href, $sample_href),
       "my_exec " .
       "--id P1 " .
       "--seq_type 1 " .
       "--barcode AAAAA " .
       "--fwd_file samples_dir/fwd_file " .
       "--fwd_linker ATCG " .
	   "--rev_linker GCTA " .
       "--fwd_primer \"GGGG\" " .
       "--rev_primer \"AAAA\" " .
       "--fwd_mer_len 9 " .
	   "--min_con_depth 2 " .
	   "--diginorm_max NA " .
	   "--con_algo MSA " .
       "--output_dir output_root/P1_AAAAA ",
       "_get_serial_command - $command_index" );
    $command_index++;
    
    ## tests for seq_type == 2
    $global_href->{seq_type} = 2;
    $sample_href->{rev_file} = "rev_file";
    $sample_href->{rev_mer_len} = 4;
    is( MTToolbox::_get_serial_command($global_href, $sample_href),
       "my_exec " .
       "--id P1 " .
       "--seq_type 2 " .
       "--barcode AAAAA " .
       "--fwd_file samples_dir/fwd_file " .
       "--rev_file samples_dir/rev_file " .
       "--fwd_linker ATCG " .
	   "--rev_linker GCTA " .
       "--fwd_primer \"GGGG\" " .
       "--rev_primer \"AAAA\" " .
       "--fwd_mer_len 9 " .
       "--rev_mer_len 4 " .
	   "--min_con_depth 2 " .
	   "--diginorm_max NA " .
	   "--con_algo MSA " .
       "--output_dir output_root/P1_AAAAA ",
       "_get_serial_command - $command_index" );
    $command_index++;
    
    ## tests for seq_type == 3
    $global_href->{seq_type} = 3;
    $sample_href->{rev_file} = "rev_file";
    $sample_href->{rev_mer_len} = 4;
    is( MTToolbox::_get_serial_command($global_href, $sample_href),
       "my_exec " .
       "--id P1 " .
       "--seq_type 3 " .
       "--barcode AAAAA " .
       "--fwd_file samples_dir/fwd_file " .
       "--rev_file samples_dir/rev_file " .
       "--fwd_linker ATCG " .
	   "--rev_linker GCTA " .
       "--fwd_primer \"GGGG\" " .
       "--rev_primer \"AAAA\" " .
       "--fwd_mer_len 9 " .
       "--rev_mer_len 4 " .
	   "--min_con_depth 2 " .
	   "--diginorm_max NA " .
	   "--con_algo MSA " .
       "--output_dir output_root/P1_AAAAA ",
       "_get_serial_command - $command_index" );
    $command_index++;
    
    # test seq_type == 3 ignore merge file
    $sample_href->{merged_file} = "merged_file";
    is( MTToolbox::_get_serial_command($global_href, $sample_href),
       "my_exec " .
       "--id P1 " .
       "--seq_type 3 " .
       "--barcode AAAAA " .
       "--fwd_file samples_dir/fwd_file " .
       "--rev_file samples_dir/rev_file " .
       "--fwd_linker ATCG " .
	   "--rev_linker GCTA " .
       "--fwd_primer \"GGGG\" " .
       "--rev_primer \"AAAA\" " .
       "--fwd_mer_len 9 " .
       "--rev_mer_len 4 " .
	   "--min_con_depth 2 " .
	   "--diginorm_max NA " .
	   "--con_algo MSA " .
       "--output_dir output_root/P1_AAAAA ",
       "_get_serial_command - $command_index" );
    $command_index++;
	
	# adding global flash_params
    $global_href->{bin_sample_exec} = "my_exec";
    $global_href->{seq_type} = 1;
    $global_href->{split_samples_root_dir} = "samples_dir";
    $global_href->{output_root} = "output_root";
    $global_href->{fwd_linker} = "ATCG";
	$global_href->{rev_linker} = "GCTA";
    $global_href->{fwd_primer} = "GGGG";
    $global_href->{rev_primer} = "AAAA";
    $global_href->{fwd_mer_len} = 9;
	$global_href->{min_con_depth} = 2;
	$global_href->{diginorm_max} = "NA";
	$global_href->{con_algo} = 'MSA';
	$global_href->{flash_params} = {
									m => 10,
									M => 70,
									x => ".20",
									p => 33,
									r => 100,
									f => 100,
									s => 20
									};
    $sample_href->{sample_id} = "P1";
    $sample_href->{barcode} = "AAAAA";
    $sample_href->{fwd_file} = "fwd_file";
    $sample_href->{output_dir} = "output_dir";
	
	my $command = MTToolbox::_get_serial_command($global_href, $sample_href);
	is( $command =~m/--flash_m 10/, 1, "flash_params - $command_index: flash_m");
	is( $command =~m/--flash_Max 70/, 1, "flash_params - $command_index: flash_M");
	is( $command =~m/--flash_x \.20/, 1, "flash_params - $command_index: flash_x");
	is( $command =~m/--flash_p 33/, 1, "flash_params - $command_index: flash_p");
	is( $command =~m/--flash_r 100/, 1, "flash_params - $command_index: flash_r");
	is( $command =~m/--flash_f 100/, 1, "flash_params - $command_index: flash_f");
	is( $command =~m/--flash_s 20/, 1, "flash_params - $command_index: flash_s");
    $command_index++;
	
	# Add local flash params
	$sample_href->{flash_params} = {
									m => 20,
									M => 80,
									x => ".30",
									p => 43,
									r => 110,
									f => 110,
									s => 30
									};
	
	$command = MTToolbox::_get_serial_command($global_href, $sample_href);
	is( $command =~m/--flash_m 20/, 1, "flash_params - $command_index: flash_m");
	is( $command =~m/--flash_Max 80/, 1, "flash_params - $command_index: flash_M");
	is( $command =~m/--flash_x \.30/, 1, "flash_params - $command_index: flash_x");
	is( $command =~m/--flash_p 43/, 1, "flash_params - $command_index: flash_p");
	is( $command =~m/--flash_r 110/, 1, "flash_params - $command_index: flash_r");
	is( $command =~m/--flash_f 110/, 1, "flash_params - $command_index: flash_f");
	is( $command =~m/--flash_s 30/, 1, "flash_params - $command_index: flash_s");
}

# test _get_parallel_command
{
    my $global_href = ();
    my $sample_href = ();
    
    my $command_index = 0;
    
    # all globally assigned when possible
    $global_href->{bin_sample_exec} = "my_exec";
    $global_href->{seq_type} = 1;
    $global_href->{split_samples_root_dir} = "samples_dir";
    $global_href->{output_root} = "output_root";
    $global_href->{fwd_linker} = "ATCG";
	$global_href->{rev_linker} = "GCTA";
    $global_href->{fwd_primer} = "GGGG";
    $global_href->{rev_primer} = "AAAA";
    $global_href->{fwd_mer_len} = 9;
	$global_href->{min_con_depth} = 2;
	$global_href->{diginorm_max} = "NA";
	$global_href->{con_algo} = 'MSA';
    $sample_href->{sample_id} = "P1";
    $sample_href->{barcode} = "AAAAA";
    $sample_href->{fwd_file} = "fwd_file";
    $sample_href->{output_dir} = "output_dir";
    
    my $expected_command = "bsub -q week -n 1 " .
                           "-J MTBIN_P1 " .
                           "-o output_root/output_dir/log.out " .
                           "-e output_root/output_dir/log.err " .
                           MTToolbox::_get_serial_command($global_href, $sample_href);
    is( MTToolbox::_get_parallel_command($global_href, $sample_href),
        $expected_command,
        "_get_paralell_command()" );
}

# test new parameter --keep_tmp_files
{
	my $global_href = ();
    my $sample_href = ();
	
	# test the new keep_tmp_file parameter in serial command
	$global_href->{keep_tmp_files} = "N";
    $global_href->{bin_sample_exec} = "my_exec";
    $global_href->{seq_type} = 1;
    $global_href->{split_samples_root_dir} = "samples_dir";
    $global_href->{output_root} = "output_root";
    $global_href->{fwd_linker} = "ATCG";
	$global_href->{rev_linker} = "GCTA";
    $global_href->{fwd_primer} = "GGGG";
    $global_href->{rev_primer} = "AAAA";
    $global_href->{fwd_mer_len} = 9;
	$global_href->{min_con_depth} = 2;
	$global_href->{diginorm_max} = "NA";
	$global_href->{con_algo} = 'MSA';
    $sample_href->{sample_id} = "P1";
    $sample_href->{barcode} = "AAAAA";
    $sample_href->{fwd_file} = "fwd_file";
    $sample_href->{output_dir} = "output_dir";
    is( MTToolbox::_get_serial_command($global_href, $sample_href),
       "my_exec " .
       "--id P1 " .
       "--seq_type 1 " .
       "--barcode AAAAA " .
       "--fwd_file samples_dir/fwd_file " .
       "--fwd_linker ATCG " .
	   "--rev_linker GCTA " .
       "--fwd_primer \"GGGG\" " .
       "--rev_primer \"AAAA\" " .
       "--fwd_mer_len 9 " .
	   "--min_con_depth 2 " .
	   "--diginorm_max NA " .
	   "--con_algo MSA " .
       "--output_dir output_root/output_dir " .
	   "--keep_tmp_files N ",
       "_get_serial_command - keep_tmp_files" );
}



### These tests are harder because it requires the output files

# test _stall
{
    ;
}

# test _cat_outputs
{
    ;
}

# test _consolidate_summaries
{
    ;
}

# test _parse_summary_file
{
    ;
}

# test _gap_and_ambig_bases_analysis
{
    ;
}

# test _print_mt_depth_hist
{
	;
}

# test _quality filtering
{
	;
}

# remember to test what happens when I forget to add rev_mer_len -- these defualt values
# are set when I run _check_fwd_mer_len in _parse_params_file.  Check this when
# you test run()