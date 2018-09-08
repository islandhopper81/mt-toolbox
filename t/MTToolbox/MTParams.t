
use strict;
use warnings;

use MTToolbox::MTParams;
use MTToolbox::MyX::MTParams;
use File::Temp qw/ tempfile tempdir /;
use XML::Simple;
use File::Basename;
use Test::Exception;
use Test::Warn;
use Test::More tests => 402;


BEGIN { use_ok( 'MTToolbox::MTParams' );
        use_ok( 'MTToolbox::MyX::MTParams' );
        }

sub get_test_href;
sub get_test_file;


# some global variables
my ($global_href, $e, $fh, $filename);

### Tests ###

# test new
my $mt_params_g;  # this object can be used globally
{
    lives_ok( sub{$mt_params_g = MTToolbox::MTParams->new()},
             "new() - lives" );
    lives_ok( sub{$mt_params_g = MTToolbox::MTParams->new({href => get_test_href()})},
             "new(href) - lives" );
    lives_ok( sub{$mt_params_g = MTToolbox::MTParams->new({xml_file => get_test_file()})},
             "new(xml_file) - lives" );
}

# test set_params_href
{
    lives_ok( sub{$mt_params_g->set_params_href()},
             "new() - lives" );
    lives_ok( sub{$mt_params_g->set_params_href({href => get_test_href()})},
             "get_params_href(href) - lives" );
    lives_ok( sub{$mt_params_g->set_params_href({href => get_test_file()})},
             "get_params_href(href) - lives" );
}

# test get_params_href
{
    my $my_href = get_test_href();
    $mt_params_g->set_params_href({href => $my_href});
    is_deeply( $my_href, $mt_params_g->get_params_href(), "get_params_href()" );
}

# test parse_xml_file
{
    # TO DO
    ;
}

# test _get_required_global_tags
{
    my %required_global_tags = map { $_ => 1 } qw(
												output_root
												fwd_primer
												rev_primer
												seq_type
												sample
												parallelize
											);
    is_deeply( MTToolbox::MTParams::_get_required_global_tags(), \%required_global_tags,
              "_get_required_global_tags" );
}

# test _get_optional_global_tags
{
    my %optional_global_tags = map { $_ => 1 } qw(
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
    is_deeply( MTToolbox::MTParams::_get_optional_global_tags(), \%optional_global_tags,
              "_get_optional_global_tags" );
}

# test _get_required_sample_tags
{
    my %required_sample_tags = map { $_ => 1 } qw(
													   sample_id
												);
    is_deeply( MTToolbox::MTParams::_get_required_sample_tags(), \%required_sample_tags,
              "_get_required_sample_tags" );
}

# test _get_optional_sample_tags
{
    my %optional_sample_tags = map { $_ => 1 } qw(
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
													   con_algo
													   output_dir
													   fwd_file
													   rev_file
													   merged_file
													   flash_params
													   qc_params
													);
    is_deeply( MTToolbox::MTParams::_get_optional_sample_tags(), \%optional_sample_tags,
              "_get_optional_sample_tags" );
}

# test _get_optional_flash_params_tags
{
    my %optional_flash_params_tags = map { $_ => 1 } qw(
														m
														M
														x
														p
														r
														f
														s
													);
    is_deeply( MTToolbox::MTParams::_get_optional_flash_params_tags(),
			  \%optional_flash_params_tags,
              "_get_optional_flash_params_tags" );
}

# test _get_optional_qc_params_tags
{
    my %optional_qc_params_tags = map { $_ => 1 } qw(
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
    is_deeply( MTToolbox::MTParams::_get_optional_qc_params_tags(),
			  \%optional_qc_params_tags,
              "_get_optional_qc_params_tags" );
}

# test _is_required_tag_present
{
    my $required_href = {tag1 => 1, tag2 => 2, tag3 => {}};
    lives_ok( sub{ MTToolbox::MTParams::_is_required_tag_present('tag1', $required_href) },
             "_is_required_tag_present - lives" );
    
    is( MTToolbox::MTParams::_is_required_tag_present('tag1', $required_href), 1,
       "_is_required_tag_present" );
    
    throws_ok( sub{ MTToolbox::MTParams::_is_required_tag_present('tag', $required_href) },
              'MTToolbox::MyX::MTParams::MissingRequiredTag',
             "_is_required_tag_present - throws MissingRequiredTag" );
    
    throws_ok( sub{ MTToolbox::MTParams::_is_required_tag_present('tag', $required_href) },
              qr/Missing tag in params href/,
             "_is_required_tag_present - throws Missing tag in params href" );
    
    eval{ MTToolbox::MTParams::_is_required_tag_present('tag', $required_href) };
    my $e = Exception::Class->caught('MTToolbox::MyX::MTParams::MissingRequiredTag');
    is( $e->tag_name, 'tag', "field test - tag" );
}

# test _is_required_tag_defined
{
    my $href = {tag1 => 1, tag2 => 2, tag3 => {}};
    
    lives_ok( sub{ MTToolbox::MTParams::_is_required_tag_defined('tag1', $href) },
             "_is_required_tag_defined - lives" );
    
    is( MTToolbox::MTParams::_is_required_tag_defined('tag1', $href), 1,
       "_is_required_tag_defined" );
    
    throws_ok( sub{ MTToolbox::MTParams::_is_required_tag_defined('tag3', $href) },
              'MTToolbox::MyX::MTParams::MissingTagValue',
              "_is_required_tag_defined - thows MissingTagValue" );
    
    throws_ok( sub{ MTToolbox::MTParams::_is_required_tag_defined('tag3', $href) },
              qr/Missing tag value in params href/,
              "_is_required_tag_defined - thows Missing tag value in params" );
    
    eval{ MTToolbox::MTParams::_is_required_tag_defined('tag3', $href) };
    my $e = Exception::Class->caught('MTToolbox::MyX::MTParams::MissingTagValue');
    is( $e->tag_name, 'tag3', "field test - tag" );
    
    throws_ok( sub{ MTToolbox::MTParams::_is_required_tag_defined('tag', $href) },
              'MTToolbox::MyX::MTParams::MissingRequiredTag',
              "_is_required_tag_defined - throws MissingRequiredTag" );
}

# test _is_optional_tag_present
{
    my $href = {tag1 => 1, tag2 => 2, tag3 => {}};
    
    lives_ok( sub{ MTToolbox::MTParams::_is_optional_tag_present('tag1', $href) },
             "_is_optional_tag_present - lives" );
    
    is( MTToolbox::MTParams::_is_optional_tag_present('tag1', $href), 1,
       "_is_optional_tag_present(tag1)" );
    
    is( MTToolbox::MTParams::_is_optional_tag_present('tag', $href), 0,
       "_is_optional_tag_present(tag) - no" );
}

# test _is_optional_tag_defined
{
    my $href = {tag1 => 1, tag2 => 2, tag3 => {}};
    
    lives_ok( sub{ MTToolbox::MTParams::_is_optional_tag_defined('tag1', $href) },
             "_is_optional_tag_defined - lives" );
    
    is( MTToolbox::MTParams::_is_optional_tag_defined('tag1', $href), 1,
       "_is_optional_tag_defined" );
    
    is( MTToolbox::MTParams::_is_optional_tag_defined('tag3', $href), 0,
            "_is_optional_tag_defined - no" );
    
    is( MTToolbox::MTParams::_is_optional_tag_defined('tag', $href), 0,
            "_is_optional_tag_defined - no" );
}

# test _check_bin_sample_exec
{
    my ($fh, $filename) = tempfile();
    $global_href = {bin_sample_exec => $filename};
    
    throws_ok( sub{ MTToolbox::MTParams::_check_bin_sample_exec($global_href) },
              'MTToolbox::MyX::MTParams::MissingExec',
              "_check_bin_sample_exec - throws MTToolbox::MyX::MTParams::MissingExec" );
    
    throws_ok( sub{ MTToolbox::MTParams::_check_bin_sample_exec($global_href) },
              qr/Can't find executable file/,
              "_check_bin_sample_exec - throws Can't find executable file" );
    
    eval{ MTToolbox::MTParams::_check_bin_sample_exec($global_href) };
    my $e = Exception::Class->caught('MTToolbox::MyX::MTParams::MissingExec');
    is( $e->exec_name, $filename, "MissingExec->exec_name eval test" );
    
    # give the temp file exe privlages
    chmod 0755, "$filename";
    lives_ok( sub{ MTToolbox::MTParams::_check_bin_sample_exec($global_href) },
             "_check_bin_sample_exec - lives" );
    is( MTToolbox::MTParams::_check_bin_sample_exec($global_href), 1,
             "_check_bin_sample_exec - 1" );
    
    $global_href = {};
    lives_ok( sub{ MTToolbox::MTParams::_check_bin_sample_exec($global_href) },
             "_check_bin_sample_exec - lives" );
    $global_href = {};
    is( MTToolbox::MTParams::_check_bin_sample_exec($global_href), 1,
             "_check_bin_sample_exec - 1" );
}

# test _check_split_samples_root_dir
{
    my $dir = tempdir( CLEANUP => 1 );
    my $global_href = {split_samples_root_dir => $dir};
    
    lives_ok( sub{ MTToolbox::MTParams::_check_split_samples_root_dir($global_href) },
             "_check_split_samples_root_dir - lives" );
    
    is( MTToolbox::MTParams::_check_split_samples_root_dir($global_href), 1,
             "_check_split_samples_root_dir - 1" );
    
    my ($fh, $filename) = tempfile();
    $global_href = {split_samples_root_dir => $filename};
    throws_ok( sub{ MTToolbox::MTParams::_check_split_samples_root_dir($global_href) },
              'MyX::Generic::DoesNotExist::Dir',
              "_check_split_samples_root_dir - throws MyX::Generic::DoesNotExist::Dir" );
    
    throws_ok( sub{ MTToolbox::MTParams::_check_split_samples_root_dir($global_href) },
              qr/Directory does not exist/,
              "_check_split_samples_root_dir - throws Directory does not exist" );
    
    eval{ MTToolbox::MTParams::_check_split_samples_root_dir($global_href) };
    my $e = Exception::Class->caught('MyX::Generic::DoesNotExist::Dir');
    is( $e->dir_name, $filename, "DoesNotExist::Dir->dir_name eval test" );
}

# test _check_output_root
{
    my $dir = tempdir( CLEANUP => 1 );
    my $global_href = {output_root => $dir};
    
    warnings_are { MTToolbox::MTParams::_check_output_root($global_href) } [],
             "_check_output_root - no warnings";
             
    is( MTToolbox::MTParams::_check_output_root($global_href), 1,
             "_check_output_root - 1" );
    
    my ($fh, $tempfile) = tempfile();
    close($fh);
    $global_href = {output_root => $tempfile};
    
    warnings_are {MTToolbox::MTParams::_check_output_root($global_href) }
        ["Making output_root directory: $tempfile"],
        "_check_output_root - Making output_root dir";

    is( MTToolbox::MTParams::_check_output_root($global_href), 1,
       "_check_output_root - 1" );
}

# test _check_fwd_linker
{
    my $global_href;
    
    # check when the forward linker is not defined
    is( MTToolbox::MTParams::_check_fwd_linker($global_href), "",
	   "fwd_linker not defined" );
    $global_href = {fwd_linker => {}};
    is( MTToolbox::MTParams::_check_fwd_linker($global_href), "",
	   "fwd_linker not defined {}" );
    
    # check warning for long linkers
    $global_href = {fwd_linker => "ATGATGATTAGTATATAGG"};
    throws_ok( sub{MTToolbox::MTParams::_check_fwd_linker($global_href)},
               'MTToolbox::MyX::MTParams::LongLinker',
               "_check_fwd_linker - throws FWD linker is longer than 10 bp" );
    throws_ok( sub{MTToolbox::MTParams::_check_fwd_linker($global_href)},
               qr/FWD linker is longer than 10 bp/,
               "_check_fwd_linker - throws FWD linker is longer than 10 bp" );
    eval{ MTToolbox::MTParams::_check_fwd_linker($global_href) };
    my $e = Exception::Class->caught('MTToolbox::MyX::MTParams::LongLinker');
    is( $e->orientation, 'fwd', "LongLinker->orientation eval test" );
    is( $e->length, 19, "LongLinker->length eval test" );
    
    # check errors for bad characters
    $global_href = {fwd_linker => "AGTN"};
    throws_ok( sub{ MTToolbox::MTParams::_check_fwd_linker($global_href) },
              'MTToolbox::MyX::MTParams::UnrecognizedChar',
              "_check_fwd_linker(AGTN) - throws UnrecognizedChar" );
    throws_ok( sub{ MTToolbox::MTParams::_check_fwd_linker($global_href) },
              qr/Unrecognized Char in fwd linker/,
              "_check_fwd_linker(AGTN) - throws Unrecognized Char in fwd linker" );
    eval{ MTToolbox::MTParams::_check_fwd_linker($global_href) };
    $e = Exception::Class->caught('MTToolbox::MyX::MTParams::UnrecognizedChar');
    is( $e->char, 'N', "UnrecognizedChar->char eval test" );
    
    # check a normal forward linker
    $global_href = {fwd_linker => "ATCG"};
    is( MTToolbox::MTParams::_check_fwd_linker($global_href), 1,
	   "_check_fwd_linker(ATCG)" );
    
    # check changes from lower to upper case
    $global_href = {fwd_linker => "atgc"};
    is( MTToolbox::MTParams::_check_fwd_linker($global_href), 1,
	   "_check_fwd_linker(atcg)" );
    is( $global_href->{fwd_linker}, "ATGC", "check change to upper case" );
}

# test _check_rev_linker
{
    my $global_href;
    
    # check when the forward linker is not defined
    is( MTToolbox::MTParams::_check_rev_linker($global_href), "",
	   "rev_linker not defined" );
    $global_href = {rev_linker => {}};
    is( MTToolbox::MTParams::_check_rev_linker($global_href), "",
	   "rev_linker not defined {}" );
    
    # check warning for long linkers
    $global_href = {rev_linker => "ATGATGATTAGTATATAGG"};
    throws_ok( sub{MTToolbox::MTParams::_check_rev_linker($global_href)},
               'MTToolbox::MyX::MTParams::LongLinker',
               "_check_rev_linker - throws rev linker is longer than 10 bp" );
    throws_ok( sub{MTToolbox::MTParams::_check_rev_linker($global_href)},
               qr/REV linker is longer than 10 bp/,
               "_check_rev_linker - throws rev linker is longer than 10 bp" );
    eval{ MTToolbox::MTParams::_check_rev_linker($global_href) };
    my $e = Exception::Class->caught('MTToolbox::MyX::MTParams::LongLinker');
    is( $e->orientation, 'rev', "LongLinker->orientation eval test" );
    is( $e->length, 19, "LongLinker->length eval test" );
    
    # check errors for bad characters
    $global_href = {rev_linker => "AGTN"};
    throws_ok( sub{ MTToolbox::MTParams::_check_rev_linker($global_href) },
              'MTToolbox::MyX::MTParams::UnrecognizedChar',
              "_check_rev_linker(AGTN) - throws UnrecognizedChar" );
    throws_ok( sub{ MTToolbox::MTParams::_check_rev_linker($global_href) },
              qr/Unrecognized Char in rev linker/,
              "_check_rev_linker(AGTN) - throws Unrecognized Char in rev linker" );
    eval{ MTToolbox::MTParams::_check_rev_linker($global_href) };
    $e = Exception::Class->caught('MTToolbox::MyX::MTParams::UnrecognizedChar');
    is( $e->char, 'N', "UnrecognizedChar->char eval test" );
    
    # check a normal forward linker
    $global_href = {rev_linker => "ATCG"};
    is( MTToolbox::MTParams::_check_rev_linker($global_href), 1,
	   "_check_rev_linker(ATCG)" );
    
    # check changes from lower to upper case
    $global_href = {rev_linker => "atgc"};
    is( MTToolbox::MTParams::_check_rev_linker($global_href), 1,
	   "_check_rev_linker(atcg)" );
    is( $global_href->{rev_linker}, "ATGC", "check change to upper case" );
}

# test _check_fwd_primer
{
    # This method doesn't do any checks right now.  It could in the future
    is( MTToolbox::MTParams::_check_fwd_primer(), 1, "_check_fwd_primer" );
}

# test _check_rev_primer
{
    # This method doesn't do any checks right now.  It could in the future
    is( MTToolbox::MTParams::_check_rev_primer(), 1, "_check_rev_primer" );
}

# test _check_fwd_mer_len
{
	# check with undef
	my $global_href = {};
	is( MTToolbox::MTParams::_check_fwd_mer_len($global_href), 1,
            "_check_fwd_mer_len() DEFAULT" );
	is( $global_href->{fwd_mer_len}, 8, "fwd_mer_len == 8" );
	$global_href = {fwd_mer_len => {}};
	is( MTToolbox::MTParams::_check_fwd_mer_len($global_href), 1,
            "_check_fwd_mer_len() DEFAULT" );
	is( $global_href->{fwd_mer_len}, 8, "fwd_mer_len == 8" );
	
    # check errors with non-digits
    $global_href = {fwd_mer_len => "a"};
    throws_ok( sub{ MTToolbox::MTParams::_check_fwd_mer_len($global_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_fwd_mer_len(a) - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_fwd_mer_len($global_href) },
              qr/fwd mer len must be a digit/,
              "_check_fwd_mer_len(a) - throws fwd mer len must be a digit" );
    eval{ MTToolbox::MTParams::_check_fwd_mer_len($global_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test" );
    
    $global_href = {fwd_mer_len => "-"};
    throws_ok( sub{ MTToolbox::MTParams::_check_fwd_mer_len($global_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_fwd_mer_len(-) - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_fwd_mer_len($global_href) },
              qr/fwd mer len must be a digit/,
              "_check_fwd_mer_len(-) - throws fwd mer len must be a digit" );
    eval{ MTToolbox::MTParams::_check_fwd_mer_len($global_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, '-', "MustBeDigit->value eval test" );
    
    # check errors with non-positive numbers
    $global_href = {fwd_mer_len => -1};
    throws_ok( sub{ MTToolbox::MTParams::_check_fwd_mer_len($global_href) },
              'MyX::Generic::Digit::TooSmall',
              "_check_fwd_mer_len(-1) - throws Digit::TooSmall" );
    throws_ok( sub{ MTToolbox::MTParams::_check_fwd_mer_len($global_href) },
              qr/fwd mer len must be > 0/,
              "_check_fwd_mer_len(-1) - throws fwd mer len must be > 0" );
    eval{ MTToolbox::MTParams::_check_fwd_mer_len($global_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, '-1', "TooSmall->value eval test" );
    is( $e->MIN, 0, "TooSmall->MIN eval test" );
    
    # check a good fwd_mer_len value
    $global_href = {fwd_mer_len => 0};
    is( MTToolbox::MTParams::_check_fwd_mer_len($global_href), 1,
            "_check_fwd_mer_len(0)" );
    $global_href = {fwd_mer_len => 9};
    is( MTToolbox::MTParams::_check_fwd_mer_len($global_href), 1,
       "_check_fwd_mer_len(9)" );
}

# test _check_rev_mer_len
{
	# check with undef
	my $global_href = {};
	is( MTToolbox::MTParams::_check_rev_mer_len($global_href), 1,
            "_check_rev_mer_len() DEFAULT" );
	is( $global_href->{rev_mer_len}, 5, "rev_mer_len == 5" );
	$global_href = {rev_mer_len => {}};
	is( MTToolbox::MTParams::_check_rev_mer_len($global_href), 1,
            "_check_rev_mer_len() DEFAULT" );
	is( $global_href->{rev_mer_len}, 5, "rev_mer_len == 5" );
	
    # check errors with non-digits
    $global_href = {rev_mer_len => "a"};
    throws_ok( sub{ MTToolbox::MTParams::_check_rev_mer_len($global_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_rev_mer_len(a) - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_rev_mer_len($global_href) },
              qr/rev mer len must be a digit/,
              "_check_rev_mer_len(a) - throws rev mer len must be a digit" );
    eval{ MTToolbox::MTParams::_check_rev_mer_len($global_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test" );
    
    $global_href = {rev_mer_len => "-"};
    throws_ok( sub{ MTToolbox::MTParams::_check_rev_mer_len($global_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_rev_mer_len(-) - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_rev_mer_len($global_href) },
              qr/rev mer len must be a digit/,
              "_check_rev_mer_len(-) - throws rev mer len must be a digit" );
    eval{ MTToolbox::MTParams::_check_rev_mer_len($global_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, '-', "MustBeDigit->value eval test" );
    
    # check errors with non-positive numbers
    $global_href = {rev_mer_len => -1};
    throws_ok( sub{ MTToolbox::MTParams::_check_rev_mer_len($global_href) },
              'MyX::Generic::Digit::TooSmall',
              "_check_rev_mer_len(-1) - throws Digit::TooSmall" );
    throws_ok( sub{ MTToolbox::MTParams::_check_rev_mer_len($global_href) },
              qr/rev mer len must be > 0/,
              "_check_rev_mer_len(-1) - throws rev mer len must be > 0" );
    eval{ MTToolbox::MTParams::_check_rev_mer_len($global_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, '-1', "TooSmall->value eval test" );
    is( $e->MIN, 0, "TooSmall->MIN eval test" );
    
    # check a good rev_mer_len value
    $global_href = {rev_mer_len => 0};
    is( MTToolbox::MTParams::_check_rev_mer_len($global_href), 1,
            "_check_rev_mer_len(0)" );
    $global_href = {rev_mer_len => 9};
    is( MTToolbox::MTParams::_check_rev_mer_len($global_href), 1,
       "_check_rev_mer_len(9)" );
}

# test _check_fwd_max_shifts
{
	# check with undef
	my $global_href = {};
	is( MTToolbox::MTParams::_check_fwd_max_shifts($global_href), 1,
            "_check_fwd_max_shifts() DEFAULT" );
	is( $global_href->{fwd_max_shifts}, 5, "fwd_max_shifts == 5" );
	$global_href = {fwd_max_shifts => {}};
	is( MTToolbox::MTParams::_check_fwd_max_shifts($global_href), 1,
            "_check_fwd_max_shifts() DEFAULT" );
	is( $global_href->{fwd_max_shifts}, 5, "fwd_max_shifts == 5" );
	
    # check errors with non-digits
    $global_href = {fwd_max_shifts => "a"};
    throws_ok( sub{ MTToolbox::MTParams::_check_fwd_max_shifts($global_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_fwd_max_shifts(a) - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_fwd_max_shifts($global_href) },
              qr/fwd max shifts must be a digit/,
              "_check_fwd_max_shifts(a) - throws fwd max shifts must be a digit" );
    eval{ MTToolbox::MTParams::_check_fwd_max_shifts($global_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test" );
    
    $global_href = {fwd_max_shifts => "-"};
    throws_ok( sub{ MTToolbox::MTParams::_check_fwd_max_shifts($global_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_fwd_max_shifts(-) - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_fwd_max_shifts($global_href) },
              qr/fwd max shifts must be a digit/,
              "_check_fwd_max_shifts(-) - throws fwd max shifts must be a digit" );
    eval{ MTToolbox::MTParams::_check_fwd_max_shifts($global_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, '-', "MustBeDigit->value eval test" );
    
    # check errors with non-positive numbers
    $global_href = {fwd_max_shifts => -1};
    throws_ok( sub{ MTToolbox::MTParams::_check_fwd_max_shifts($global_href) },
              'MyX::Generic::Digit::TooSmall',
              "_check_fwd_max_shifts(-1) - throws Digit::TooSmall" );
    throws_ok( sub{ MTToolbox::MTParams::_check_fwd_max_shifts($global_href) },
              qr/fwd max shifts must be > 0/,
              "_check_fwd_max_shifts(-1) - throws fwd max shifts must be > 0" );
    eval{ MTToolbox::MTParams::_check_fwd_max_shifts($global_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, '-1', "TooSmall->value eval test" );
    is( $e->MIN, 0, "TooSmall->MIN eval test" );
    
    # check a good rev_mer_len value
    $global_href = {fwd_max_shifts => 0};
    is( MTToolbox::MTParams::_check_fwd_max_shifts($global_href), 1,
            "_check_rev_mer_len(0)" );
    $global_href = {fwd_max_shifts => 9};
    is( MTToolbox::MTParams::_check_fwd_max_shifts($global_href), 1,
       "_check_rev_mer_len(9)" );
}

# test _check_rev_max_shifts
{
	# check with undef
	my $global_href = {};
	is( MTToolbox::MTParams::_check_rev_max_shifts($global_href), 1,
            "_check_rev_max_shifts() DEFAULT" );
	is( $global_href->{rev_max_shifts}, 5, "rev_max_shifts == 5" );
	$global_href = {rev_max_shifts => {}};
	is( MTToolbox::MTParams::_check_rev_max_shifts($global_href), 1,
            "_check_rev_max_shifts() DEFAULT" );
	is( $global_href->{rev_max_shifts}, 5, "rev_max_shifts == 5" );
	
    # check errors with non-digits
    $global_href = {rev_max_shifts => "a"};
    throws_ok( sub{ MTToolbox::MTParams::_check_rev_max_shifts($global_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_rev_max_shifts(a) - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_rev_max_shifts($global_href) },
              qr/rev max shifts must be a digit/,
              "_check_rev_max_shifts(a) - throws rev max shifts must be a digit" );
    eval{ MTToolbox::MTParams::_check_rev_max_shifts($global_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test" );
    
    $global_href = {rev_max_shifts => "-"};
    throws_ok( sub{ MTToolbox::MTParams::_check_rev_max_shifts($global_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_rev_max_shifts(-) - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_rev_max_shifts($global_href) },
              qr/rev max shifts must be a digit/,
              "_check_rev_max_shifts(-) - throws rev max shifts must be a digit" );
    eval{ MTToolbox::MTParams::_check_rev_max_shifts($global_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, '-', "MustBeDigit->value eval test" );
    
    # check errors with non-positive numbers
    $global_href = {rev_max_shifts => -1};
    throws_ok( sub{ MTToolbox::MTParams::_check_rev_max_shifts($global_href) },
              'MyX::Generic::Digit::TooSmall',
              "_check_rev_max_shifts(-1) - throws Digit::TooSmall" );
    throws_ok( sub{ MTToolbox::MTParams::_check_rev_max_shifts($global_href) },
              qr/rev max shifts must be > 0/,
              "_check_rev_max_shifts(-1) - throws rev max shifts must be > 0" );
    eval{ MTToolbox::MTParams::_check_rev_max_shifts($global_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, '-1', "TooSmall->value eval test" );
    is( $e->MIN, 0, "TooSmall->MIN eval test" );
    
    # check a good rev_mer_len value
    $global_href = {rev_max_shifts => 0};
    is( MTToolbox::MTParams::_check_rev_max_shifts($global_href), 1,
            "_check_rev_mer_len(0)" );
    $global_href = {rev_max_shifts => 9};
    is( MTToolbox::MTParams::_check_rev_max_shifts($global_href), 1,
       "_check_rev_mer_len(9)" );
}

# test _check_seq_type
{
    # check errors with non-digits
    my $global_href = {seq_type => "a"};
    throws_ok( sub{ MTToolbox::MTParams::_check_seq_type($global_href) },
              'MTToolbox::MyX::MTParams::BadSeqType',
              "_check_seq_type(a) - thows BadSeqType" );
    throws_ok( sub{ MTToolbox::MTParams::_check_seq_type($global_href) },
              qr/seq_type must be either 1, 2, or 3/,
              "_check_seq_type(a) - thows seq_type must be either 1, 2, or 3" );
    eval{ MTToolbox::MTParams::_check_seq_type($global_href) };
    $e = Exception::Class->caught('MTToolbox::MyX::MTParams::BadSeqType');
    is( $e->value, 'a', "BadSeqType->value eval test" );
    
    # check errors with wrong numbers
    $global_href = {seq_type => 0};
    throws_ok( sub{ MTToolbox::MTParams::_check_seq_type($global_href) },
              'MTToolbox::MyX::MTParams::BadSeqType',
              "_check_seq_type(0) - thows BadSeqType" );
    throws_ok( sub{ MTToolbox::MTParams::_check_seq_type($global_href) },
              qr/seq_type must be either 1, 2, or 3/,
              "_check_seq_type(0) - thows seq_type must be either 1, 2, or 3" );
    eval{ MTToolbox::MTParams::_check_seq_type($global_href) };
    $e = Exception::Class->caught('MTToolbox::MyX::MTParams::BadSeqType');
    is( $e->value, '0', "BadSeqType->value eval test" );
    
    $global_href = {seq_type => -1};
    throws_ok( sub{ MTToolbox::MTParams::_check_seq_type($global_href) },
              'MTToolbox::MyX::MTParams::BadSeqType',
              "_check_seq_type(-1) - thows BadSeqType" );
    throws_ok( sub{ MTToolbox::MTParams::_check_seq_type($global_href) },
              qr/seq_type must be either 1, 2, or 3/,
              "_check_seq_type(-1) - thows seq_type must be either 1, 2, or 3" );
    eval{ MTToolbox::MTParams::_check_seq_type($global_href) };
    $e = Exception::Class->caught('MTToolbox::MyX::MTParams::BadSeqType');
    is( $e->value, '-1', "BadSeqType->value eval test" );
    
    $global_href = {seq_type => 4};
    throws_ok( sub{ MTToolbox::MTParams::_check_seq_type($global_href) },
              'MTToolbox::MyX::MTParams::BadSeqType',
              "_check_seq_type(4) - thows BadSeqType" );
    throws_ok( sub{ MTToolbox::MTParams::_check_seq_type($global_href) },
              qr/seq_type must be either 1, 2, or 3/,
              "_check_seq_type(4) - thows seq_type must be either 1, 2, or 3" );
    eval{ MTToolbox::MTParams::_check_seq_type($global_href) };
    $e = Exception::Class->caught('MTToolbox::MyX::MTParams::BadSeqType');
    is( $e->value, '4', "BadSeqType->value eval test" );
    
    # check with correct values
    $global_href = {seq_type => 1};
    is( MTToolbox::MTParams::_check_seq_type($global_href), 1, "_check_seq_type(1)" );
    $global_href = {seq_type => 2};
    is( MTToolbox::MTParams::_check_seq_type($global_href), 1, "_check_seq_type(2)" );
    $global_href = {seq_type => 3};
    is( MTToolbox::MTParams::_check_seq_type($global_href), 1, "_check_seq_type(3)" );
}

# test _check_parallelize -- NOTE: Not fully tested!!
{
    # check for bad characters
    my $global_href = {parallelize => "NO dude"};
    throws_ok( sub{ MTToolbox::MTParams::_check_parallelize($global_href) },
              'MyX::Generic::BadValue',
              "_check_parallelize(NO dude) - throws BadValue" );
    throws_ok( sub{ MTToolbox::MTParams::_check_parallelize($global_href) },
              qr/Unrecognized value in parallelize/,
              "_check_parallelize(NO dude) - throws Unrecognized value in parallelize" );
    eval{ MTToolbox::MTParams::_check_parallelize($global_href) };
    $e = Exception::Class->caught('MyX::Generic::BadValue');
    is( $e->value, 'NO DUDE', "BadValue->value eval test" );
    
    $global_href = {parallelize => "M"};
    throws_ok( sub{ MTToolbox::MTParams::_check_parallelize($global_href) },
              'MyX::Generic::BadValue',
              "_check_parallelize(M) - throws BadValue" );
    throws_ok( sub{ MTToolbox::MTParams::_check_parallelize($global_href) },
              qr/Unrecognized value in parallelize/,
              "_check_parallelize(M) - throws Unrecognized value in parallelize" );
    eval{ MTToolbox::MTParams::_check_parallelize($global_href) };
    $e = Exception::Class->caught('MyX::Generic::BadValue');
    is( $e->value, 'M', "BadValue->value eval test" );
    
    # check some of the good values
    $global_href = {parallelize => "N"};
    MTToolbox::MTParams::_check_parallelize($global_href);
    is( $global_href->{parallelize}, "N", "_check_parallelize(N)" );
    
    $global_href = {parallelize => "NO"};
    MTToolbox::MTParams::_check_parallelize($global_href);
    is( $global_href->{parallelize}, "NO", "_check_parallelize(NO)" );
    
    $global_href = {parallelize => "n"};
    MTToolbox::MTParams::_check_parallelize($global_href);
    is( $global_href->{parallelize}, "n", "_check_parallelize(n)" );
    
    $global_href = {parallelize => "no"};
    MTToolbox::MTParams::_check_parallelize($global_href);
    is( $global_href->{parallelize}, "no", "_check_parallelize(no)" );
    
    # it is much harder to test this because in order to get a 1 (ie true values)
    # you have to have a system with the bsub command.  So _check_parallelize
    # is not fully tested
}

# test _check_keep_tmp_files
{
	# check for bad characters
    my $global_href = {keep_tmp_files => "NO dude"};
    throws_ok( sub{ MTToolbox::MTParams::_check_keep_tmp_files($global_href) },
              'MyX::Generic::BadValue',
              "_check_keep_tmp_files(NO dude) - throws BadValue" );
    throws_ok( sub{ MTToolbox::MTParams::_check_keep_tmp_files($global_href) },
              qr/Unrecognized value in keep_tmp_files/,
              "_check_keep_tmp_files(NO dude) - throws Unrecognized value in keep_tmp_files" );
    eval{ MTToolbox::MTParams::_check_keep_tmp_files($global_href) };
    $e = Exception::Class->caught('MyX::Generic::BadValue');
    is( $e->value, 'NO DUDE', "BadValue->value eval test" );
    
    $global_href = {keep_tmp_files => "M"};
    throws_ok( sub{ MTToolbox::MTParams::_check_keep_tmp_files($global_href) },
              'MyX::Generic::BadValue',
              "_check_keep_tmp_files(M) - throws BadValue" );
    throws_ok( sub{ MTToolbox::MTParams::_check_keep_tmp_files($global_href) },
              qr/Unrecognized value in keep_tmp_files/,
              "_check_keep_tmp_files(M) - throws Unrecognized value in keep_tmp_files" );
    eval{ MTToolbox::MTParams::_check_keep_tmp_files($global_href) };
    $e = Exception::Class->caught('MyX::Generic::BadValue');
    is( $e->value, 'M', "BadValue->value eval test" );
    
    # check some of the good values
    $global_href = {keep_tmp_files => "N"};
    MTToolbox::MTParams::_check_keep_tmp_files($global_href);
    is( $global_href->{keep_tmp_files}, "N", "_check_keep_tmp_files(N)" );
    
    $global_href = {keep_tmp_files => "NO"};
    MTToolbox::MTParams::_check_keep_tmp_files($global_href);
    is( $global_href->{keep_tmp_files}, "NO", "_check_keep_tmp_files(NO)" );
    
    $global_href = {keep_tmp_files => "n"};
    MTToolbox::MTParams::_check_keep_tmp_files($global_href);
    is( $global_href->{keep_tmp_files}, "n", "_check_keep_tmp_files(n)" );
    
    $global_href = {keep_tmp_files => "no"};
    MTToolbox::MTParams::_check_keep_tmp_files($global_href);
    is( $global_href->{keep_tmp_files}, "no", "_check_keep_tmp_files(no)" );
}

# test _check_min_con_depth
{
    # check errors with non-digits
    my $global_href = {min_con_depth => "a"};
    throws_ok( sub{ MTToolbox::MTParams::_check_min_con_depth($global_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_min_con_depth(a) - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_min_con_depth($global_href) },
              qr/min_con_depth must be an int > 1/,
              "_check_min_con_depth(a) - throws min_con_depth must be an int > 1" );
    eval{ MTToolbox::MTParams::_check_min_con_depth($global_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test" );
    
    $global_href = {min_con_depth => "-"};
    throws_ok( sub{ MTToolbox::MTParams::_check_min_con_depth($global_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_min_con_depth(-) - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_min_con_depth($global_href) },
              qr/min_con_depth must be an int > 1/,
              "_check_min_con_depth(-) - throws min_con_depth must be an int > 1" );
    eval{ MTToolbox::MTParams::_check_min_con_depth($global_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, '-', "MustBeDigit->value eval test" );
    
    # check errors with non-positive numbers
    $global_href = {min_con_depth => -1};
    throws_ok( sub{ MTToolbox::MTParams::_check_min_con_depth($global_href) },
              'MyX::Generic::Digit::TooSmall',
              "_check_min_con_depth(-1) - throws Digit::TooSmall" );
    throws_ok( sub{ MTToolbox::MTParams::_check_min_con_depth($global_href) },
              qr/min_con_depth must be an int > 1/,
              "_check_min_con_depth(-1) - throws min_con_depth must be an int > 1" );
    eval{ MTToolbox::MTParams::_check_min_con_depth($global_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, '-1', "TooSmall->value eval test" );
    is( $e->MIN, 2, "TooSmall->MIN eval test" );
    
    # check a good min_con_depth value
    $global_href = {min_con_depth => 2};
    is( MTToolbox::MTParams::_check_min_con_depth($global_href), 1,
        "_check_min_con_depth(0)" );
    $global_href = {min_con_depth => 9};
    is( MTToolbox::MTParams::_check_min_con_depth($global_href), 1,
       "_check_min_con_depth(9)" );
}

# test _check_diginorm_max
{
    # check errors with non-digits
    my $global_href = {diginorm_max => "a"};
    throws_ok( sub{ MTToolbox::MTParams::_check_diginorm_max($global_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_diginorm_max(a) - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_diginorm_max($global_href) },
              qr/diginorm_max must be an int > 1 or NA/,
              "_check_diginorm_max(a) - throws diginorm_max must be an int > 1 or NA" );
    eval{ MTToolbox::MTParams::_check_diginorm_max($global_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test" );
    
    $global_href = {diginorm_max => "-"};
    throws_ok( sub{ MTToolbox::MTParams::_check_diginorm_max($global_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_diginorm_max(-) - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_diginorm_max($global_href) },
              qr/diginorm_max must be an int > 1 or NA/,
              "_check_diginorm_max(-) - throws diginorm_max must be an int > 1 or NA" );
    eval{ MTToolbox::MTParams::_check_diginorm_max($global_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, '-', "MustBeDigit->value eval test" );
    
    # check errors with non-positive numbers
    $global_href = {diginorm_max => -1};
    throws_ok( sub{ MTToolbox::MTParams::_check_diginorm_max($global_href) },
              'MyX::Generic::Digit::TooSmall',
              "_check_diginorm_max(-1) - throws Digit::TooSmall" );
    throws_ok( sub{ MTToolbox::MTParams::_check_diginorm_max($global_href) },
              qr/diginorm_max must be an int > 1 or NA/,
              "_check_diginorm_max(-1) - throws diginorm_max must be an int > 1 or NA" );
    eval{ MTToolbox::MTParams::_check_diginorm_max($global_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, '-1', "TooSmall->value eval test" );
    is( $e->MIN, 2, "TooSmall->MIN eval test" );
    
    # check a good diginorm_max value
    $global_href = {diginorm_max => 2};
    is( MTToolbox::MTParams::_check_diginorm_max($global_href), 1,
        "_check_diginorm_max(0)" );
    $global_href = {diginorm_max => 9};
    is( MTToolbox::MTParams::_check_diginorm_max($global_href), 1,
       "_check_diginorm_max(9)" );
    $global_href = {diginorm_max => "NA"};
    is( MTToolbox::MTParams::_check_diginorm_max($global_href), 1,
       "_check_diginorm_max(NA)" );
}

# test _check_con_algo
{
	# check errors with non-acceptable values
    my $global_href = {con_algo => "blah"};
    throws_ok( sub{ MTToolbox::MTParams::_check_con_algo($global_href) },
              'MyX::Generic::BadValue',
              "_check_con_algo(blah) - throws BadValue" );
    throws_ok( sub{ MTToolbox::MTParams::_check_con_algo($global_href) },
              qr/Unrecognized value in con_algo.  Options: Muscle | Clustalw | NoMSA'/,
              '_check_con_algo(blah) - Unrecognized value in con_algo.  Options: Muscle | Clustalw | NoMSA' );
    eval{ MTToolbox::MTParams::_check_con_algo($global_href) };
    $e = Exception::Class->caught('MyX::Generic::BadValue');
    is( $e->value, 'blah', "BadValue->value eval test" );
    
    # check a good con_algo value
    $global_href = {con_algo => 'Muscle'};
    is( MTToolbox::MTParams::_check_con_algo($global_href), 1,
        "_check_con_algo(Muscle)" );
    $global_href = {con_algo => 'NoMSA'};
    is( MTToolbox::MTParams::_check_con_algo($global_href), 1,
       "_check_con_algo(NoMSA)" );
    $global_href = {con_algo => "Clustalw"};
    is( MTToolbox::MTParams::_check_con_algo($global_href), 1,
       "_check_con_algo(clustalw)" );
}

# test _check_min_overlap
{
	my $flash_params_href = {};
	
	lives_ok( sub{ MTToolbox::MTParams::_check_min_overlap($flash_params_href) },
            "_check_min_overlap(does not exist) - lives (it is optional)" );
	
    $flash_params_href = {m => 'a'};
	throws_ok( sub{ MTToolbox::MTParams::_check_min_overlap($flash_params_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_min_overlap(m => a) - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_min_overlap($flash_params_href) },
              qr/FLASH param 'm' is not a digit/,
              "_check_min_overlap(m => a) - throws FLASH param 'm' is not a digit" );
    eval{ MTToolbox::MTParams::_check_min_overlap($flash_params_href) };
    my $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test" );
    
	$flash_params_href = {m => -1};
	throws_ok( sub{ MTToolbox::MTParams::_check_min_overlap($flash_params_href) },
              'MyX::Generic::Digit::TooSmall',
              "_check_min_overlap(m => -1) - throws TooSmall" );
    throws_ok( sub{ MTToolbox::MTParams::_check_min_overlap($flash_params_href) },
              qr/FLASH param 'm' must be > 0/,
              "_check_min_overlap(m => -1) - throws FLASH param 'm' must be > 0" );
    eval{ MTToolbox::MTParams::_check_min_overlap($flash_params_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, '-1', "TooSmall->value eval test" );
    is( $e->MIN, 0, "TooSmall->MIN eval test" );
	
	$flash_params_href = {m => 1};
	lives_ok( sub{ MTToolbox::MTParams::_check_min_overlap($flash_params_href) },
            "_check_min_overlap(m => 1) - lives" );
}

# test _check_max_overlap
{
	my $flash_params_href = {};
	
	lives_ok( sub{ MTToolbox::MTParams::_check_max_overlap($flash_params_href) },
            "_check_max_overlap(does not exist) - lives (it is optional)" );
	
    $flash_params_href = {M => 'a'};
	throws_ok( sub{ MTToolbox::MTParams::_check_max_overlap($flash_params_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_max_overlap(M => a) - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_max_overlap($flash_params_href) },
              qr/FLASH param 'M' is not a digit/,
              "_check_max_overlap(M => a) - throws FLASH param 'M' is not a digit" );
    eval{ MTToolbox::MTParams::_check_max_overlap($flash_params_href) };
    my $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test" );
    
	$flash_params_href = {M => -1};
	throws_ok( sub{ MTToolbox::MTParams::_check_max_overlap($flash_params_href) },
              'MyX::Generic::Digit::TooSmall',
              "_check_max_overlap(M => -1) - throws TooSmall" );
    throws_ok( sub{ MTToolbox::MTParams::_check_max_overlap($flash_params_href) },
              qr/FLASH param 'M' must be > 0/,
              "_check_max_overlap(M => -1) - throws FLASH param 'M' must be > 0" );
    eval{ MTToolbox::MTParams::_check_max_overlap($flash_params_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, '-1', "TooSmall->value eval test" );
    is( $e->MIN, 0, "TooSmall->MIN eval test" );
	
	$flash_params_href = {M => 1};
	lives_ok( sub{ MTToolbox::MTParams::_check_max_overlap($flash_params_href) },
            "_check_max_overlap(m => 1) - lives" );
}

# test _check_ratio
{
	my $flash_params_href = {};
	
	lives_ok( sub{ MTToolbox::MTParams::_check_ratio($flash_params_href) },
            "_check_ratio(does not exist) - lives (it is optional)" );
	
	$flash_params_href = {x => 'a'};
	throws_ok( sub{ MTToolbox::MTParams::_check_ratio($flash_params_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_ratio(x => a) - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_ratio($flash_params_href) },
              qr/FLASH param 'x' is not a digit/,
              "_check_ratio(x => a) - throws FLASH param 'x' is not a digit" );
    eval{ MTToolbox::MTParams::_check_ratio($flash_params_href) };
    my $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test" );
    
    # OOB == Out Of Bounds
    $flash_params_href = {x => '-1'};
	throws_ok( sub{ MTToolbox::MTParams::_check_ratio($flash_params_href) },
              'MyX::Generic::Digit::OOB',
              "_check_ratio(x => -1) - throws OOB" );
    throws_ok( sub{ MTToolbox::MTParams::_check_ratio($flash_params_href) },
              qr/FLASH param \'x\' must be > 0 and < 1/,
              "_check_ratio(x => -1) - throws FLASH param \'x\' must be > 0 and < 1" );
    eval{ MTToolbox::MTParams::_check_ratio($flash_params_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::OOB');
    is( $e->value, '-1', "OOB->value eval test" );
    is( $e->MIN, 0, "OOB->MIN eval test" );
    is( $e->MAX, 1, "OOB->MAX eval test" );
    
    $flash_params_href = {x => '2'};
	throws_ok( sub{ MTToolbox::MTParams::_check_ratio($flash_params_href) },
              'MyX::Generic::Digit::OOB',
              "_check_ratio(x => 2) - throws OOB" );
    throws_ok( sub{ MTToolbox::MTParams::_check_ratio($flash_params_href) },
              qr/FLASH param \'x\' must be > 0 and < 1/,
              "_check_ratio(x => 2) - throws FLASH param \'x\' must be > 0 and < 1" );
    eval{ MTToolbox::MTParams::_check_ratio($flash_params_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::OOB');
    is( $e->value, '2', "OOB->value eval test" );
    is( $e->MIN, 0, "OOB->MIN eval test" );
    is( $e->MAX, 1, "OOB->MAX eval test" );
	
	$flash_params_href = {x => 1};
	lives_ok( sub{ MTToolbox::MTParams::_check_ratio($flash_params_href) },
            "_check_ratio(x => 1) - lives" );
}

# test _check_phred_offset
{
	my $flash_params_href = {};
	
	lives_ok( sub{ MTToolbox::MTParams::_check_phred_offset($flash_params_href) },
            "_check_phred_offset(does not exist) - lives (it is optional)" );
	
	$flash_params_href = {p => 'a'};
	throws_ok( sub{ MTToolbox::MTParams::_check_phred_offset($flash_params_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_phred_offset(p => 'a') - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_phred_offset($flash_params_href) },
              qr/FLASH param \'p\' is not a digit/,
              "_check_phred_offset(p => 'a') - throws FLASH param \'p\' is not a digit" );
    eval{ MTToolbox::MTParams::_check_phred_offset($flash_params_href) };
    my $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test" );
    
    $flash_params_href = {p => -1};
	throws_ok( sub{ MTToolbox::MTParams::_check_phred_offset($flash_params_href) },
              'MyX::Generic::Digit::TooSmall',
              "_check_phred_offset(p => -1) - throws TooSmall" );
    throws_ok( sub{ MTToolbox::MTParams::_check_phred_offset($flash_params_href) },
              qr/FLASH param \'p\' must be > 0/,
              "_check_phred_offset(p => -1) - throws FLASH param \'p\' must be > 0" );
    eval{ MTToolbox::MTParams::_check_phred_offset($flash_params_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, -1, "TooSmall->value eval test" );
    is( $e->MIN, 0, "TooSmall->MIN eval test" );

	$flash_params_href = {p => 1};
	lives_ok( sub{ MTToolbox::MTParams::_check_phred_offset($flash_params_href) },
            "_check_phred_offset(p => 1) - lives" );
}

# test _check_read_len
{
	my $flash_params_href = {};
	
	lives_ok( sub{ MTToolbox::MTParams::_check_read_len($flash_params_href) },
            "_check_read_len(does not exist) - lives (it is optional)" );
	
	$flash_params_href = {r => 'a'};
	throws_ok( sub{ MTToolbox::MTParams::_check_read_len($flash_params_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_read_len(r => 'a') - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_read_len($flash_params_href) },
              qr/FLASH param \'r\' is not a digit/,
              "_check_read_len(r => 'a') - throws FLASH param \'r\' is not a digit" );
    eval{ MTToolbox::MTParams::_check_read_len($flash_params_href) };
    my $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test" );
    
    $flash_params_href = {r => -1};
	throws_ok( sub{ MTToolbox::MTParams::_check_read_len($flash_params_href) },
              'MyX::Generic::Digit::TooSmall',
              "_check_read_len(r => -1) - throws TooSmall" );
    throws_ok( sub{ MTToolbox::MTParams::_check_read_len($flash_params_href) },
              qr/FLASH param \'r\' must be > 0/,
              "_check_read_len(r => -1) - throws FLASH param \'r\' must be > 0" );
    eval{ MTToolbox::MTParams::_check_read_len($flash_params_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, -1, "TooSmall->value eval test" );
    is( $e->MIN, 0, "TooSmall->MIN eval test" );

	$flash_params_href = {r => 1};
	lives_ok( sub{ MTToolbox::MTParams::_check_read_len($flash_params_href) },
            "_check_read_len(p => 1) - lives" );
}

# test _check_frag_len
{
	my $flash_params_href = {};
	
	lives_ok( sub{ MTToolbox::MTParams::_check_frag_len($flash_params_href) },
            "_check_frag_len(does not exist) - lives (it is optional)" );
	
	$flash_params_href = {f => 'a'};
	throws_ok( sub{ MTToolbox::MTParams::_check_frag_len($flash_params_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_frag_len(f => 'a') - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_frag_len($flash_params_href) },
              qr/FLASH param \'f\' is not a digit/,
              "_check_frag_len(f => 'a') - throws FLASH param \'f\' is not a digit" );
    eval{ MTToolbox::MTParams::_check_frag_len($flash_params_href) };
    my $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test" );
    
    $flash_params_href = {f => -1};
	throws_ok( sub{ MTToolbox::MTParams::_check_frag_len($flash_params_href) },
              'MyX::Generic::Digit::TooSmall',
              "_check_frag_len(f => -1) - throws TooSmall" );
    throws_ok( sub{ MTToolbox::MTParams::_check_frag_len($flash_params_href) },
              qr/FLASH param \'f\' must be > 0/,
              "_check_frag_len(f => -1) - throws FLASH param \'f\' must be > 0" );
    eval{ MTToolbox::MTParams::_check_frag_len($flash_params_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, -1, "TooSmall->value eval test" );
    is( $e->MIN, 0, "TooSmall->MIN eval test" );

	$flash_params_href = {f => 1};
	lives_ok( sub{ MTToolbox::MTParams::_check_frag_len($flash_params_href) },
            "_check_frag_len(p => 1) - lives" );
}

# test _check_sd
{
	my $flash_params_href = {};
	
	lives_ok( sub{ MTToolbox::MTParams::_check_sd($flash_params_href) },
            "_check_sd(does not exist) - lives (it is optional)" );
	
	$flash_params_href = {s => 'a'};
	throws_ok( sub{ MTToolbox::MTParams::_check_sd($flash_params_href) },
              'MyX::Generic::Digit::MustBeDigit',
              "_check_sd(s => 'a') - throws MustBeDigit" );
    throws_ok( sub{ MTToolbox::MTParams::_check_sd($flash_params_href) },
              qr/FLASH param \'s\' is not a digit/,
              "_check_sd(s => 'a') - throws FLASH param \'s\' is not a digit" );
    eval{ MTToolbox::MTParams::_check_sd($flash_params_href) };
    my $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test" );
    
    $flash_params_href = {s => -1};
	throws_ok( sub{ MTToolbox::MTParams::_check_sd($flash_params_href) },
              'MyX::Generic::Digit::TooSmall',
              "_check_sd(s => -1) - throws TooSmall" );
    throws_ok( sub{ MTToolbox::MTParams::_check_sd($flash_params_href) },
              qr/FLASH param \'s\' must be > 0/,
              "_check_sd(s => -1) - throws FLASH param \'s\' must be > 0" );
    eval{ MTToolbox::MTParams::_check_sd($flash_params_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, -1, "TooSmall->value eval test" );
    is( $e->MIN, 0, "TooSmall->MIN eval test" );

	$flash_params_href = {s => 1};
	lives_ok( sub{ MTToolbox::MTParams::_check_sd($flash_params_href) },
            "_check_sd(p => 1) - lives" );
}

# test _check_trim_to_base
{
	my $qc_params_href = {};
	
	lives_ok( sub{ MTToolbox::MTParams::_check_trim_to_base($qc_params_href) },
			 "_check_trim_to_base(does not exist) - lives (it is optional)" );
	
	$qc_params_href = {con_trim_to_base => 'a'};
    throws_ok( sub{ MTToolbox::MTParams::_check_trim_to_base($qc_params_href)},
              'MyX::Generic::Digit::MustBeDigit',
              "_check_trim_to_base - throws MustBeDigit");
	throws_ok( sub{ MTToolbox::MTParams::_check_trim_to_base($qc_params_href)},
               qr/Con filter param trim_to_base is not a digit/,
               "_check_trim_to_base - throw - not a digit - con" );
    eval{ MTToolbox::MTParams::_check_trim_to_base($qc_params_href) };
    my $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test");
	
	$qc_params_href = {con_trim_to_base => '-1'};
	throws_ok( sub{ MTToolbox::MTParams::_check_trim_to_base($qc_params_href)},
              'MyX::Generic::Digit::TooSmall',
              "_check_trim_to_base(-1) - throws TooSamll");
	throws_ok( sub{ MTToolbox::MTParams::_check_trim_to_base($qc_params_href)},
               qr/Con filter param trim_to_base must be > 0/,
               "_check_trim_to_base(-1) - throw - Con filter param trim_to_base must be > 0" );
    eval{ MTToolbox::MTParams::_check_trim_to_base($qc_params_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, -1, "TooSamll->value eval test" );
    is( $e->MIN, 0, "TooSmall->MIN eval test" );
	
	$qc_params_href = {con_trim_to_base => '100'};
	lives_ok( sub{ MTToolbox::MTParams::_check_trim_to_base($qc_params_href)},
		"_check_trim_to_base - lives - con" );
	
	$qc_params_href = {con_trim_to_base => 'NA'};
	lives_ok( sub{ MTToolbox::MTParams::_check_trim_to_base($qc_params_href)},
		"_check_trim_to_base(NA) - lives - con" );
	
	$qc_params_href = {SRC_trim_to_base => 'a'};
    throws_ok( sub{ MTToolbox::MTParams::_check_trim_to_base($qc_params_href)},
              'MyX::Generic::Digit::MustBeDigit',
              "_check_trim_to_base - throws MustBeDigit");
	throws_ok( sub{ MTToolbox::MTParams::_check_trim_to_base($qc_params_href)},
               qr/SRC filter param trim_to_base is not a digit/,
               "_check_trim_to_base - throws SRC filter param trim_to_base is not a digit" );
    eval{ MTToolbox::MTParams::_check_trim_to_base($qc_params_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test");
	
	$qc_params_href = {SRC_trim_to_base => '-1'};
	throws_ok( sub{ MTToolbox::MTParams::_check_trim_to_base($qc_params_href)},
              'MyX::Generic::Digit::TooSmall',
              "_check_trim_to_base(-1) - throws TooSamll");
	throws_ok( sub{ MTToolbox::MTParams::_check_trim_to_base($qc_params_href)},
               qr/SRC filter param trim_to_base must be > 0/,
               "_check_trim_to_base(-1) - throw - Con filter param trim_to_base must be > 0" );
    eval{ MTToolbox::MTParams::_check_trim_to_base($qc_params_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, -1, "TooSamll->value eval test" );
    is( $e->MIN, 0, "TooSmall->MIN eval test" );
	
	$qc_params_href = {SRC_trim_to_base => '100'};
	lives_ok( sub{ MTToolbox::MTParams::_check_trim_to_base($qc_params_href)},
		"_check_trim_to_base - lives - SRC" );
	
	$qc_params_href = {SRC_trim_to_base => 'na'};
	lives_ok( sub{ MTToolbox::MTParams::_check_trim_to_base($qc_params_href)},
		"_check_trim_to_base(na) - lives - SRC" );
}

# test _check_min_len
{
	my $qc_params_href = {};
	
	lives_ok( sub{ MTToolbox::MTParams::_check_min_len($qc_params_href) },
			 "_check_min_len(does not exist) - lives (it is optional)" );
	
	$qc_params_href = {con_min_len => 'a'};
    throws_ok( sub{ MTToolbox::MTParams::_check_min_len($qc_params_href)},
              'MyX::Generic::Digit::MustBeDigit',
              "_check_min_len - throws MustBeDigit");
	throws_ok( sub{ MTToolbox::MTParams::_check_min_len($qc_params_href)},
               qr/Con filter param min_len is not a digit/,
               "_check_min_len - throw - not a digit - con" );
    eval{ MTToolbox::MTParams::_check_min_len($qc_params_href) };
    my $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test");
	
	$qc_params_href = {con_min_len => '-1'};
	throws_ok( sub{ MTToolbox::MTParams::_check_min_len($qc_params_href)},
              'MyX::Generic::Digit::TooSmall',
              "_check_min_len(-1) - throws TooSamll");
	throws_ok( sub{ MTToolbox::MTParams::_check_min_len($qc_params_href)},
               qr/Con filter param min_len must be > 0/,
               "_check_min_len(-1) - throw - Con filter param min_len must be > 0" );
    eval{ MTToolbox::MTParams::_check_min_len($qc_params_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, -1, "TooSamll->value eval test" );
    is( $e->MIN, 0, "TooSmall->MIN eval test" );
	
	$qc_params_href = {con_min_len => '100'};
	lives_ok( sub{ MTToolbox::MTParams::_check_min_len($qc_params_href)},
		"_check_min_len - lives - con" );
	
	$qc_params_href = {con_min_len => 'NA'};
	lives_ok( sub{ MTToolbox::MTParams::_check_min_len($qc_params_href)},
		"_check_min_len(NA) - lives - con" );
	
	$qc_params_href = {SRC_min_len => 'a'};
    throws_ok( sub{ MTToolbox::MTParams::_check_min_len($qc_params_href)},
              'MyX::Generic::Digit::MustBeDigit',
              "_check_min_len - throws MustBeDigit");
	throws_ok( sub{ MTToolbox::MTParams::_check_min_len($qc_params_href)},
               qr/SRC filter param min_len is not a digit/,
               "_check_min_len - throws SRC filter param min_len is not a digit" );
    eval{ MTToolbox::MTParams::_check_min_len($qc_params_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test");
	
	$qc_params_href = {SRC_min_len => '-1'};
	throws_ok( sub{ MTToolbox::MTParams::_check_min_len($qc_params_href)},
              'MyX::Generic::Digit::TooSmall',
              "_check_min_len(-1) - throws TooSamll");
	throws_ok( sub{ MTToolbox::MTParams::_check_min_len($qc_params_href)},
               qr/SRC filter param min_len must be > 0/,
               "_check_min_len(-1) - throw - Con filter param min_len must be > 0" );
    eval{ MTToolbox::MTParams::_check_min_len($qc_params_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, -1, "TooSamll->value eval test" );
    is( $e->MIN, 0, "TooSmall->MIN eval test" );
	
	$qc_params_href = {SRC_min_len => '100'};
	lives_ok( sub{ MTToolbox::MTParams::_check_min_len($qc_params_href)},
		"_check_min_len - lives - SRC" );
	
	$qc_params_href = {SRC_min_len => 'na'};
	lives_ok( sub{ MTToolbox::MTParams::_check_min_len($qc_params_href)},
		"_check_min_len(na) - lives - SRC" );
}


# test _check_min_avg_qual
{
	my $qc_params_href = {};
	
	lives_ok( sub{ MTToolbox::MTParams::_check_min_avg_qual($qc_params_href) },
			 "_check_min_avg_qual(does not exist) - lives (it is optional)" );
	
	$qc_params_href = {con_min_avg_qual => 'a'};
    throws_ok( sub{ MTToolbox::MTParams::_check_min_avg_qual($qc_params_href)},
              'MyX::Generic::Digit::MustBeDigit',
              "_check_min_avg_qual - throws MustBeDigit");
	throws_ok( sub{ MTToolbox::MTParams::_check_min_avg_qual($qc_params_href)},
               qr/Con filter param min_avg_qual is not a digit/,
               "_check_min_avg_qual - throw - not a digit - con" );
    eval{ MTToolbox::MTParams::_check_min_avg_qual($qc_params_href) };
    my $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test");
	
	$qc_params_href = {con_min_avg_qual => '-1'};
	throws_ok( sub{ MTToolbox::MTParams::_check_min_avg_qual($qc_params_href)},
              'MyX::Generic::Digit::TooSmall',
              "_check_min_avg_qual(-1) - throws TooSamll");
	throws_ok( sub{ MTToolbox::MTParams::_check_min_avg_qual($qc_params_href)},
               qr/Con filter param min_avg_qual must be > 0/,
               "_check_min_avg_qual(-1) - throw - Con filter param min_avg_qual must be > 0" );
    eval{ MTToolbox::MTParams::_check_min_avg_qual($qc_params_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, -1, "TooSamll->value eval test" );
    is( $e->MIN, 0, "TooSmall->MIN eval test" );
	
	$qc_params_href = {con_min_avg_qual => '40'};
	lives_ok( sub{ MTToolbox::MTParams::_check_min_avg_qual($qc_params_href)},
		"_check_min_avg_qual - lives - con" );
	
	$qc_params_href = {SRC_min_avg_qual => 'a'};
    throws_ok( sub{ MTToolbox::MTParams::_check_min_avg_qual($qc_params_href)},
              'MyX::Generic::Digit::MustBeDigit',
              "_check_min_avg_qual - throws MustBeDigit");
	throws_ok( sub{ MTToolbox::MTParams::_check_min_avg_qual($qc_params_href)},
               qr/SRC filter param min_avg_qual is not a digit/,
               "_check_min_avg_qual - throws SRC filter param min_avg_qual is not a digit" );
    eval{ MTToolbox::MTParams::_check_min_avg_qual($qc_params_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test");
	
	$qc_params_href = {SRC_min_avg_qual => '-1'};
	throws_ok( sub{ MTToolbox::MTParams::_check_min_avg_qual($qc_params_href)},
              'MyX::Generic::Digit::TooSmall',
              "_check_min_avg_qual(-1) - throws TooSamll");
	throws_ok( sub{ MTToolbox::MTParams::_check_min_avg_qual($qc_params_href)},
               qr/SRC filter param min_avg_qual must be > 0/,
               "_check_min_avg_qual(-1) - throw - Con filter param min_avg_qual must be > 0" );
    eval{ MTToolbox::MTParams::_check_min_avg_qual($qc_params_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, -1, "TooSamll->value eval test" );
    is( $e->MIN, 0, "TooSmall->MIN eval test" );
	
	$qc_params_href = {SRC_min_avg_qual => '40'};
	lives_ok( sub{ MTToolbox::MTParams::_check_min_avg_qual($qc_params_href)},
		"_check_min_avg_qual - lives - SRC" );
}

# test _check_allow_gaps
{
	my $qc_params_href = {};
	
	lives_ok( sub{ MTToolbox::MTParams::_check_allow_gaps($qc_params_href) },
			 "_check_allow_gaps(does not exist) - lives (it is optional)" );
	
	$qc_params_href = {con_allow_gaps => 'a'};
    throws_ok( sub{ MTToolbox::MTParams::_check_allow_gaps($qc_params_href)},
               'MyX::Generic::BadValue',
               "_check_allow_gaps(a) - throws BadValue" );
	throws_ok( sub{ MTToolbox::MTParams::_check_allow_gaps($qc_params_href)},
               qr/Con allow_gaps must be Y or N/,
               "_check_allow_gaps(a) - throws Con allow_gaps must be Y or N" );
    eval{ MTToolbox::MTParams::_check_allow_gaps($qc_params_href) };
    my $e = Exception::Class->caught('MyX::Generic::BadValue');
    is( $e->value, 'a', "BadValue->value eval test" );
	
	$qc_params_href = {con_allow_gaps => '2'};
    throws_ok( sub{ MTToolbox::MTParams::_check_allow_gaps($qc_params_href)},
               'MyX::Generic::BadValue',
               "_check_allow_gaps(2) - throws BadValue" );
	throws_ok( sub{ MTToolbox::MTParams::_check_allow_gaps($qc_params_href)},
               qr/Con allow_gaps must be Y or N/,
               "_check_allow_gaps(2) - throws Con allow_gaps must be Y or N" );
    eval{ MTToolbox::MTParams::_check_allow_gaps($qc_params_href) };
    $e = Exception::Class->caught('MyX::Generic::BadValue');
    is( $e->value, 2, "BadValue->value eval test" );
	
	$qc_params_href = {con_min_avg_qual => 'Y'};
	lives_ok( sub{ MTToolbox::MTParams::_check_allow_gaps($qc_params_href)},
		"_check_allow_gaps - lives - con" );
	
	$qc_params_href = {con_min_avg_qual => 'y'};
	lives_ok( sub{ MTToolbox::MTParams::_check_allow_gaps($qc_params_href)},
		"_check_allow_gaps - lives - con" );
	
	$qc_params_href = {SRC_allow_gaps => 'a'};
	throws_ok( sub{ MTToolbox::MTParams::_check_allow_gaps($qc_params_href)},
               'MyX::Generic::BadValue',
               "_check_allow_gaps(a) - throws BadValue" );
	throws_ok( sub{ MTToolbox::MTParams::_check_allow_gaps($qc_params_href)},
               qr/SRC allow_gaps must be Y or N/,
               "_check_allow_gaps(a) - throws SRC allow_gaps must be Y or N" );
    eval{ MTToolbox::MTParams::_check_allow_gaps($qc_params_href) };
    $e = Exception::Class->caught('MyX::Generic::BadValue');
    is( $e->value, 'a', "BadValue->value eval test" );
	
	$qc_params_href = {SRC_allow_gaps => '2'};
	throws_ok( sub{ MTToolbox::MTParams::_check_allow_gaps($qc_params_href)},
               'MyX::Generic::BadValue',
               "_check_allow_gaps(2) - throws BadValue" );
	throws_ok( sub{ MTToolbox::MTParams::_check_allow_gaps($qc_params_href)},
               qr/SRC allow_gaps must be Y or N/,
               "_check_allow_gaps(2) - throws SRC allow_gaps must be Y or N" );
    eval{ MTToolbox::MTParams::_check_allow_gaps($qc_params_href) };
    $e = Exception::Class->caught('MyX::Generic::BadValue');
    is( $e->value, 2, "BadValue->value eval test" );
	
	$qc_params_href = {SRC_allow_gaps => 'N'};
	lives_ok( sub{ MTToolbox::MTParams::_check_allow_gaps($qc_params_href)},
		"_check_allow_gaps - lives - SRC" );
	
	$qc_params_href = {SRC_allow_gaps => 'n'};
	lives_ok( sub{ MTToolbox::MTParams::_check_allow_gaps($qc_params_href)},
		"_check_allow_gaps - lives - SRC" );
}


# test _check_allow_ambig
{
	my $qc_params_href = {};
	
	lives_ok( sub{ MTToolbox::MTParams::_check_allow_ambig($qc_params_href) },
			 "_check_allow_ambig(does not exist) - lives (it is optional)" );
	
	$qc_params_href = {con_allow_ambig => 'a'};
    throws_ok( sub{ MTToolbox::MTParams::_check_allow_ambig($qc_params_href)},
               'MyX::Generic::BadValue',
               "_check_allow_ambig(a) - throws BadValue" );
	throws_ok( sub{ MTToolbox::MTParams::_check_allow_ambig($qc_params_href)},
               qr/Con allow_ambig must be Y or N/,
               "_check_allow_ambig(a) - throws Con allow_ambig must be Y or N" );
    eval{ MTToolbox::MTParams::_check_allow_ambig($qc_params_href) };
    my $e = Exception::Class->caught('MyX::Generic::BadValue');
    is( $e->value, 'a', "BadValue->value eval test" );
	
	$qc_params_href = {con_allow_ambig => '2'};
    throws_ok( sub{ MTToolbox::MTParams::_check_allow_ambig($qc_params_href)},
               'MyX::Generic::BadValue',
               "_check_allow_ambig(2) - throws BadValue" );
	throws_ok( sub{ MTToolbox::MTParams::_check_allow_ambig($qc_params_href)},
               qr/Con allow_ambig must be Y or N/,
               "_check_allow_ambig(2) - throws Con allow_ambig must be Y or N" );
    eval{ MTToolbox::MTParams::_check_allow_ambig($qc_params_href) };
    $e = Exception::Class->caught('MyX::Generic::BadValue');
    is( $e->value, 2, "BadValue->value eval test" );
	
	$qc_params_href = {con_min_avg_qual => 'Y'};
	lives_ok( sub{ MTToolbox::MTParams::_check_allow_ambig($qc_params_href)},
		"_check_allow_ambig - lives - con" );
	
	$qc_params_href = {con_min_avg_qual => 'y'};
	lives_ok( sub{ MTToolbox::MTParams::_check_allow_ambig($qc_params_href)},
		"_check_allow_ambig - lives - con" );
	
	$qc_params_href = {SRC_allow_ambig => 'a'};
	throws_ok( sub{ MTToolbox::MTParams::_check_allow_ambig($qc_params_href)},
               'MyX::Generic::BadValue',
               "_check_allow_ambig(a) - throws BadValue" );
	throws_ok( sub{ MTToolbox::MTParams::_check_allow_ambig($qc_params_href)},
               qr/SRC allow_ambig must be Y or N/,
               "_check_allow_ambig(a) - throws SRC allow_ambig must be Y or N" );
    eval{ MTToolbox::MTParams::_check_allow_ambig($qc_params_href) };
    $e = Exception::Class->caught('MyX::Generic::BadValue');
    is( $e->value, 'a', "BadValue->value eval test" );
	
	$qc_params_href = {SRC_allow_ambig => '2'};
	throws_ok( sub{ MTToolbox::MTParams::_check_allow_ambig($qc_params_href)},
               'MyX::Generic::BadValue',
               "_check_allow_ambig(2) - throws BadValue" );
	throws_ok( sub{ MTToolbox::MTParams::_check_allow_ambig($qc_params_href)},
               qr/SRC allow_ambig must be Y or N/,
               "_check_allow_ambig(2) - throws SRC allow_ambig must be Y or N" );
    eval{ MTToolbox::MTParams::_check_allow_ambig($qc_params_href) };
    $e = Exception::Class->caught('MyX::Generic::BadValue');
    is( $e->value, 2, "BadValue->value eval test" );
	
	$qc_params_href = {SRC_allow_ambig => 'N'};
	lives_ok( sub{ MTToolbox::MTParams::_check_allow_ambig($qc_params_href)},
		"_check_allow_ambig - lives - SRC" );
	
	$qc_params_href = {SRC_allow_ambig => 'n'};
	lives_ok( sub{ MTToolbox::MTParams::_check_allow_ambig($qc_params_href)},
		"_check_allow_ambig - lives - SRC" );
}

# test _check_min_c_score
{
	my $qc_params_href = {};
	
	lives_ok( sub{ MTToolbox::MTParams::_check_min_c_score($qc_params_href) },
			 "_check_min_c_score(does not exist) - lives (it is optional)" );
	
	$qc_params_href = {con_min_c_score => 'a'};
    throws_ok( sub{ MTToolbox::MTParams::_check_min_c_score($qc_params_href)},
              'MyX::Generic::Digit::MustBeDigit',
              "_check_min_c_score - throws MustBeDigit");
	throws_ok( sub{ MTToolbox::MTParams::_check_min_c_score($qc_params_href)},
               qr/Con filter param min_c_score is not a digit/,
               "_check_min_c_score - throw - not a digit - con" );
    eval{ MTToolbox::MTParams::_check_min_c_score($qc_params_href) };
    my $e = Exception::Class->caught('MyX::Generic::Digit::MustBeDigit');
    is( $e->value, 'a', "MustBeDigit->value eval test");
	
	$qc_params_href = {con_min_c_score => '-1'};
	throws_ok( sub{ MTToolbox::MTParams::_check_min_c_score($qc_params_href)},
              'MyX::Generic::Digit::TooSmall',
              "_check_min_c_score(-1) - throws TooSamll");
	throws_ok( sub{ MTToolbox::MTParams::_check_min_c_score($qc_params_href)},
               qr/Con filter param min_c_score must be > 0/,
               "_check_min_c_score(-1) - throw - Con filter param min_c_score must be > 0" );
    eval{ MTToolbox::MTParams::_check_min_c_score($qc_params_href) };
    $e = Exception::Class->caught('MyX::Generic::Digit::TooSmall');
    is( $e->value, -1, "TooSamll->value eval test" );
    is( $e->MIN, 0, "TooSmall->MIN eval test" );
	
	$qc_params_href = {con_min_c_score => '40'};
	lives_ok( sub{ MTToolbox::MTParams::_check_min_c_score($qc_params_href)},
		"_check_min_c_score - lives - con" );
}

# test _check_merge_params_file
{
    # check error when given an non-existant file
    my $global_href = {merge_params_file => "", seq_type => 2};
    throws_ok( sub{ MTToolbox::MTParams::_check_merge_params_file($global_href) },
              'MyX::Generic::DoesNotExist::File',
              "_check_merge_params_file() - throws DoesNotExist::File" );
    throws_ok( sub{ MTToolbox::MTParams::_check_merge_params_file($global_href) },
              qr/Merge params file does not exist/,
              "_check_merge_params_file() - throws Merge params file does not exist" );
    eval{ MTToolbox::MTParams::_check_merge_params_file($global_href) };
    my $e = Exception::Class->caught('MyX::Generic::DoesNotExist::File');
    is( $e->file_name, "", "DoesNotExist::File->file_name eval test" );
    
    # check error when the given file is empty
    ($fh, $filename) = tempfile();
    $global_href = {merge_params_file => $filename, seq_type => 2};
    throws_ok( sub{ MTToolbox::MTParams::_check_merge_params_file($global_href) },
              'MyX::Generic::File::Empty',
              "_check_merge_params_file() - throws File::Empty" );
    throws_ok( sub{ MTToolbox::MTParams::_check_merge_params_file($global_href) },
              qr/Merge params file is empty/,
              "_check_merge_params_file() - throws Merge params file is empty" );
    eval{ MTToolbox::MTParams::_check_merge_params_file($global_href) };
    $e = Exception::Class->caught('MyX::Generic::File::Empty');
    is( $e->file_name, $filename, "File::Empty->file_name eval test" );
    
    # check a good params file
    print $fh "My params file";
    close($fh);
    is( MTToolbox::MTParams::_check_merge_params_file($global_href), 1,
       "_check_merge_params_file - good" );
}

# test _check_samples
{
	# check when passing an array ref with no entries
	my @tmp_arr = ();
	my $global_href = {sample => \@tmp_arr};
	throws_ok( sub{ MTToolbox::MTParams::_check_samples($global_href) },
              'MTToolbox::MyX::MTParams::MissingTagValue',
              "_check_samples(empty aref) - throws MissingTagValue" );
	
	# Create a sample that will pass
	my $sample_href = {};
	# now create a working example with seq_type == 1 (SE)
	my $split_samples_root_dir = tempdir( CLEANUP => 1 );
	($fh, $filename) = tempfile( DIR => $split_samples_root_dir,
								 SUFFIX => ".fastq" );
	print $fh ">a fastq file line";
	close($fh);
	my $filename_no_path = fileparse($filename);
	$sample_href->{fwd_file} = $filename_no_path;
	$sample_href->{sample_id} = "P0";
	my @sample_arr = ($sample_href);
	
	# add all previously created and other needed variables
	$global_href->{sample} = \@sample_arr;
	$global_href->{split_samples_root_dir} = $split_samples_root_dir;
	$global_href->{seq_type} = 1;
	
	# test for correct execution of _check_samples when passing an array ref
	is( MTToolbox::MTParams::_check_samples($global_href), 1,
	   "_check_samples(array ref) -- should work" );
	
	# check when passing a hash element with just one sample
	$global_href->{sample} = $sample_href;
	is( MTToolbox::MTParams::_check_samples($global_href), 1,
	   "_check_samples(hash ref) -- should work" );
	is( ref $global_href->{sample}, "ARRAY",
	   "_check_samples(href) -- change to aref" );
	
}

# test check_sample_entry
{
    # test the required tags
    my $sample_href = {};
    throws_ok( sub{ MTToolbox::MTParams::check_sample_entry($sample_href, undef, undef) },
              'MTToolbox::MyX::MTParams::MissingRequiredTag',
              "_check_smaple_entry({}, {}) - throws MissingRequiredTag" );
    
    # test die because of unknown tag
    $sample_href = {sample_id => "P1", barcode => "AAAAA", unknown => 1};
    throws_ok( sub{ MTToolbox::MTParams::check_sample_entry($sample_href, undef, undef) },
              'MTToolbox::MyX::MTParams::UnrecognizedTag',
              "_check_smaple_entry({}, {unknown}) - throws UnrecognizedTag" );
    delete $sample_href->{unknown};
    
    # tests for seq_type and file type agreement
    {
        # the seq_type is defined in the sample INCORRECTLY
        $sample_href->{seq_type} = "A";
        throws_ok( sub{ MTToolbox::MTParams::check_sample_entry($sample_href, undef, undef) },
                  'MTToolbox::MyX::MTParams::BadSeqType',
                  "_check_smaple_entry({}, {seq_type = A}) - throws BadSeqType" );
        $sample_href->{seq_type} = "4";
        throws_ok( sub{ MTToolbox::MTParams::check_sample_entry($sample_href, undef, undef) },
                  'MTToolbox::MyX::MTParams::BadSeqType',
                  "_check_smaple_entry({}, {seq_type = 4}) - throws BadSeqType" );
        
        ## the seq_type is defined in the sample CORRECTLY as SE reads
        # but the fwd_file tag is missing
        $sample_href->{seq_type} = 1;
        throws_ok( sub{ MTToolbox::MTParams::check_sample_entry($sample_href, undef, undef) },
                  'MTToolbox::MyX::MTParams::MissingRequiredTag',
                  "check_sample_entry() - no fwd_file - throws MissingRequiredTag" );
        
        # now create a fwd_file but it does not exist in the split_samples_root_dir
        my $split_samples_root_dir = tempdir( CLEANUP => 1 );
        ($fh, $filename) = tempfile();
        my $filename_no_path = fileparse($filename);
        $sample_href->{fwd_file} = $filename_no_path;
        throws_ok( sub{ MTToolbox::MTParams::check_sample_entry($sample_href,
                                                      undef,
                                                      $split_samples_root_dir) },
                  'MyX::Generic::DoesNotExist::File',
                  "check_sample_entry() - throws DoesNotExist::File" );
        
        # now create a fwd_file but it is empty in the split_samples_root_dir
        ($fh, $filename) = tempfile( DIR => $split_samples_root_dir );
        $filename_no_path = fileparse($filename);
        $sample_href->{fwd_file} = $filename_no_path;
        throws_ok( sub{ MTToolbox::MTParams::check_sample_entry($sample_href,
                                                      undef,
                                                      $split_samples_root_dir) },
                  'MyX::Generic::File::Empty',
                  "check_sample_entry() - File::Empty" );
        
        # now create a fwd_file with something in it, but a wrong extension
        ($fh, $filename) = tempfile( DIR => $split_samples_root_dir,
                                     SUFFIX => ".blah" );
        print $fh ">a fastq file line";
        close($fh);
        $filename_no_path = fileparse($filename);
        $sample_href->{fwd_file} = $filename_no_path;
        throws_ok( sub{ MTToolbox::MTParams::check_sample_entry($sample_href,
                                                      undef,
                                                      $split_samples_root_dir) },
                  'MyX::Generic::File::BadExtension',
                  "check_sample_entry() - fwd_file has bad suffix" );
        
        # now create a working example with seq_type == 1 (SE)
        ($fh, $filename) = tempfile( DIR => $split_samples_root_dir,
                                     SUFFIX => ".fastq" );
        print $fh ">a fastq file line";
        close($fh);
        $filename_no_path = fileparse($filename);
        $sample_href->{fwd_file} = $filename_no_path;
        is( MTToolbox::MTParams::check_sample_entry($sample_href,
                                          undef,
                                          $split_samples_root_dir), 1,
                  "check_sample_entry() - fwd_file should work!" );
        
        # now test a globally defined seq_type
        delete $sample_href->{seq_type};
        #$global_href->{seq_type} = 1;
        ($fh, $filename) = tempfile( DIR => $split_samples_root_dir,
                                     SUFFIX => ".fastq" );
        print $fh ">a fastq file line";
        close($fh);
        $filename_no_path = fileparse($filename);
        $sample_href->{fwd_file} = $filename_no_path;
        is( MTToolbox::MTParams::check_sample_entry($sample_href,
                                          1,
                                          $split_samples_root_dir), 1,
                  "check_sample_entry() - fwd_file should work w/ global seq_type!" );
        
        ## Now test when seq_type is CORRECTLY set as PE w/ overlap (2)
        delete $global_href->{seq_type};
        $sample_href->{seq_type} = 2;
        delete $sample_href->{fwd_file};
        
        # test that seq_type == 2 works when providing a merged file
        ($fh, $filename) = tempfile( DIR => $split_samples_root_dir,
                                     SUFFIX => ".fastq" );
        print $fh ">a fastq merged file line";
        close($fh);
        $filename_no_path = fileparse($filename);
        $sample_href->{merged_file} = $filename_no_path;
        is( MTToolbox::MTParams::check_sample_entry($sample_href,
                                          undef,
                                          $split_samples_root_dir), 1,
                  "check_sample_entry() - merged_file should work!" );
        
        # dies when no merged file or fwd file
        delete $sample_href->{merged_file};
        throws_ok( sub{ MTToolbox::MTParams::check_sample_entry($sample_href,
                                                      undef,
                                                      $split_samples_root_dir) },
                  'MTToolbox::MyX::MTParams::MissingRequiredTag',
                  "check_sample_entry(seq_type = 2) - no fwd_file - throws MissingRequiredTag" );
        
        # dies when only a fwd file is provided
        ($fh, $filename) = tempfile();
        $sample_href->{fwd_file} = $filename;
        throws_ok( sub{ MTToolbox::MTParams::check_sample_entry($sample_href,
                                                      undef,
                                                      $split_samples_root_dir) },
                  'MTToolbox::MyX::MTParams::MissingRequiredTag',
                  "check_sample_entry(seq_type = 2) - only fwd_file - throws MissingRequiredTag" );
        
        # dies when only a rev file is provided
        $sample_href->{rev_file} = $sample_href->{fwd_file};
        delete $sample_href->{fwd_file};
        throws_ok( sub{ MTToolbox::MTParams::check_sample_entry($sample_href,
                                                      undef,
                                                      $split_samples_root_dir) },
                  'MTToolbox::MyX::MTParams::MissingRequiredTag',
                  "check_sample_entry(seq_type = 2) - only rev_file - throws MissingRequiredTag" );
        
        # lives when good fwd and rev files are provided
        ($fh, $filename) = tempfile( DIR => $split_samples_root_dir,
                                     SUFFIX => ".fastq" );
        print $fh ">a fastq fwd file line";
        close($fh);
        $filename_no_path = fileparse($filename);
        $sample_href->{fwd_file} = $filename_no_path;
        ($fh, $filename) = tempfile( DIR => $split_samples_root_dir,
                                     SUFFIX => ".fastq" );
        print $fh ">a fastq fwd file line";
        close($fh);
        $filename_no_path = fileparse($filename);
        $sample_href->{rev_file} = $filename_no_path;
        is( MTToolbox::MTParams::check_sample_entry($sample_href,
                                          undef,
                                          $split_samples_root_dir), 1,
                  "check_sample_entry(seq_type = 2) - should work!" );
        
        ## Now test when seq_type is CORRECTLY set as PE w/OUT overlap (3)
        delete $global_href->{seq_type};
        $sample_href->{seq_type} = 3;
        delete $sample_href->{fwd_file};
        delete $sample_href->{rev_file};
        
        # dies when no fwd or rev file are provided
        throws_ok( sub{ MTToolbox::MTParams::check_sample_entry($sample_href,
                                                      undef,
                                                      $split_samples_root_dir) },
                  'MTToolbox::MyX::MTParams::MissingRequiredTag',
                  "check_sample_entry(seq_type = 3) - no fwd or rev - throws MissingRequiredTag" );
        
        # dies when only a fwd file is provided
        ($fh, $filename) = tempfile();
        $sample_href->{fwd_file} = $filename;
        throws_ok( sub{ MTToolbox::MTParams::check_sample_entry($sample_href,
                                                      undef,
                                                      $split_samples_root_dir) },
                  'MTToolbox::MyX::MTParams::MissingRequiredTag',
                  "check_sample_entry(seq_type = 3) - only fwd_file - throws MissingRequiredTag" );
        
        # dies when only a rev file is provided
        $sample_href->{rev_file} = $sample_href->{fwd_file};
        delete $sample_href->{fwd_file};
        throws_ok( sub{ MTToolbox::MTParams::check_sample_entry($sample_href,
                                                      undef,
                                                      $split_samples_root_dir) },
                  'MTToolbox::MyX::MTParams::MissingRequiredTag',
                  "check_sample_entry(seq_type = 3) - only rev_file - throws MissingRequiredTag" );
        
        # lives when good fwd and rev files are provided
        ($fh, $filename) = tempfile( DIR => $split_samples_root_dir,
                                     SUFFIX => ".fastq" );
        print $fh ">a fastq fwd file line";
        close($fh);
        $filename_no_path = fileparse($filename);
        $sample_href->{fwd_file} = $filename_no_path;
        ($fh, $filename) = tempfile( DIR => $split_samples_root_dir,
                                     SUFFIX => ".fastq" );
        print $fh ">a fastq fwd file line";
        close($fh);
        $filename_no_path = fileparse($filename);
        $sample_href->{rev_file} = $filename_no_path;
        is( MTToolbox::MTParams::check_sample_entry($sample_href,
                                          undef,
                                          $split_samples_root_dir), 1,
                  "check_sample_entry(seq_type = 3) - should work!" );
    }
}





### Subroutines ###
sub get_test_href {
    my $tempdir = tempdir();
    
    # I need to add samples to this!!
    
    my $xml_href = {
        output_root => $tempdir,
		split_samples_root_dir => '',
        parallelize => 'Y',
        fwd_primer => "GTGCCAGC[AC]GCCGCGGTAA",
        rev_primer => "GGACTAC[ACT][ACG]GGGT[AT]TCTAAT",
        seq_type => 2,
		keep_tmp_files => 'N',
        fwd_linker => "GA",
        rev_linker => "AC",
        fwd_mer_len => 8,
        rev_mer_len => 5,
        fwd_max_shifts => 5,
        rev_max_shifts => 5,
        min_con_depth => 2,
        diginorm_max => "NA",
		con_algo => 'MSA',
        qc_params => {
            con_min_len => "NA",
            con_min_avg_qual => 0,
            con_allow_gaps => 0,
            con_allow_ambig => 0,
            SRC_min_len => "NA",
            SRC_min_avg_qual => 0,
            SRC_allow_gaps => 0,
            SRC_allow_ambig => 0,
        },
        filter_params => {
            blast_db => "",
            min_eval => 0.01,
            min_perc_iden => 94,
        },
        flash_params => {
            m => 30,
            M => 250,
            x => 0.25,
            p => 33,
            r => 250,
            f => 310,
            s => 20,
        },
        #sample => \@sample_arr,
    };
    
    return $xml_href;
}

sub get_test_file {
    my ($fh, $filename) = tempfile();
    
    my $xml_href = get_test_href();
    
    # write the hrep
    my $xml = XMLout($xml_href);
    print $fh $xml;
    
    close($fh);
    
    return $filename;
}
