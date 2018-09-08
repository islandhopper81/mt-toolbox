
use strict;
use warnings;

use MTToolbox::SampleSummary;
use File::Temp qw/ tempfile/;
use Test::More tests => 27;

BEGIN { use_ok( 'MTToolbox::SampleSummary'); }

### Create a MTToolbox::SampleSummary object to test
my $ssi = MTToolbox::SampleSummary->new();

# set some values
is( $ssi->set_sample_name('sample1'), 1, "set_sample_name(sample1)" );
is( $ssi->set_sample_id('P1_ATCG'), 1, "set_sample_id(P1_ATCG)" );
is( $ssi->set_file_format('fastq'), 1, "set_file_format" );
is( $ssi->set_seq_count(110), 1, "set_seq_count(110)" );
is( $ssi->set_match_count(100), 1, "set_match_count(100)" );
is( $ssi->set_mismatch_count(10), 1, "set_mismatch_count(10)" );
is( $ssi->set_MT_count(25), 1, "set_MT_count(25)" );
is( $ssi->set_SRC_count(20), 1, "set_SRC_count(20)" );

is( $ssi->get_sample_name(), "sample1", "get_sample_name()" );
is( $ssi->get_sample_id(), "P1_ATCG", "get_sample_id()" );
is( $ssi->get_file_format(), "fastq", "get_file_format()" );
is( $ssi->get_match_count(), 100, "get_match_count()" );
is( $ssi->get_mismatch_count(), 10, "get_mismatch_count()" );
is( $ssi->get_seq_count(), 110, "get_seq_count()" );
is( $ssi->get_MT_count(), 25, "get_MT_count()" );
is( $ssi->get_SRC_count(), 20, "get_SRC_count()" );

is( $ssi->to_string(),
"Sample name: sample1
ID: P1_ATCG
File format: fastq
Seq count: 110
Match count: 100
Mismatch count: 10
Molecule Tag (MT) count: 25
Single read category count: 20
", "to_string()"
);

# test the load_summary_file
{
    my ($fh, $file_name) = tempfile();
    print $fh, $ssi->to_string();
    is( $ssi->load_summary_file($file_name), 1, "load_summary_file" );
    is( $ssi->get_sample_name(), "sample1", "get_sample_name()" );
    is( $ssi->get_sample_id(), "P1_ATCG", "get_sample_id()" );
    is( $ssi->get_file_format(), "fastq", "get_file_format()" );
    is( $ssi->get_match_count(), 100, "get_match_count()" );
    is( $ssi->get_mismatch_count(), 10, "get_mismatch_count()" );
    is( $ssi->get_seq_count(), 110, "get_seq_count()" );
    is( $ssi->get_MT_count(), 25, "get_MT_count()" );
    is( $ssi->get_SRC_count(), 20, "get_SRC_count()" );
}

