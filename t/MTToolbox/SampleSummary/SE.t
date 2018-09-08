
use strict;
use warnings;

use MTToolbox::SampleSummary::SE;
use Test::More tests => 16;

BEGIN { use_ok( 'MTToolbox::SampleSummary::SE'); }

### Create a MTToolbox::SampleSummary::SE object to test
my $ssi = MTToolbox::SampleSummary::SE->new();

# set some values
is( $ssi->set_sample_name("sample1"), 1, "set_sample_name(sample1)" );
is( $ssi->set_sample_id("P1_ATCG"), 1, "set_sample_id(P1_ATCG)" );
is( $ssi->set_file_format('fastq'), 1, "set_file_format" );
is( $ssi->set_match_count(100), 1, "set_match_count(100)" );
is( $ssi->set_mismatch_count(10), 1, "set_mismatch_count(10)" );
is( $ssi->set_seq_count(110), 1, "set_seq_count(110)" );
is( $ssi->set_MT_count(25), 1, "set_MT_count(25)" );

is( $ssi->get_sample_name(), "sample1", "get_sample_name()" );
is( $ssi->get_sample_id(), "P1_ATCG", "get_sample_id()" );
is( $ssi->get_file_format(), "fastq", "get_file_format()" );
is( $ssi->get_match_count(), 100, "get_match_count()" );
is( $ssi->get_mismatch_count(), 10, "get_mismatch_count()" );
is( $ssi->get_seq_count(), 110, "get_seq_count()" );
is( $ssi->get_MT_count(), 25, "get_MT_count()" );

is( $ssi->to_string(),
   "Sample type: SE
Sample name: sample1
ID: P1_ATCG
File format: fastq
Seq count: 110
Match count: 100
Mismatch count: 10
Molecule Tag (MT) count: 25
Single read category count: --
", "to_string()"
);
