
use strict;
use warnings;

use MTToolbox::SampleSummary::PE::Overlap;
use File::Temp qw/ tempfile /;
use Test::More tests => 47;

BEGIN { use_ok( 'MTToolbox::SampleSummary::PE::Overlap'); }

### Create a MTToolbox::SampleSummary::PE::Overlap object to test
my $ssi = MTToolbox::SampleSummary::PE::Overlap->new();

# set some values -- these are also tested in the parent class
is( $ssi->set_sample_name("sample1"), 1, "set_sample_name(sample1)" );
is( $ssi->set_sample_id("P1_ATCG"), 1, "set_sample_id(P1_ATCG)" );
is( $ssi->set_file_format('fastq'), 1, "set_file_format" );
is( $ssi->set_match_count(100), 1, "set_match_count(100)" );
is( $ssi->set_mismatch_count(10), 1, "set_mismatch_count(10)" );
is( $ssi->set_seq_count(110), 1, "set_seq_count(110)" );
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
   "Sample type: PE w/ Overlap
Sample name: sample1
ID: P1_ATCG
File format: fastq
Seq count: 110
Merged count: --
Not merged count: --
Percent merged: --
Match count: 100
Mismatch count: 10
Molecule Tag (MT) count: 25
Single read category count: 20
", "to_string()"
);


# Add some merge information
is( $ssi->set_merged_count(90), 1, "set_merged_count(110)" );
is( $ssi->set_not_merged_count(20), 1, "set_not_merged_count(20)" );
is( $ssi->set_match_count(80), 1, "set_match_count(80)" );
is( $ssi->get_merged_count(), 90, "get_merged_count()" );
is( $ssi->get_not_merged_count(), 20, "get_not_merged_count()" );
is( $ssi->get_perc_merged(), "81%", "get_perc_merged()" );
is( $ssi->to_string(),
   "Sample type: PE w/ Overlap
Sample name: sample1
ID: P1_ATCG
File format: fastq
Seq count: 110
Merged count: 90
Not merged count: 20
Percent merged: 81%
Match count: 80
Mismatch count: 10
Molecule Tag (MT) count: 25
Single read category count: 20
", "to_string()"
);

# test the load_summary_file
{
    my ($fh, $file_name) = tempfile();
    my $ssi = MTToolbox::SampleSummary::PE::Overlap->new();
    
    # set some new values
    is( $ssi->set_sample_name("sample2"), 1, "set_sample_name(sample2)" );
    is( $ssi->set_sample_id("P2_ATCG"), 1, "set_sample_id(P2_ATCG)" );
    is( $ssi->set_file_format('fastq'), 1, "set_file_format" );
    is( $ssi->set_seq_count(220), 1, "set_seq_count(220)" );
    is( $ssi->set_merged_count(220), 1, "set_merged_count(220)" );
    is( $ssi->set_not_merged_count(0), 1, "set_not_merged_count(0)" );
    is( $ssi->set_match_count(200), 1, "set_match_count(200)" );
    is( $ssi->set_mismatch_count(20), 1, "set_mismatch_count(20)" );
    is( $ssi->set_MT_count(35), 1, "set_MT_count(30)" );
    is( $ssi->set_SRC_count(25), 1, "set_SRC_count(25)" );
    
    # print to temp file
    print $fh $ssi->to_string();
    
    # test
    is( $ssi->load_summary_file($file_name), 1, "load_summary_file" );
    is( $ssi->get_sample_name(), "sample2", "get_sample_name()" );
    is( $ssi->get_sample_id(), "P2_ATCG", "get_sample_id()" );
    is( $ssi->get_file_format(), "fastq", "get_file_format()" );
    is( $ssi->get_seq_count(), 220, "get_seq_count()" );
    is( $ssi->get_merged_count(), 220, "get_merged_count()" );
    is( $ssi->get_not_merged_count(), 0, "get_not_merged_count()" );
    is( $ssi->get_perc_merged(), "100%", "get_perc_merged()" );
    is( $ssi->get_match_count(), 200, "get_match_count()" );
    is( $ssi->get_mismatch_count(), 20, "get_mismatch_count()" );
    is( $ssi->get_MT_count(), 35, "get_MT_count()" );
    is( $ssi->get_SRC_count(), 25, "get_SRC_count()" );
}