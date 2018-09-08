
use strict;
use warnings;

use MTToolbox::RunSummary;
use MTToolbox::SampleSummary::SE;
use MTToolbox::SampleSummary::PE::NoOverlap;
use MTToolbox::SampleSummary::PE::Overlap;
use File::Temp qw/ tempfile tempdir /;
use Test::More tests => 57;

BEGIN { use_ok( 'MTToolbox::RunSummary'); }

### Create a MTToolbox::SampleSummary::SE object to test
my $summ = MTToolbox::RunSummary->new();

# test the simple setter methods
{
    is( $summ->set_sample_count(5), 1, "set_sample_count(5)" );
    is( $summ->set_total_seqs_count(100), 1, "set_total-seqs_count(100)" );
    is( $summ->set_merged_seqs_count(50), 1, "set_merged_seqs_count(50)" );
    is( $summ->set_match_count(30), 1, "set_match_count(30)" );
    is( $summ->set_mismatch_count(20), 1, "set_mismatch_count(20)" );
    is( $summ->set_mt_count(5), 1, "set_mt_count(5)" );
    is( $summ->set_SRC_count(10), 1, "set_SRC_count(10)" );
}

# test the simple getter methods
{
    is( $summ->get_sample_count(), 5, "get_sample_count() - 5" );
    is( $summ->get_total_seqs_count(), 100, "get_total-seqs_count() - 100" );
    is( $summ->get_merged_seqs_count(), 50, "get_merged_seqs_count() - 50" );
    is( $summ->get_not_merged_seqs_count(), 50, "get_not_merged_seqs_count() - 50" );
    is( $summ->get_merged_percent(), "50%", "get_merged_percent() - 50" );
    is( $summ->get_match_count(), 30, "get_match_count() - 30" );
    is( $summ->get_mismatch_count(), 20, "get_mismatch_count() - 20" );
    is( $summ->get_mt_count(), 5, "get_mt_count() - 5" );
    is( $summ->get_SRC_count(), 10, "get_SRC_count() - 10" );
}

# test set/get_samples_href
{
    # create two samples
    my $sample0_summ = MTToolbox::SampleSummary::PE::Overlap->new();
    $sample0_summ->set_sample_name('sample0');
    $sample0_summ->set_sample_id('P0_ATCG');
    $sample0_summ->set_file_format('fastq');
    $sample0_summ->set_seq_count(221120);
    $sample0_summ->set_merged_count(23816);
    $sample0_summ->set_not_merged_count(197304);
    $sample0_summ->set_match_count(18979);
    $sample0_summ->set_mismatch_count(4837);
    $sample0_summ->set_MT_count(5123);
    $sample0_summ->set_SRC_count(18979);
    
    my $sample1_summ = MTToolbox::SampleSummary::PE::Overlap->new();
    $sample1_summ->set_sample_name('sample1');
    $sample1_summ->set_sample_id('P1_ATCG');
    $sample1_summ->set_file_format('fastq');
    $sample1_summ->set_seq_count(147203);
    $sample1_summ->set_merged_count(9534);
    $sample1_summ->set_not_merged_count(137669);
    $sample1_summ->set_match_count(6822);
    $sample1_summ->set_mismatch_count(2712);
    $sample1_summ->set_MT_count(1718);
    $sample1_summ->set_SRC_count(6822);
    
    my %hash;
    $hash{P0_ATCG} = $sample0_summ;
    $hash{P1_ATCG} = $sample1_summ;
    
    is( $summ->set_samples_href(\%hash), 1, "set_sample_href()" );
    
    is_deeply( $summ->get_samples_href(), \%hash, "get_samples_href()" );
}

# test _get_sample_name helper method
{
    my $name = "output_20130207/P1_CGTCGGTA/summary.txt";
    is( MTToolbox::RunSummary::_get_sample_name($name), "P1_CGTCGGTA", "_get_sample_name()" );
}

# test load_sample_summaries with normal PE w/ Overlap output
{
    # reset summ -- I save this as the global so I can use it to test to_string
    $summ = MTToolbox::RunSummary->new();
    
    # make a directory structure similar to the output of MTToolbox
    my $root_dir = tempdir();
    mkdir "$root_dir/P0_ATCG";
    mkdir "$root_dir/P1_ATCG";
    
    my $file0 = "$root_dir/P0_ATCG/summary.txt";
    my $file1 = "$root_dir/P1_ATCG/summary.txt";
    
    my $summ0_str = "Sample type: PE w/ Overlap\n" .
                    "Sample name: sample0\n" .
                    "ID: P0_ATCG\n" . 
                    "File format: fastq\n" .
                    "Seq count: 221120\n" .
                    "Merged count: 23816\n" .
                    "Not merged count: 197304\n" .
                    "Percent merged: 10%\n" .
                    "Match count: 18979\n" .
                    "Mismatch count: 4837\n" .
                    "Molecule Tag (MT) count: 5123\n" .
                    "Single read category count: 18979\n";
    
    my $summ1_str = "Sample type: PE w/ Overlap\n" .
                    "Sample name: sample1\n" .
                    "ID: P1_ATCG\n" . 
                    "File format: fastq\n" .
                    "Seq count: 147203\n" .
                    "Merged count: 9534\n" .
                    "Not merged count: 137669\n" .
                    "Percent merged: 6%\n" .
                    "Match count: 6822\n" .
                    "Mismatch count: 2712\n" .
                    "Molecule Tag (MT) count: 1718\n" .
                    "Single read category count: 6822\n";
    
    open (my $OUT0, ">", $file0) or
        croak("Cannot open tmp file: $file0");
    open (my $OUT1, ">", $file1) or
        croak("Cannot open tmp file: $file1");
    
    print $OUT0 $summ0_str;
    print $OUT1 $summ1_str;
    
    close($OUT0);
    close($OUT1);
    
    is( $summ->load_sample_summaries($root_dir), 1, "load_sample_summaries()" );
    is( $summ->get_sample_count(), 2, "get_sample_count() - 2" );
    is( $summ->get_total_seqs_count(), 368323,
       "get_total-seqs_count() - 368323" );
    is( $summ->get_merged_seqs_count(), 33350,
       "get_merged_seqs_count() - 33350" );
    is( $summ->get_not_merged_seqs_count(), 334973,
       "get_not_merged_seqs_count() - 334973" );
    is( $summ->get_merged_percent(), "9%", "get_merged_percent() - 9%" );
    is( $summ->get_match_count(), 25801, "get_match_count() - 25801" );
    is( $summ->get_mismatch_count(), 7549, "get_mismatch_count() - 7549" );
    is( $summ->get_mt_count(), 6841, "get_mt_count() - 6841" );
    is( $summ->get_SRC_count(), 25801, "get_SRC_count() - 25801" );
}

# test load_sample_summaries with PE w/o Overlap output
{
    # reset summ 
    my $summ = MTToolbox::RunSummary->new();
    
    # make a directory structure similar to the output of MTToolbox
    my $root_dir = tempdir();
    mkdir "$root_dir/P0_ATTGTGAG";
    mkdir "$root_dir/P1_CGTCGGTA";
    
    my $file0 = "$root_dir/P0_ATTGTGAG/summary.txt";
    my $file1 = "$root_dir/P1_CGTCGGTA/summary.txt";
    
    my $summ0_str = "Sample type: PE w/o Overlap\n" .
                    "Sample name: sample0\n" .
                    "ID: P0_ATCG\n" .
                    "File format: fastq\n" .
                    "Seq count: 221120\n" .
                    "Match count: 18979\n" .
                    "Mismatch count: 4837\n" .
                    "Molecule Tag (MT) count: 5123\n" .
                    "Single read category count: 18979\n";
    
    my $summ1_str = "Sample type: PE w/o Overlap\n" .
                    "Sample name: sample1\n" .
                    "ID: P1_ATCG\n" .
                    "File format: fastq\n" .
                    "Seq count: 147203\n" .
                    "Match count: 6822\n" .
                    "Mismatch count: 2712\n" .
                    "Molecule Tag (MT) count: 1718\n" .
                    "Single read category count: 6822\n";
    
    open (my $OUT0, ">", $file0) or
        croak("Cannot open tmp file: $file0");
    open (my $OUT1, ">", $file1) or
        croak("Cannot open tmp file: $file1");
    
    print $OUT0 $summ0_str;
    print $OUT1 $summ1_str;
    
    close($OUT0);
    close($OUT1);
    
    is( $summ->load_sample_summaries($root_dir), 1, "load_sample_summaries()" );
    is( $summ->get_sample_count(), 2, "get_sample_count() - 2" );
    is( $summ->get_total_seqs_count(), 368323,
       "get_total-seqs_count() - 368323" );
    is( $summ->get_match_count(), 25801, "get_match_count() - 25801" );
    is( $summ->get_mismatch_count(), 7549, "get_mismatch_count() - 7549" );
    is( $summ->get_mt_count(), 6841, "get_mt_count() - 6841" );
    is( $summ->get_SRC_count(), 25801, "get_SRC_count() - 25801" );
}

# test load_sample_summaries with SE output
{
    # reset summ
    my $summ = MTToolbox::RunSummary->new();
    
    # make a directory structure similar to the output of MTToolbox
    my $root_dir = tempdir();
    mkdir "$root_dir/P0_ATTGTGAG";
    mkdir "$root_dir/P1_CGTCGGTA";
    
    my $file0 = "$root_dir/P0_ATTGTGAG/summary.txt";
    my $file1 = "$root_dir/P1_CGTCGGTA/summary.txt";
    
    my $summ0_str = "Sample type: SE\n" .
                    "Sample name: sample0\n" .
                    "ID: P0_ATCG\n" .
                    "File format: fastq\n" .
                    "Seq count: 221120\n" .
                    "Match count: 18979\n" .
                    "Mismatch count: 4837\n" .
                    "Molecule Tag (MT) count: 5123\n" .
                    "Single read category count: 18979\n";
    
    my $summ1_str = "Sample type: SE\n" .
                    "Sample name: sample1\n" .
                    "ID: P1_ATCG\n" .
                    "File format: fastq\n" .
                    "Seq count: 147203\n" .
                    "Match count: 6822\n" .
                    "Mismatch count: 2712\n" .
                    "Molecule Tag (MT) count: 1718\n" .
                    "Single read category count: 6822\n";
    
    open (my $OUT0, ">", $file0) or
        croak("Cannot open tmp file: $file0");
    open (my $OUT1, ">", $file1) or
        croak("Cannot open tmp file: $file1");
    
    print $OUT0 $summ0_str;
    print $OUT1 $summ1_str;
    
    close($OUT0);
    close($OUT1);
    
    is( $summ->load_sample_summaries($root_dir), 1, "load_sample_summaries()" );
    is( $summ->get_sample_count(), 2, "get_sample_count() - 2" );
    is( $summ->get_total_seqs_count(), 368323,
       "get_total-seqs_count() - 368323" );
    is( $summ->get_match_count(), 25801, "get_match_count() - 25801" );
    is( $summ->get_mismatch_count(), 7549, "get_mismatch_count() - 7549" );
    is( $summ->get_mt_count(), 6841, "get_mt_count() - 6841" );
    is( $summ->get_SRC_count(), 25801, "get_SRC_count() - 25801" );
}

# test load_sample_summaries with mix output types
{
    # reset summ
    my $summ = MTToolbox::RunSummary->new();
    
    # make a directory structure similar to the output of MTToolbox
    my $root_dir = tempdir();
    mkdir "$root_dir/P0_ATTGTGAG";
    mkdir "$root_dir/P1_CGTCGGTA";
    mkdir "$root_dir/P2_GAAGTGAC";
    
    my $file0 = "$root_dir/P0_ATTGTGAG/summary.txt";
    my $file1 = "$root_dir/P1_CGTCGGTA/summary.txt";
    my $file2 = "$root_dir/P2_GAAGTGAC/summary.txt";
    
    my $summ0_str = "Sample type: SE\n" .
                    "Sample name: sample0\n" .
                    "ID: P0_ATCG\n" .
                    "File format: fastq\n" .
                    "Seq count: 221120\n" .
                    "Match count: 18979\n" .
                    "Mismatch count: 4837\n" .
                    "Molecule Tag (MT) count: 5123\n" .
                    "Single read category count: 18979\n";
    
    my $summ1_str = "Sample type: PE w/o Overlap\n" .
                    "Sample name: sample1\n" .
                    "ID: P1_ATCG\n" .
                    "File format: fastq\n" .
                    "Seq count: 147203\n" .
                    "Match count: 6822\n" .
                    "Mismatch count: 2712\n" .
                    "Molecule Tag (MT) count: 1718\n" .
                    "Single read category count: 6822\n";
    
    my $summ2_str = "Sample type: PE w/ Overlap\n" .
                    "Sample name: sample2\n" .
                    "ID: P2_ATCG\n" .
                    "File format: fastq\n" .
                    "Seq count: 221120\n" .
                    "Merged count: 23816\n" .
                    "Not merged count: 197304\n" .
                    "Percent merged: 10%\n" .
                    "Match count: 18979\n" .
                    "Mismatch count: 4837\n" .
                    "Molecule Tag (MT) count: 5123\n" .
                    "Single read category count: 18979\n";
    
    open (my $OUT0, ">", $file0) or
        croak("Cannot open tmp file: $file0");
    open (my $OUT1, ">", $file1) or
        croak("Cannot open tmp file: $file1");
    open (my $OUT2, ">", $file2) or
        croak("Cannot open tmp file: $file2");
    
    print $OUT0 $summ0_str;
    print $OUT1 $summ1_str;
    print $OUT2 $summ2_str;
    
    close($OUT0);
    close($OUT1);
    close($OUT2);
    
    is( $summ->load_sample_summaries($root_dir), 1, "load_sample_summaries()" );
    is( $summ->get_sample_count(), 3, "get_sample_count() - 2" );
    is( $summ->get_total_seqs_count(), 589443,
       "get_total-seqs_count() - 589443" );
    is( $summ->get_merged_seqs_count(), 23816,
       "get_merged_seqs_count() - 33350" );
    is( $summ->get_not_merged_seqs_count(), 565627,
       "get_not_merged_seqs_count() - 565627" );
    is( $summ->get_merged_percent(), "4%", "get_merged_percent() - 4%" );
    is( $summ->get_match_count(), 44780, "get_match_count() - 44780" );
    is( $summ->get_mismatch_count(), 12386, "get_mismatch_count() - 12386" );
    is( $summ->get_mt_count(), 11964, "get_mt_count() - 11964" );
    is( $summ->get_SRC_count(), 44780, "get_SRC_count() - 44780" );
    
    # check the to_string method on this
    my $str = "Sample count: 3\n" .
              "Seq count: 589443\n" .
              "Merged count: 23816\n" .
              "Not merged count: 565627\n" .
              "Percent merged: 4%\n" .
              "Match count: 44780\n" .
              "Mismatch count: 12386\n" .
              "Molecule Tag (MT) count: 11964\n" .
              "Single read category count: 44780\n" .
              "\n\n" .
              "# Name\tID\tTotal_Seq_Count\tMerged_Count\tNot_Merged_Count\t" .
              "Percent_Merged\tMatch_Count\tMismatch_Count\tMT_Count\t" .
              "SRC_Count\n" .
              "sample0\tP0_ATCG\t221120\tNA\tNA\tNA\t18979\t4837\t5123\t18979\n" .
              "sample1\tP1_ATCG\t147203\tNA\tNA\tNA\t6822\t2712\t1718\t6822\n" .
              "sample2\tP2_ATCG\t221120\t23816\t197304\t10%\t18979\t4837\t5123\t" .
              "18979\n";
    
    is( $summ->to_string(), $str, "to_string() -- mixed sample types" );
}

# test to_string method
{
    # NOTE: I am using the $summ object created in the PE w/ Overlap test
    # So this test the to_string method coming from load_sample_summaries().
    
    my $str = "Sample count: 2\n" .
              "Seq count: 368323\n" .
              "Merged count: 33350\n" .
              "Not merged count: 334973\n" .
              "Percent merged: 9%\n" .
              "Match count: 25801\n" .
              "Mismatch count: 7549\n" .
              "Molecule Tag (MT) count: 6841\n" .
              "Single read category count: 25801\n" .
              "\n\n" .
              "# Name\tID\tTotal_Seq_Count\tMerged_Count\tNot_Merged_Count\t" .
              "Percent_Merged\tMatch_Count\tMismatch_Count\tMT_Count\t" .
              "SRC_Count\n" .
              "sample0\tP0_ATCG\t221120\t23816\t197304\t10%\t18979\t4837\t" .
              "5123\t18979\n" .
              "sample1\tP1_ATCG\t147203\t9534\t137669\t6%\t6822\t2712\t" .
              "1718\t6822\n";
    
    is( $summ->to_string(), $str, "to_string() - PE w/ overlap samples" );
}

# test print_summary_files
{
    my $root_dir = tempdir();
    $summ->print_summary_files($root_dir);
    
    cmp_ok( -s "$root_dir/all_summary.ps", '>', 0, "print_summary_files()" );
}


