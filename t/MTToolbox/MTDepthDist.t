
use strict;
use warnings;

use MTToolbox::MTDepthDist;
use File::Temp qw/ tempfile tempdir /;
use Test::More tests => 20;
use Test::Exception;

sub _create_temp_file;

BEGIN { use_ok( 'MTToolbox::MTDepthDist' ); }

### Create a SampleSummary object to test
my $depth_dist;
lives_ok( sub{ $depth_dist = MTToolbox::MTDepthDist->new(); }, "new() - lives" );

# create a temp seqs_per_mt_dist.txt file
my $temp = _create_temp_file;

# test set_dist_file
{
    is( $depth_dist->get_dist_file(), undef, "get_dist_file() - undef" );
    is( $depth_dist->set_dist_file($temp), 1, "set_dist_file(temp)" );
}

# test get_dist_file
{
    is( $depth_dist->get_dist_file(), $temp, "get_dist_file() - temp" );
}

# test set_hist_href
{
    my $href = {1 => 10,
                2 => 5,
                3 => 3,
                4 => 1,
               };
    
    is( $depth_dist->set_hist_href($href), 1, "set_hist_href(href)" ) ;
}

# test get_hist_href
{
    my $href = {1 => 10,
                2 => 5,
                3 => 3,
                4 => 1,
               };
    lives_ok( sub{ $depth_dist->set_hist_href($href) },
             "set_hist_href() - lives" );
    
    my $got;
    lives_ok( sub{ $got = $depth_dist->get_hist_href() },
             "get_hist_href() - lives" );
    is_deeply($got, $href, "get_hist_href()" );
}

# test build_dist
{
    # make MoleculeTagCategory objects to test this
    ;
}

# test build_hist
{
    my $expect = {1 => 14,
                  2 => 5,
                  3 => 6,
                  4 => 5,
                  5 => 9,
                  6 => 1,
                  7 => 6,
                  8 => 1,
                  14 => 1,
                };
    
    my $got;
    lives_ok( sub{ $got = $depth_dist->build_hist() }, "build_hist() - lives" );
    is_deeply( $got, $expect, "build_hist()" );
}

# test hist_to_str
{
    my $href = {1 => 10,
                2 => 5,
                3 => 3,
                4 => 1,
               };
    
    lives_ok( sub{ $depth_dist->set_hist_href($href) },
             "set_hist_href() - lives" );
    
    my $str = "1\t10\n" .
              "2\t5\n" .
              "3\t3\n" .
              "4\t1\n";
    
    is( $depth_dist->hist_to_string(), $str, "hist_to_string()" );
}

# test print_hist
{
    my $href = {1 => 10,
                2 => 5,
                3 => 3,
                4 => 1,
               };
    
    lives_ok( sub{ $depth_dist->set_hist_href($href) },
             "set_hist_href() - lives" );
    
    my ($fh, $filename) = tempfile();
    close($fh);
    
    lives_ok( sub{ $depth_dist->print_hist($filename) },
             "print_hist() - lives");
    
    cmp_ok( -s $filename, '>', 0, "print_hist()" );
}

# test print_hist_graph
{
    lives_ok( sub{ $depth_dist->set_dist_file($temp) },
             "set_dist_file - lives" );
    lives_ok( sub{ $depth_dist->build_hist() }, "build_hist() - lives" );
    
    my $temp_dir = tempdir();
    my $file = "$temp_dir/all_seqs_per_mt_hist.ps";

    lives_ok( sub{ $depth_dist->print_hist_graph($file, "MT Depth") },
             "print_hist_graph() - lives" );
    cmp_ok( -s $file, '>', 0, "print_hist_graph()" );
}



# subroutines
sub _create_temp_file {
    
    my ($fh, $filename) = tempfile();
    
    my $str =
"P0_AAACACGTCA-AAAA      6
P0_GGGTTCGCC-TATGTC     3
P0_TGTCTTCCG-TGGAA      4
P0_GGTACTTCTA-TCAGC     4
P0_GGAACTAAGA-GCTT      1
P0_GAGTAGGAATA-TCTAT    1
P0_TGCAAGGTATA-TGCGT    7
P0_TTATTCTATTA-TCCCC    7
P0_ATTCCTCGATA-GGGC     7
P0_GTACAATGTA-TCTTA     7
P0_AATCTATTCA-TAAGAC    1
P0_GGGAAGATAA-TAACTG    1
P0_CGTAAACTTTA-ATAA     5
P0_CATTAAACTTA-TTTCG    4
P0_GCTAAGTACA-TCCCA     2
P0_AAGTAGCGGTA-ACTA     2
P0_ATGTGGTGC-TAAATC     1
P0_GAGGGACTGTA-ATGA     1
P0_TGAATCGTC-TAGCAC     5
P0_ATCTGTTTCTA-GAAG     3
P0_GTAGTGAGGA-ATGA      4
P0_CGGTTTCATTA-TTCA     5
P0_ATCAAAGAATA-TGAGA    2
P0_CTTCTTTACA-TCGT      1
P0_CTGTCTTGATA-TAGCAC   3
P0_TGAAAGTGATA-TACCTT   1
P0_GTATATAAG-GTTA       1
P0_ATACCATGTA-TAAG      1
P0_GCCGAGTCA-TACGAT     7
P0_GTATAGGAAA-TGAC      1
P0_TTCTTAAGCA-TGCAT     3
P0_TCGCGTAAAA-TAAGAT    8
P0_AGGCCAGCGA-TTGAT     5
P0_TGTCCTTGG-TCGCT      4
P0_GTCAGGTAGTA-TTTGT    5
P0_TTTTTGTCAA-GGAC      5
P0_TTGGCGCCGTA-TAACTG   3
P0_CCTAGACGTA-TACTCC    5
P0_TATGCGGTG-TGAGA      3
P0_TATTATCGA-TAGATA     1
P0_ACCCGCAAC-TGACG      2
P0_AGTCCGGTCTA-TATTC    14
P0_GTATATGCAA-TGCGA     5
P0_TAAACATAT-CCGG       1
P0_GCGAACATTTA-GATT     7
P0_GTATTCTGT-TACACT     5
P0_TTAGGGTGTTA-TGGA     1
P0_TAGGACTTGA-TGACT     2";

    print $fh $str;
    
    close($fh);
    
    return $filename;
}

