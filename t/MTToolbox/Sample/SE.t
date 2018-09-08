
use strict;
use warnings;

use MTToolbox::Sample::SE;
use BioUtils::FastqSeq 1.0.0;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;;
use Test::More tests => 54;
use Test::Exception;
use Test::Warn;

BEGIN { use_ok( 'MTToolbox::Sample::SE'); }

# Create a FastqSample to test the methods
my $id = "P1";
my $name = "sample1";
my $barcode = "ATCAC";
my $fwd_file = "";
my $fwdLinkerSeq = "CAGT";
my $fwdPrimerSeq = "GTGCCAGC[AC]GCCGCGGTAA";
my $revPrimerSeq = "GGACTAC[ACT][ACG]GGGT[AT]TCTAAT";
my $fwd_mer_len = 9;
my $fwd_max_shifts = 2;  # now this is normally 5, but for these test seqs it is 2
my $min_con_depth = 2;
my $diginorm_max = "NA";
my $outputDir = "";
my $SE = MTToolbox::Sample::SE->new({id => $id,
                          name => $name,
                          barcode => $barcode,
                          reads_file => $fwd_file,
                          fwd_linker => $fwdLinkerSeq,
                          fwd_primer => $fwdPrimerSeq,
                          rev_primer => $revPrimerSeq,
                          fwd_mer_len => $fwd_mer_len,
                          fwd_max_shifts => $fwd_max_shifts,
                          min_con_depth => $min_con_depth,
                          diginorm_max => $diginorm_max,
                          output_dir => $outputDir,
                          });

# Getters and Setters are tested in Sample.t (the parent class test file)

### The biggest thing that I need to test is categorize_by_MTs().  I am going to make a few sequences and their matching
# qual arrays to use for testing.
my $fwdSeq1 =      'GTCGTGCAGCAGTGTGCCAGCAGCCGCGGTAATACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACC';
my $fwdSeq1Quals = '@@CDFDFFHHHHHIJJIJGIIGGHIGGIIHIDGGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC###';

my ($isFwdMatch1, $fwdTag1, $fwd_linker1, $fwdPrimer1, $fwdAmplicon1, $fwdTail1) = $SE->_parse_fwd_read($fwdSeq1, $fwdSeq1Quals);

# Now I can start writing the tests for the above information (_parse_fwd_read)
# These are testing the sequences returns
is( $isFwdMatch1, 1, "_parse_fwd_read fwdSeq1 isFwdMatch" );
is( $fwdTag1->get_seq(), "GTCGTGCAG", "_parse_fwd_read fwdSeq1 fwdTag" );
is( $fwd_linker1->get_seq(), $fwdLinkerSeq, "_parse_fwd_read fwdSeq1 fwd_linker" );
is( $fwdPrimer1->get_seq(), "GTGCCAGCAGCCGCGGTAA", "_parse_fwd_read fwdSeq1 fwdPrimer" );
is( $fwdAmplicon1->get_seq(), "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTA", "_parse_fwd_read fwdSeq1 fwdAmplicon" );
is( $fwdTail1->get_seq(), "CC", "_parse_fwd_read fwdSeq1 fwdTail");

# These are testing the array returns
is ( $fwdTag1->get_quals_str(), '@@CDFDFFH', "_parse_fwd_read fwdSeq1 fwdTagQuals");
is( $fwd_linker1->get_quals_str(), 'HHHH', "_parse_fwd_read fwdSeq1 fwdLinkerQuals" );
is( $fwdPrimer1->get_quals_str(), 'IJJIJGIIGGHIGGIIHID', "_parse_fwd_read fwdSeq1 fwdPrimerQuals" );
is( $fwdAmplicon1->get_quals_str(), 'GGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC#', "_parse_fwd_read fwdSeq1 fwdAmpliconQuals" );
is( $fwdTail1->get_quals_str(), '##', "_parse_fwd_read fwdSeq1 fwdTailQuals");


# Seqeunce 2
my $fwdSeq2 =      'GTCGTGCAGACAGTGTGCCAGCAGCCGCGGTAATACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTAC';
my $fwdSeq2Quals = '@@CDFDFFHHHHHIJJIJGIIGGHIGGIIHIDGGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC###';

my ($isFwdMatch2, $fwdTag2, $fwd_linker2, $fwdPrimer2, $fwdAmplicon2, $fwdTail2) = $SE->_parse_fwd_read($fwdSeq2, $fwdSeq2Quals);

is( $isFwdMatch2, 1, "_parse_fwd_read fwdSeq2 isFwdMatch" );
is( $fwdTag2->get_seq(), "GTCGTGCAGA", "_parse_fwd_read fwdSeq2 fwdTag" );
is( $fwd_linker2->get_seq(), $fwdLinkerSeq, "_parse_fwd_read fwdSeq2 fwd_linker" );
is( $fwdPrimer2->get_seq(), "GTGCCAGCAGCCGCGGTAA", "_parse_fwd_read fwdSeq2 fwdPrimer" );
is( $fwdAmplicon2->get_seq(), "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTA", "_parse_fwd_read fwdSeq2 fwdAmplicon" );
is( $fwdTail2->get_seq(), "C", "_parse_fwd_read fwdSeq2 fwdTail");

# These are testing the array returns
is ( $fwdTag2->get_quals_str(), '@@CDFDFFHH', "_parse_fwd_read fwdSeq2 fwdTagQuals");
is( $fwd_linker2->get_quals_str(), 'HHHI', "_parse_fwd_read fwdSeq2 fwdLinkerQuals" );
is( $fwdPrimer2->get_quals_str(), 'JJIJGIIGGHIGGIIHIDG', "_parse_fwd_read fwdSeq2 fwdPrimerQuals" );
is( $fwdAmplicon2->get_quals_str(), 'GHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC##', "_parse_fwd_read fwdSeq2 fwdAmpliconQuals" );
is( $fwdTail2->get_quals_str(), '#', "_parse_fwd_read fwdSeq2 fwdTailQuals");


# Seqeunce 3
my $fwdSeq3 =      'GTCGTGCAGTACAGTGTGCCAGCAGCCGCGGTAATACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTA';
my $fwdSeq3Quals = '@@CDFDFFHHHHHIJJIJGIIGGHIGGIIHIDGGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC###';

my ($isFwdMatch3, $fwdTag3, $fwd_linker3,
    $fwdPrimer3, $fwdAmplicon3, $fwdTail3) =
    $SE->_parse_fwd_read($fwdSeq3, $fwdSeq3Quals);

is( $isFwdMatch3, 1, "_parse_fwd_read fwdSeq3 isFwdMatch" );
is( $fwdTag3->get_seq(), "GTCGTGCAGTA", "_parse_fwd_read fwdSeq3 fwdTag" );
is( $fwd_linker3->get_seq(), $fwdLinkerSeq, "_parse_fwd_read fwdSeq3 fwd_linker" );
is( $fwdPrimer3->get_seq(), "GTGCCAGCAGCCGCGGTAA", "_parse_fwd_read fwdSeq3 fwdPrimer" );
is( $fwdAmplicon3->get_seq(), "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTA", "_parse_fwd_read fwdSeq3 fwdAmplicon" );
is( $fwdTail3->get_seq(), "", "_parse_fwd_read fwdSeq3 fwdTail");

# These are testing the array returns
is ( $fwdTag3->get_quals_str(), '@@CDFDFFHHH', "_parse_fwd_read fwdSeq3 fwdTagQuals");
is( $fwd_linker3->get_quals_str(), 'HHIJ', "_parse_fwd_read fwdSeq3 fwdLinkerQuals" );
is( $fwdPrimer3->get_quals_str(), 'JIJGIIGGHIGGIIHIDGG', "_parse_fwd_read fwdSeq3 fwdPrimerQuals" );
is( $fwdAmplicon3->get_quals_str(), 'HFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC###', "_parse_fwd_read fwdSeq3 fwdAmpliconQuals" );
is( $fwdTail3->get_quals_str(), q{}, "_parse_fwd_read fwdSeq3 fwdTailQuals");


# Testing the event that we read into the reverse read primer and beyond
my $fwdSeq4 =      'GTCGTGCAGCAGTGTGCCAGCAGCCGCGGTAATACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACCGGACTACAAGGGTATCTAATATG';
my $fwdSeq4Quals = '@@CDFDFFHHHHHIJJIJGIIGGHIGGIIHIDGGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC###BABCBDBDD>BB8::AC######';

my ($isFwdMatch4, $fwdTag4, $fwd_linker4,
    $fwdPrimer4, $fwdAmplicon4, $fwdTail4) =
    $SE->_parse_fwd_read($fwdSeq4, $fwdSeq4Quals);

# Now I can start writing the tests for the above information (_parse_fwd_read)
# These are testing the sequences returns
is( $isFwdMatch4, 1, "_parse_fwd_read fwdSeq4 isFwdMatch" );
is( $fwdTag4->get_seq(), "GTCGTGCAG", "_parse_fwd_read fwdSeq4 fwdTag" );
is( $fwd_linker4->get_seq(), $fwdLinkerSeq, "_parse_fwd_read fwdSeq4 fwd_linker" );
is( $fwdPrimer4->get_seq(), "GTGCCAGCAGCCGCGGTAA", "_parse_fwd_read fwdSeq4 fwdPrimer" );
is( $fwdAmplicon4->get_seq(), "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTA", "_parse_fwd_read fwdSeq4 fwdAmplicon" );
is( $fwdTail4->get_seq(), "CC", "_parse_fwd_read fwdSeq4 fwdTail");

# These are testing the array returns
is( $fwdTag4->get_quals_str(), '@@CDFDFFH', "_parse_fwd_read fwdSeq4 fwdTagQuals");
is( $fwd_linker4->get_quals_str(), 'HHHH', "_parse_fwd_read fwdSeq4 fwdLinkerQuals" );
is( $fwdPrimer4->get_quals_str(), 'IJJIJGIIGGHIGGIIHID', "_parse_fwd_read fwdSeq4 fwdPrimerQuals" );
is( $fwdAmplicon4->get_quals_str(), 'GGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC#', "_parse_fwd_read fwdSeq4 fwdAmpliconQuals" );
is( $fwdTail4->get_quals_str(), '##', "_parse_fwd_read fwdSeq4 fwdTailQuals");


# Try creating a MTToolbox::Sample::SE with an undefined file -- the error is caugth in the constructor
dies_ok( sub { MTToolbox::Sample::SE->new({id => $id,
                       name => "sample1",
                       barcode => $barcode,
                       reads_file => undef,
                       fwd_linker => $fwdLinkerSeq,
                       fwd_primer => $fwdPrimerSeq,
                       rev_primer => $revPrimerSeq,
                       fwd_mer_len => $fwd_mer_len,
                       fwd_max_shifts => $fwd_max_shifts,
                       min_con_depth => $min_con_depth,
                       diginorm_max => $diginorm_max,
                       output_dir => $outputDir,
                       })
              },
        'Expected to die'
        );

# test categorize_by_MTs
{
    # make a MTToolbox::Sample::SE object with the correct variables
    my $in_file = dirname($0) . "/test_files/SE_seqs.fastq";
    my $dir = tempdir();
    my $SE = MTToolbox::Sample::SE->new({
                                id => $id,
                                name => "sample1",
                                barcode => $barcode,
                                reads_file => $in_file,
                                fwd_linker => "GA",
                                fwd_primer => "GTGCCAGC[AC]GCCGCGGTAA",
                                rev_primer => "GGACTAC[ACT][ACG]GGGT[AT]TCTAAT",
                                fwd_mer_len => 8,
                                fwd_max_shifts => 5,
                                min_con_depth => 6,
                                diginorm_max => 15,
                                output_dir => $dir,
                            });
    
    # run categorize_by_MTs
    lives_ok( sub{ $SE->categorize_by_MTs() }, "lives - categorize_by_MTs()" );
    
    # get the resulting MT objects
    my $MT_objs_href;
    lives_ok( sub{ $MT_objs_href = $SE->get_MT_objs_href() },
             "lives - get_MT_objs_href()" );
    
    # get the MoleculeTagCategory::SE object in the MT_objs_href
    my $MT_obj;
    lives_ok( sub{ $MT_obj = $MT_objs_href->{GGAATGACTAGGT} },
             "lives- get the MT_obj" );
    
    # make sure the MT obj only has 15 seqs in it (because of the diginorm step)
    is( $MT_obj->get_seq_count(), 15, "check diginorm" );
    
    # make sure that min_con_depth works
    $MT_obj = $MT_objs_href->{GGAAACTGCAG};
    is( $MT_obj->get_seq_count(), 5, "check min_con_depth" );
    lives_ok( sub{ $SE->build_MSAs() },
             "lives - build_MSAs()" );
    lives_ok( sub{ $SE->build_consensi() },
             "lives - build_consensi()" );
    warnings_exist{ $MT_obj->get_con_fastq_str() }
              [qr/Consensus is not defined/],
              "throws - Consensus is not defined";
}
