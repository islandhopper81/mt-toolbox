
use strict;
use warnings;

use MTToolbox::Sample::PE::NoOverlap;
use BioUtils::FastqSeq 1.0.0;
use Test::More tests => 103;
use Test::Exception;
use Test::Warn;
use File::Basename;
use File::Temp qw(tempfile tempdir);

my $test_count = 0;

BEGIN { use_ok( 'MTToolbox::Sample::PE::NoOverlap'); }

# Create a FastqSample to test the methods
my $id = "P1";
my $name = "sample1";
my $barcode = "ATCAC";
my $fwd_reads_file = "";
my $rev_reads_file = "";
my $fwd_linker_seq = "CAGT";
my $rev_linker_seq = "TGAC";
my $fwd_primer_seq = "GTGCCAGC[AC]GCCGCGGTAA";
my $rev_primer_seq = "GGACTAC[ACT][ACG]GGGT[AT]TCTAAT";
my $fwd_mer_len = 9;
my $rev_mer_len = 4;
my $fwd_max_shifts = 2;
my $rev_max_shifts = 2;
my $min_con_depth = 2;
my $diginorm_max = "NA";
my $PENoOverSample = MTToolbox::Sample::PE::NoOverlap->new({
                            id => $id,
                            name => $name,
                            barcode => $barcode,
                            fwd_reads_file => $fwd_reads_file,
                            rev_reads_file => $rev_reads_file,
                            fwd_linker => $fwd_linker_seq,
                            rev_linker => $rev_linker_seq,
                            fwd_primer => $fwd_primer_seq,
                            rev_primer => $rev_primer_seq,
                            fwd_mer_len => $fwd_mer_len,
                            rev_mer_len => $rev_mer_len,
                            fwd_max_shifts => $fwd_max_shifts,
                            rev_max_shifts => $rev_max_shifts,
                            min_con_depth => $min_con_depth,
                            diginorm_max => $diginorm_max,
                            output_dir => '',  # leave this empty because I temp files
                        });


# Getters and Setters are tested in Sample.t (the parent class test file)

# Test the rev_mer_len getter and setter methods
is( $PENoOverSample->get_rev_mer_len(), 4, "get_rev_mer_len()" );
is( $PENoOverSample->set_rev_mer_len(3), 1, "set_rev_mer_len()" );
is( $PENoOverSample->get_rev_mer_len(), 3, "get_rev_mer_len()" );
$PENoOverSample->set_rev_mer_len(4);  # reset the rev_mer_len attribute



### The biggest thing that I need to test is categorize_by_MTs().  I am planning on testing that in basically
# two pieces--_parse_fwd_tag and _parse_rev_tag.  So I am going to make a few sequences and their matching
# qual arrays to use for testing.
my $fwdSeq1 =      'GTCGTGCAGCAGTGTGCCAGCAGCCGCGGTAATACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACC';
my $fwdSeq1Quals = '@@CDFDFFHHHHHIJJIJGIIGGHIGGIIHIDGGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC###';

my ($is_fwd_match1, $fwd_tag1, $fwd_linker1, $fwd_primer1, $fwd_amplicon1, $fwd_tail1) = $PENoOverSample->_parse_fwd_read($fwdSeq1, $fwdSeq1Quals);

# Now I can start writing the tests for the above information (_parse_fwd_read)
# These are testing the sequences returns
is( $is_fwd_match1, 1, "_parse_fwd_read fwdSeq1 is_fwd_match" );
is( $fwd_tag1->get_seq(), "GTCGTGCAG", "_parse_fwd_read fwdSeq1 fwd_tag" );
is( $fwd_linker1->get_seq(), $fwd_linker_seq, "_parse_fwd_read fwdSeq1 fwd_linker" );
is( $fwd_primer1->get_seq(), "GTGCCAGCAGCCGCGGTAA", "_parse_fwd_read fwdSeq1 fwd_primer" );
is( $fwd_amplicon1->get_seq(), "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTA", "_parse_fwd_read fwdSeq1 fwd_amplicon" );
is( $fwd_tail1->get_seq(), "CC", "_parse_fwd_read fwdSeq1 fwd_tail");

# These are testing the array returns
is ( $fwd_tag1->get_quals_str(), '@@CDFDFFH', "_parse_fwd_read fwdSeq1 fwd_tagQuals");
is( $fwd_linker1->get_quals_str(), 'HHHH', "_parse_fwd_read fwdSeq1 fwd_linkerQuals" );
is( $fwd_primer1->get_quals_str(), 'IJJIJGIIGGHIGGIIHID', "_parse_fwd_read fwdSeq1 fwd_primerQuals" );
is( $fwd_amplicon1->get_quals_str(), 'GGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC#', "_parse_fwd_read fwdSeq1 fwd_ampliconQuals" );
is( $fwd_tail1->get_quals_str(), '##', "_parse_fwd_read fwdSeq1 fwd_tailQuals");


# Seqeunce 2
my $fwdSeq2 = "GTCGTGCAGACAGTGTGCCAGCAGCCGCGGTAATACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTAC";
my $fwdSeq2Quals = '@@CDFDFFHHHHHIJJIJGIIGGHIGGIIHIDGGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC###';

my ($is_fwd_match2, $fwd_tag2, $fwd_linker2, $fwd_primer2, $fwd_amplicon2, $fwd_tail2) = $PENoOverSample->_parse_fwd_read($fwdSeq2, $fwdSeq2Quals);

is( $is_fwd_match2, 1, "_parse_fwd_read fwdSeq2 is_fwd_match" );
is( $fwd_tag2->get_seq(), "GTCGTGCAGA", "_parse_fwd_read fwdSeq2 fwd_tag" );
is( $fwd_linker2->get_seq(), $fwd_linker_seq, "_parse_fwd_read fwdSeq2 fwd_linker" );
is( $fwd_primer2->get_seq(), "GTGCCAGCAGCCGCGGTAA", "_parse_fwd_read fwdSeq2 fwd_primer" );
is( $fwd_amplicon2->get_seq(), "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTA", "_parse_fwd_read fwdSeq2 fwd_amplicon" );
is( $fwd_tail2->get_seq(), "C", "_parse_fwd_read fwdSeq2 fwd_tail");

# These are testing the array returns
is ( $fwd_tag2->get_quals_str(), '@@CDFDFFHH', "_parse_fwd_read fwdSeq2 fwd_tagQuals");
is( $fwd_linker2->get_quals_str(), "HHHI", "_parse_fwd_read fwdSeq2 fwd_linkerQuals" );
is( $fwd_primer2->get_quals_str(), "JJIJGIIGGHIGGIIHIDG", "_parse_fwd_read fwdSeq2 fwd_primerQuals" );
is( $fwd_amplicon2->get_quals_str(), 'GHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC##', "_parse_fwd_read fwdSeq2 fwd_ampliconQuals" );
is( $fwd_tail2->get_quals_str(), "#", "_parse_fwd_read fwdSeq2 fwd_tailQuals");


# Seqeunce 3
my $fwdSeq3 = "GTCGTGCAGTACAGTGTGCCAGCAGCCGCGGTAATACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTA";
my $fwdSeq3Quals = '@@CDFDFFHHHHHIJJIJGIIGGHIGGIIHIDGGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC###';

my ($is_fwd_match3, $fwd_tag3, $fwd_linker3, $fwd_primer3, $fwd_amplicon3, $fwd_tail3) = $PENoOverSample->_parse_fwd_read($fwdSeq3, $fwdSeq3Quals);

is( $is_fwd_match3, 1, "_parse_fwd_read fwdSeq3 is_fwd_match" );
is( $fwd_tag3->get_seq(), "GTCGTGCAGTA", "_parse_fwd_read fwdSeq3 fwd_tag" );
is( $fwd_linker3->get_seq(), $fwd_linker_seq, "_parse_fwd_read fwdSeq3 fwd_linker" );
is( $fwd_primer3->get_seq(), "GTGCCAGCAGCCGCGGTAA", "_parse_fwd_read fwdSeq3 fwd_primer" );
is( $fwd_amplicon3->get_seq(), "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTA", "_parse_fwd_read fwdSeq3 fwd_amplicon" );
is( $fwd_tail3->get_seq(), "", "_parse_fwd_read fwdSeq3 fwd_tail");

# These are testing the array returns
is ( $fwd_tag3->get_quals_str(), '@@CDFDFFHHH', "_parse_fwd_read fwdSeq3 fwd_tagQuals");
is( $fwd_linker3->get_quals_str(), "HHIJ", "_parse_fwd_read fwdSeq3 fwd_linkerQuals" );
is( $fwd_primer3->get_quals_str(), "JIJGIIGGHIGGIIHIDGG", "_parse_fwd_read fwdSeq3 fwd_primerQuals" );
is( $fwd_amplicon3->get_quals_str(), 'HFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC###', "_parse_fwd_read fwdSeq3 fwd_ampliconQuals" );
is( $fwd_tail3->get_quals_str(), "", "_parse_fwd_read fwdSeq3 fwd_tailQuals");




### Test the _parse_rev_read function similarly to the above tests
my $revSeq1 = "AGAATGACGGACTACCGGGGTTTCTAATCCCATTCGCTCCCCTAGCTTTCGTCTCTCAGTGTCAGTGTCGGCCCAGCAGAGTGCTTTCGCCGTTGGTGTTCTTTCCGATCTCTACGCATTTCACCGCTCCACCGGAAATTCCCTCTGCCC";
my $revSeq1Quals = '#1:B????DDFDHHHBHIIIHIIIIIIIIIIEHIIIGIIIIIIIIIIIIIIIIIIIHHHHGHHDEFFFFCCBBBBBCCB<ACCCECCCCBBBBCBB?@BBCECCCCBBBBBCEECBB?@BCDAA(8>B><@CA>B><B:3>:@ACCCCA?';

my ($is_rev_match1, $rev_tag1, $rev_linker1, $rev_primer1, $rev_amplicon1, $rev_tail1) = $PENoOverSample->_parse_rev_read($revSeq1, $revSeq1Quals);

# Now I can start writing the tests for the above information (parserevRead)
# These are testing the sequences returns
is( $is_rev_match1, 1, "parserevRead revSeq1 is_rev_match" );
is( $rev_tag1->get_seq(), "AGAA", "parserevRead revSeq1 rev_tag" );
is( $rev_linker1->get_seq(), $rev_linker_seq, "_parse_rev_read revSeq1 rev_linker" );
is( $rev_primer1->get_seq(), "GGACTACCGGGGTTTCTAAT", "parserevRead revSeq1 rev_primer" );
is( $rev_amplicon1->get_seq(), "CCCATTCGCTCCCCTAGCTTTCGTCTCTCAGTGTCAGTGTCGGCCCAGCAGAGTGCTTTCGCCGTTGGTGTTCTTTCCGATCTCTACGCATTTCACCGCTCCACCGGAAATTCCCTCTGC", "parserevRead revSeq1 rev_amplicon" );
is( $rev_tail1->get_seq(), "CC", "parserevRead revSeq1 rev_tail");

# These are testing the array returns
is( $rev_tag1->get_quals_str(), '#1:B', "parserevRead revSeq1 rev_tagQuals");
is( $rev_linker1->get_quals_str(), '????', "parseRevRead revSeq1 rev_linker" );
is( $rev_primer1->get_quals_str(), "DDFDHHHBHIIIHIIIIIII", "parserevRead revSeq1 rev_primerQuals" );
is( $rev_amplicon1->get_quals_str(), 'IIIEHIIIGIIIIIIIIIIIIIIIIIIIHHHHGHHDEFFFFCCBBBBBCCB<ACCCECCCCBBBBCBB?@BBCECCCCBBBBBCEECBB?@BCDAA(8>B><@CA>B><B:3>:@ACCCC', "parserevRead revSeq1 rev_ampliconQuals" );
is( $rev_tail1->get_quals_str(), "A?", "parserevRead revSeq1 rev_tailQuals");


# Seqeunce 2
### Test the _parse_rev_read function similarly to the above tests
my $revSeq2 = "AGAAATGACGGACTACCGGGGTTTCTAATCCCATTCGCTCCCCTAGCTTTCGTCTCTCAGTGTCAGTGTCGGCCCAGCAGAGTGCTTTCGCCGTTGGTGTTCTTTCCGATCTCTACGCATTTCACCGCTCCACCGGAAATTCCCTCTGCC";
my $revSeq2Quals = '#1:BD????DFDHHHBHIIIHIIIIIIIIIIEHIIIGIIIIIIIIIIIIIIIIIIIHHHHGHHDEFFFFCCBBBBBCCB<ACCCECCCCBBBBCBB?@BBCECCCCBBBBBCEECBB?@BCDAA(8>B><@CA>B><B:3>:@ACCCCA?';

my ($is_rev_match2, $rev_tag2, $rev_linker2, $rev_primer2, $rev_amplicon2, $rev_tail2) = $PENoOverSample->_parse_rev_read($revSeq2, $revSeq2Quals);

# Now I can start writing the tests for the above information (parserevRead)
# These are testing the sequences returns
is( $is_rev_match2, 1, "parserevRead revSeq2 is_rev_match" );
is( $rev_tag2->get_seq(), "AGAAA", "parserevRead revSeq2 rev_tag" );
is( $rev_linker2->get_seq(), "TGAC", "parseRevRead revSeq2 rev_linker" );
is( $rev_primer2->get_seq(), "GGACTACCGGGGTTTCTAAT", "parserevRead revSeq2 rev_primer" );
is( $rev_amplicon2->get_seq(), "CCCATTCGCTCCCCTAGCTTTCGTCTCTCAGTGTCAGTGTCGGCCCAGCAGAGTGCTTTCGCCGTTGGTGTTCTTTCCGATCTCTACGCATTTCACCGCTCCACCGGAAATTCCCTCTGC", "parserevRead revSeq2 rev_amplicon" );
is( $rev_tail2->get_seq(), "C", "parserevRead revSeq2 rev_tail");

# These are testing the array returns
is( $rev_tag2->get_quals_str(), '#1:BD', "parserevRead revSeq2 rev_tagQuals");
is( $rev_linker2->get_quals_str(), '????', "parseRevRead revSeq2 rev_linker" );
is( $rev_primer2->get_quals_str(), "DFDHHHBHIIIHIIIIIIII", "parserevRead revSeq2 rev_primerQuals" );
is( $rev_amplicon2->get_quals_str(), 'IIEHIIIGIIIIIIIIIIIIIIIIIIIHHHHGHHDEFFFFCCBBBBBCCB<ACCCECCCCBBBBCBB?@BBCECCCCBBBBBCEECBB?@BCDAA(8>B><@CA>B><B:3>:@ACCCCA', "parserevRead revSeq2 rev_ampliconQuals" );
is( $rev_tail2->get_quals_str(), "?", "parserevRead revSeq2 rev_tailQuals");


# Seqeunce 3
### Test the _parse_rev_read function similarly to the above tests
my $revSeq3 = "AGAATATGACGGACTACCGGGGTTTCTAATCCCATTCGCTCCCCTAGCTTTCGTCTCTCAGTGTCAGTGTCGGCCCAGCAGAGTGCTTTCGCCGTTGGTGTTCTTTCCGATCTCTACGCATTTCACCGCTCCACCGGAAATTCCCTCTGC";
my $revSeq3Quals = '#1:BDD????FDHHHBHIIIHIIIIIIIIIIEHIIIGIIIIIIIIIIIIIIIIIIIHHHHGHHDEFFFFCCBBBBBCCB<ACCCECCCCBBBBCBB?@BBCECCCCBBBBBCEECBB?@BCDAA(8>B><@CA>B><B:3>:@ACCCCA?';

my ($is_rev_match3, $rev_tag3, $rev_linker3, $rev_primer3, $rev_amplicon3, $rev_tail3) = $PENoOverSample->_parse_rev_read($revSeq3, $revSeq3Quals);

# Now I can start writing the tests for the above information (parserevRead)
# These are testing the sequences returns
is( $is_rev_match3, 1, "parserevRead revSeq3 is_rev_match" );
is( $rev_tag3->get_seq(), "AGAATA", "parserevRead revSeq3 rev_tag" );
is( $rev_linker3->get_seq(), "TGAC", "parseRevRead revSeq3 rev_linker" );
is( $rev_primer3->get_seq(), "GGACTACCGGGGTTTCTAAT", "parserevRead revSeq3 rev_primer" );
is( $rev_amplicon3->get_seq(), "CCCATTCGCTCCCCTAGCTTTCGTCTCTCAGTGTCAGTGTCGGCCCAGCAGAGTGCTTTCGCCGTTGGTGTTCTTTCCGATCTCTACGCATTTCACCGCTCCACCGGAAATTCCCTCTGC", "parserevRead revSeq3 rev_amplicon" );
is( $rev_tail3->get_seq(), "", "parserevRead revSeq3 rev_tail");

# These are testing the array returns
is( $rev_tag3->get_quals_str(), '#1:BDD', "parserevRead revSeq3 rev_tagQuals");
is( $rev_linker3->get_quals_str(), '????', "parseRevRead revSeq3 rev_linker" );
is( $rev_primer3->get_quals_str(), "FDHHHBHIIIHIIIIIIIII", "parserevRead revSeq3 rev_primerQuals" );
is( $rev_amplicon3->get_quals_str(), 'IEHIIIGIIIIIIIIIIIIIIIIIIIHHHHGHHDEFFFFCCBBBBBCCB<ACCCECCCCBBBBCBB?@BBCECCCCBBBBBCEECBB?@BCDAA(8>B><@CA>B><B:3>:@ACCCCA?', "parserevRead revSeq3 rev_ampliconQuals" );
is( $rev_tail3->get_quals_str(), "", "parserevRead revSeq3 rev_tailQuals");

### Test the Constructor ###
# I need to add in more test here!  What happens when the file is empty? (ie. '')
# Try creating a MTToolbox::Sample::PE::NoOverlap with an undefined fwd file
my ($fh, $filename) = tempfile();
dies_ok( sub { $PENoOverSample = MTToolbox::Sample::PE::NoOverlap->new({
                                        id => $id,
                                        name => $name,
                                        barcode => $barcode,
                                        fwd_reads_file => undef,
                                        rev_reads_file => $filename,
                                        fwd_linker => $fwd_linker_seq,
                                        rev_linker => $rev_linker_seq,
                                        fwd_primer => $fwd_primer_seq,
                                        rev_primer => $rev_primer_seq,
                                        output_dir => '',  # leave this empty because I temp files
                                    })
              }, 'Expected to die' );
throws_ok( sub { $PENoOverSample = MTToolbox::Sample::PE::NoOverlap->new({
                                        id => $id,
                                        name => $name,
                                        barcode => $barcode,
                                        fwd_reads_file => undef,
                                        rev_reads_file => $filename,
                                        fwd_linker => $fwd_linker_seq,
                                        rev_linker => $rev_linker_seq,
                                        fwd_primer => $fwd_primer_seq,
                                        rev_primer => $rev_primer_seq,
                                        output_dir => '',  # leave this empty because I temp files
                                    })
              }, qr/Undefined parameter/, 'Throws Undefined parameter' );

# Try creating a PE::NoOverlap with an undefined rev file
#$PENoOverSample = MTToolbox::Sample::PE::NoOverlap->new($id, $barcode, $filename, undef, $linker_seq, $fwd_primer_seq, $rev_primer_seq);
#dies_ok( sub { $PENoOverSample->binByMTs() }, 'Expected to die' );
#throws_ok( sub { $PENoOverSample->binByMTs() }, qr/Undefined parameter/, 'Throws Undefined parameter' );

dies_ok( sub { $PENoOverSample = MTToolbox::Sample::PE::NoOverlap->new({
                                        id => $id,
                                        name => $name,
                                        barcode => $barcode,
                                        fwd_reads_file => $filename,
                                        rev_reads_file => undef,
                                        fwd_linker => $fwd_linker_seq,
                                        rev_linker => $rev_linker_seq,
                                        fwd_primer => $fwd_primer_seq,
                                        rev_primer => $rev_primer_seq,
                                        output_dir => '',  # leave this empty because I temp files
                                    })
              }, 'Expected to die' );
throws_ok( sub { $PENoOverSample = MTToolbox::Sample::PE::NoOverlap->new({
                                        id => $id,
                                        name => $name,
                                        barcode => $barcode,
                                        fwd_reads_file => $filename,
                                        rev_reads_file => undef,
                                        fwd_linker => $fwd_linker_seq,
                                        rev_linker => $rev_linker_seq,
                                        fwd_primer => $fwd_primer_seq,
                                        rev_primer => $rev_primer_seq,
                                        output_dir => '',  # leave this empty because I temp files
                                    })
              }, qr/Undefined parameter/, 'Throws Undefined parameter' );

# test categorize_by_MTs
{
    # make a Sample::SE object with the correct variables
    my $fwd_file = dirname($0) . "/test_files/fwd.fastq";
    my $rev_file = dirname($0) . "/test_files/rev.fastq";
    my $dir = tempdir();
    my $PE = MTToolbox::Sample::PE::NoOverlap->new({
                id => $id,
                name => $name,
                barcode => $barcode,
                fwd_reads_file => $fwd_file,
                rev_reads_file => $rev_file,
                fwd_linker => "GA",
                rev_linker => "AC",
                fwd_primer => "GTGCCAGC[AC]GCCGCGGTAA",
                rev_primer => "GGACTAC[ACT][ACG]GGGT[AT]TCTAAT",
                fwd_mer_len => 8,
                rev_mer_len => 5,
                fwd_max_shifts => 5,
                rev_max_shifts => 5,
                min_con_depth => 6,
                diginorm_max => 15,
                output_dir => $dir,
            });
    
    # run categorize_by_MTs
    #lives_ok( sub{ $PE->categorize_by_MTs() }, "lives - categorize_by_MTs()" );
    $PE->categorize_by_MTs();
    
    # get the resulting MT objects
    my $MT_objs_href;
    lives_ok( sub{ $MT_objs_href = $PE->get_MT_objs_href() },
             "lives - get_MT_objs_href()" );
    
    # get the MoleculeTagCategory::SE object in the MT_objs_href
    my $MT_obj;
    lives_ok( sub{ $MT_obj = $MT_objs_href->{"GGAATGACTAGGT-TACTCTG"} },
             "lives- get the MT_obj" );
    
    foreach ( keys %{$MT_objs_href} ) {
        print "$_\n";
    }
    
    # make sure the MT obj only has 15 seqs in it (because of the diginorm step)
    is( $MT_obj->get_seq_count(), 15, "check diginorm" );
    
    # make sure that min_con_depth works
    $MT_obj = $MT_objs_href->{"GGAAACTGCAG-CCCTTAG"};
    is( $MT_obj->get_seq_count(), 5, "check min_con_depth" );
    lives_ok( sub{ $PE->build_MSAs("Muscle") },
             "lives - build_MSAs(Muscle)" );
    lives_ok( sub{ $PE->build_consensi("Muscle") },
             "lives - build_consensi(muscle)" );
    throws_ok( sub{ $MT_obj->get_con_fastq_str($id, 'fwd') },
              qr/Undefined attribute value/,
              "throws- Undefined attribute value" );
}


# Testing a problematic sequence
# I found out that the problematic part was that the fwd and rev primers were
# the same
{
    $PENoOverSample = MTToolbox::Sample::PE::NoOverlap->new({
                                id => $id,
                                name => $name,
                                barcode => $barcode,
                                fwd_reads_file => $fwd_reads_file,
                                rev_reads_file => $rev_reads_file,
                                fwd_linker => "",
                                rev_linker => "",
                                fwd_primer => "GAGCGTAGGCTGT",
                                rev_primer => "GAGCGTAGGCTGT",
                                fwd_mer_len => 8,
                                rev_mer_len => 8,
                                fwd_max_shifts => 0,
                                rev_max_shifts => 0,
                                min_con_depth => $min_con_depth,
                                diginorm_max => $diginorm_max,
                                output_dir => '',  # leave this empty because I temp files
                            });
    
    $fwdSeq1 =      'GATAAAACGAGCGTAGGCTGTGACTGGAGTTCAGACGTCGCTCTTCCGATCTCACGCTGGGAGCGTAGGCTGTGTCTGAACTCCAGTCACAGCCTACGCT';
    $fwdSeq1Quals = 'CCCCCGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGEEGGGGGEGEGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG:';
    
    ($is_fwd_match1, $fwd_tag1, $fwd_linker1, $fwd_primer1, $fwd_amplicon1, $fwd_tail1) = $PENoOverSample->_parse_fwd_read($fwdSeq1, $fwdSeq1Quals);
    
    # Now I can start writing the tests for the above information (_parse_fwd_read)
    # These are testing the sequences returns
    is( $is_fwd_match1, 1, "_parse_fwd_read fwdSeq1 is_fwd_match" );
    is( $fwd_tag1->get_seq(), "GATAAAAC", "_parse_fwd_read fwdSeq1 fwd_tag" );
    is( $fwd_linker1->get_seq(), "", "_parse_fwd_read fwdSeq1 fwd_linker" );
    is( $fwd_primer1->get_seq(), "GAGCGTAGGCTGT", "_parse_fwd_read fwdSeq1 fwd_primer" );
    is( $fwd_amplicon1->get_seq(), "GACTGGAGTTCAGACGTCGCTCTTCCGATCTCACGCTGG", "_parse_fwd_read fwdSeq1 fwd_amplicon" );
    is( $fwd_tail1->get_seq(), "", "_parse_fwd_read fwdSeq1 fwd_tail");
    
    # These are testing the array returns
    is ( $fwd_tag1->get_quals_str(), 'CCCCCGGG', "_parse_fwd_read fwdSeq1 fwd_tagQuals");
    is( $fwd_linker1->get_quals_str(), '', "_parse_fwd_read fwdSeq1 fwd_linkerQuals" );
    is( $fwd_primer1->get_quals_str(), 'GGGGGGGGGGGGG', "_parse_fwd_read fwdSeq1 fwd_primerQuals" );
    is( $fwd_amplicon1->get_quals_str(), 'GGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGEE', "_parse_fwd_read fwdSeq1 fwd_ampliconQuals" );
    is( $fwd_tail1->get_quals_str(), '', "_parse_fwd_read fwdSeq1 fwd_tailQuals");
    
    # testing the corresponding rev function
    ### Test the _parse_rev_read function similarly to the above tests
    $revSeq1 =      "CCGCGAACGAGCGTAGGCTGTGACTGGAGTTCAGACGTGTGCTCTTCCGA";
    $revSeq1Quals = 'CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGDFGGGDGGGGGEC';
    
    ($is_rev_match1, $rev_tag1, $rev_linker1, $rev_primer1, $rev_amplicon1, $rev_tail1) = $PENoOverSample->_parse_rev_read($revSeq1, $revSeq1Quals);
    
    # Now I can start writing the tests for the above information (parserevRead)
    # These are testing the sequences returns
    is( $is_rev_match1, 1, "parserevRead revSeq1 is_rev_match" );
    is( $rev_tag1->get_seq(), "CCGCGAAC", "parserevRead revSeq1 rev_tag" );
    is( $rev_linker1->get_seq(), "", "_parse_rev_read revSeq1 rev_linker" );
    is( $rev_primer1->get_seq(), "GAGCGTAGGCTGT", "parserevRead revSeq1 rev_primer" );
    is( $rev_amplicon1->get_seq(), "GACTGGAGTTCAGACGTGTGCTCTTCCGA", "parserevRead revSeq1 rev_amplicon" );
    is( $rev_tail1->get_seq(), "", "parserevRead revSeq1 rev_tail");
    
    # These are testing the array returns
    is( $rev_tag1->get_quals_str(), 'CCCCCGGG', "parserevRead revSeq1 rev_tagQuals");
    is( $rev_linker1->get_quals_str(), '', "parseRevRead revSeq1 rev_linker" );
    is( $rev_primer1->get_quals_str(), "GGGGGGGGGGGGG", "parserevRead revSeq1 rev_primerQuals" );
    is( $rev_amplicon1->get_quals_str(), 'GGGGGGGGGGGGGGGGDFGGGDGGGGGEC', "parserevRead revSeq1 rev_ampliconQuals" );
    is( $rev_tail1->get_quals_str(), "", "parserevRead revSeq1 rev_tailQuals");
}
