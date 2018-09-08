
use strict;
use warnings;

use MTToolbox::Match::PE::NoOverlap;
use BioUtils::FastqSeq 1.0.0;
use Test::More tests => 37;

BEGIN { use_ok( 'MTToolbox::Match::PE::NoOverlap'); }

# Create the pieces of a MTToolbox::Match::PE::NoOverlap object to test.
my $id = "P1_1";
my $desc = "UNC20:4:000000000-A0T9D:1:1:17685:2032";

my $fwd_tag_quals_str = '@@CDFDFFH';
my $fwd_tag = BioUtils::FastqSeq->new( {seq => "ATTTGAGTG", quals_str => $fwd_tag_quals_str} );

my $fwd_linker_quals_str = "HHHH";
my $fwd_linker = BioUtils::FastqSeq->new( {seq => "CAGT", quals_str => $fwd_linker_quals_str} );

my $fwd_primer_quals_str = "IJJIJGIIGGHIGGIIHID";
my $fwd_primer = BioUtils::FastqSeq->new( {seq => "GTGCCAGCAGCCGCGGTAA", quals_str => $fwd_primer_quals_str} );

my $fwd_amplicon_quals_str = 'GGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC#';
my $fwd_amplicon = BioUtils::FastqSeq->new( {seq => "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTA", quals_str => $fwd_amplicon_quals_str} );

my $fwd_tail_quals_str = '##';
my $fwd_tail = BioUtils::FastqSeq->new( {seq => "CC", quals_str => $fwd_tail_quals_str} );

my $rev_tag_quals_str = '#1:B';
my $rev_tag = BioUtils::FastqSeq->new( {seq => "AGAA", quals_str => $rev_tag_quals_str} );

my $rev_primer_quals_str = 'DDFDHHHBHIIIHIIIIIII';
my $rev_primer = BioUtils::FastqSeq->new( {seq => "GGACTACCGGGGTTTCTAAT", quals_str => $rev_primer_quals_str} );

my $rev_linker_quals_str = '????';
my $rev_linker = BioUtils::FastqSeq->new( {seq => "TGCA", quals_str => $rev_linker_quals_str} );

my $rev_amplicon_quals_str = 'IIIEHIIIGIIIIIIIIIIIIIIIIIIIHHHHGHHDEFFFFCCBBBBBCCB<ACCCECCCCBBBBCBB?@BBCECCCCBBBBBCEECBB?@BCDAA(8>B><@CA>B><B:3>:@ACCCC';
my $rev_amplicon = BioUtils::FastqSeq->new( {seq => "CCCATTCGCTCCCCTAGCTTTCGTCTCTCAGTGTCAGTGTCGGCCCAGCAGAGTGCTTTCGCCGTTGGTGTTCTTTCCGATCTCTACGCATTTCACCGCTCCACCGGAAATTCCCTCTGC", quals_str => $rev_amplicon_quals_str} );

my $rev_tail_quals_str = 'A?';
my $rev_tail = BioUtils::FastqSeq->new( {seq => "CC", quals_str => $rev_tail_quals_str} );

### Create the FastqSEMatch object
my $PE_NO = MTToolbox::Match::PE::NoOverlap->new({id => $id,
                                       desc => $desc,
                                       fwd_tag => $fwd_tag,
                                       fwd_linker => $fwd_linker,
                                       fwd_primer => $fwd_primer,
                                       fwd_amplicon => $fwd_amplicon,
                                       fwd_tail => $fwd_tail,
                                       rev_tag => $rev_tag,
                                       rev_linker => $rev_linker,
                                       rev_primer => $rev_primer,
                                       rev_amplicon => $rev_amplicon,
                                       rev_tail => $rev_tail,
                                      });


### Tests the simple getter subroutines
is( $PE_NO->get_id(), $id, "get_id()" );
is( $PE_NO->get_desc(), $desc, "get_desc()" );
is( $PE_NO->get_header(), "$id $desc", "get_header()" );
is( $PE_NO->get_fwd_tag(), $fwd_tag, "get_fwd_tag()" );
is( $PE_NO->get_fwd_linker(), $fwd_linker, "get_fwd_linker()" );
is( $PE_NO->get_fwd_primer(), $fwd_primer, "get_fwd_primer()" );
is( $PE_NO->get_fwd_amplicon(), $fwd_amplicon, "get_fwd_amplicon()" );
is( $PE_NO->get_fwd_tail(), $fwd_tail, "get_fwd_tail()" );
is( $PE_NO->get_rev_tag() , $rev_tag, "get_rev_tag()" );
is( $PE_NO->get_rev_linker(), $rev_linker, "get_rev_linker()" );
is( $PE_NO->get_rev_primer(), $rev_primer, "get_rev_primer()" );
is( $PE_NO->get_rev_amplicon(), $rev_amplicon, "get_rev_amplicon()" );
is( $PE_NO->get_rev_tail(), $rev_tail, "get_rev_tail()" );

### Test the simple setter routines (I need to make new elements).
# I made small changes to the values from the above objects
my $fwd_tag_quals_str2 = 'C@CDFDFFH';
my $fwd_tag2 = BioUtils::FastqSeq->new( {seq => "TTTTGAGTG", quals_str => $fwd_tag_quals_str2} );

my $fwd_linker_quals_str2 = "IHHH";
my $fwd_linker2 = BioUtils::FastqSeq->new( {seq => "CAGT", quals_str => $fwd_linker_quals_str2} );

my $fwd_primer_quals_str2 = "HJJIJGIIGGHIGGIIHID";
my $fwd_primer2 = BioUtils::FastqSeq->new( {seq => "GTGCCAGCAGCCGCGGTAA", quals_str => $fwd_primer_quals_str2} );

my $fwd_amplicon_quals_str2 = 'HGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC#';
my $fwd_amplicon2 = BioUtils::FastqSeq->new( {seq => "AACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTA", quals_str => $fwd_amplicon_quals_str2} );

my $fwd_tail_quals_str2 = 'H#';
my $fwd_tail2 = BioUtils::FastqSeq->new( {seq => "TC", quals_str => $fwd_tail_quals_str2} );

my $rev_tag_quals_str2 = '@1:B';
my $rev_tag2 = BioUtils::FastqSeq->new( {seq => "AGAA", quals_str => $rev_tag_quals_str2} );

my $rev_primer_quals_str2 = 'FDFDHHHBHIIIHIIIIIII';
my $rev_primer2 = BioUtils::FastqSeq->new( {seq =>"GGACTACCGGGGTTTCTAAT", quals_str => $rev_primer_quals_str2} );

my $rev_linker_quals_str2 = '????';
my $rev_linker2 = BioUtils::FastqSeq->new( {seq => "TGCA", quals_str => $rev_linker_quals_str2} );

my $rev_amplicon_quals_str2 = 'HIIEHIIIGIIIIIIIIIIIIIIIIIIIHHHHGHHDEFFFFCCBBBBBCCB<ACCCECCCCBBBBCBB?@BBCECCCCBBBBBCEECBB?@BCDAA(8>B><@CA>B><B:3>:@ACCCC';
my $rev_amplicon2 = BioUtils::FastqSeq->new( {seq => "GCCATTCGCTCCCCTAGCTTTCGTCTCTCAGTGTCAGTGTCGGCCCAGCAGAGTGCTTTCGCCGTTGGTGTTCTTTCCGATCTCTACGCATTTCACCGCTCCACCGGAAATTCCCTCTGC", quals_str => $rev_amplicon_quals_str2} );

my $rev_tail_quals_str2 = '@?';
my $rev_tail2 = BioUtils::FastqSeq->new( {seq => "CC", quals_str => $rev_tail_quals_str2} );

# Set the new values with the set subroutines
is( $PE_NO->set_fwd_tag($fwd_tag2), 1, "set_fwd_tag()" );
is( $PE_NO->set_fwd_linker($fwd_linker2), 1, "set_fwd_linker()" );
is( $PE_NO->set_fwd_primer($fwd_primer2), 1, "set_fwd_primer()" );
is( $PE_NO->set_fwd_amplicon($fwd_amplicon2), 1, "set_fwd_amplicon()" );
is( $PE_NO->set_fwd_tail($fwd_tail2), 1, "set_fwd_tail()" );
is( $PE_NO->set_rev_tag($rev_tag2), 1, "set_rev_tag()" );
is( $PE_NO->set_rev_linker($rev_linker2), 1, "set_rev_linker()" );
is( $PE_NO->set_rev_primer($rev_primer2), 1, "set_rev_primer()" );
is( $PE_NO->set_rev_amplicon($rev_amplicon2), 1, "set_rev_amplicon()" );
is( $PE_NO->set_rev_tail($rev_tail2), 1, "set_rev_tail()" );

# Test that the setters work correctly with the get subroutines
is( $PE_NO->get_fwd_tag(), $fwd_tag2, "get_fwd_tag() -- 2" );
is( $PE_NO->get_fwd_linker(), $fwd_linker2, "get_fwd_linker() -- 2" );
is( $PE_NO->get_fwd_primer(), $fwd_primer2, "get_fwd_primer() -- 2" );
is( $PE_NO->get_fwd_amplicon(), $fwd_amplicon2, "get_fwd_amplicon() -- 2" );
is( $PE_NO->get_fwd_tail(), $fwd_tail2, "get_fwd_tail() -- 2" );
is( $PE_NO->get_rev_tag() , $rev_tag2, "get_rev_tag() -- 2" );
is( $PE_NO->get_rev_linker(), $rev_linker2, "get_rev_linker() -- 2" );
is( $PE_NO->get_rev_primer(), $rev_primer2, "get_rev_primer() -- 2" );
is( $PE_NO->get_rev_amplicon(), $rev_amplicon2, "get_rev_amplicon() -- 2" );
is( $PE_NO->get_rev_tail(), $rev_tail2, "get_rev_tail() -- 2" );

### Before going on I want to reset the PE_NO object to use the first set of sequence objects
$PE_NO = MTToolbox::Match::PE::NoOverlap->new({id => $id,
                                    desc => $desc,
                                    fwd_tag => $fwd_tag,
                                    fwd_linker => $fwd_linker,
                                    fwd_primer => $fwd_primer,
                                    fwd_amplicon => $fwd_amplicon,
                                    fwd_tail => $fwd_tail,
                                    rev_tag => $rev_tag,
                                    rev_linker => $rev_linker,
                                    rev_primer => $rev_primer,
                                    rev_amplicon => $rev_amplicon,
                                    rev_tail => $rev_tail,
                                    });

### Test the more obscure get subroutines
is( $PE_NO->get_fasta_str(), ">P1_1 UNC20:4:000000000-A0T9D:1:1:17685:2032\n"
   . "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTA\n"
   . ">P1_1 UNC20:4:000000000-A0T9D:1:1:17685:2032\n"
   . "CCCATTCGCTCCCCTAGCTTTCGTCTCTCAGTGTCAGTGTCGGCCCAGCAGAGTGCTTTCGCCGTTGGTGTTCTTTCCGATCTCTACGCATTTCACCGCTCCACCGGAAATTCCCTCTGC\n"
   , "get_fasta_str()"
);
is( $PE_NO->get_fastq_str(), '@P1_1 UNC20:4:000000000-A0T9D:1:1:17685:2032'
   . "\n"
   . 'TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTA'
   . "\n+\n"
   . $fwd_amplicon_quals_str
   . "\n"
   . '@P1_1 UNC20:4:000000000-A0T9D:1:1:17685:2032'
   . "\n"
   . "CCCATTCGCTCCCCTAGCTTTCGTCTCTCAGTGTCAGTGTCGGCCCAGCAGAGTGCTTTCGCCGTTGGTGTTCTTTCCGATCTCTACGCATTTCACCGCTCCACCGGAAATTCCCTCTGC\n"
   . "+\n"
   . $rev_amplicon_quals_str
   . "\n"
   , "get_fastq_str()"
);

is_deeply( $PE_NO->get_quals_aref(), [ 
                                      (split //, $fwd_amplicon_quals_str),
                                      (split //, $rev_amplicon_quals_str),
                                      ],
                                      "get_quals_aref()"
                                      );


