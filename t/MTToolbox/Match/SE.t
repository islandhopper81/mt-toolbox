
use strict;
use warnings;

use MTToolbox::Match::SE;
use BioUtils::FastqSeq 1.0.0;
use Test::More tests => 7;


BEGIN { use_ok( 'MTToolbox::Match::SE'); }

# Create a the pieces of a MTToolbox::Match::SE object to test.
my $id = "P1_1";
my $desc = "UNC20:4:000000000-A0T9D:1:1:17685:2032";

my $fwdTagQualsStr = '@@CDFDFFH';
my $fwdTag = BioUtils::FastqSeq->new( {seq => "ATTTGAGTG", quals_str => $fwdTagQualsStr} );

my $fwdLinkerQualsStr = "HHHH";
my $fwd_linker = BioUtils::FastqSeq->new( {seq => "CAGT", quals_str => $fwdLinkerQualsStr} );

my $fwdPrimerQualsStr = "IJJIJGIIGGHIGGIIHID";
my $fwdPrimer = BioUtils::FastqSeq->new( {seq => "GTGCCAGCAGCCGCGGTAA", quals_str => $fwdPrimerQualsStr} );

my $fwdAmpliconQualStr = 'GGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC#';
my $fwdAmplicon = BioUtils::FastqSeq->new( {seq => "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACC", quals_str => $fwdAmpliconQualStr} );

my $fwdTailQualsStr = '##';
my $fwdTail = BioUtils::FastqSeq->new( {seq => "CC", quals_str => $fwdTailQualsStr} );


### Create the MTToolbox::Match::SE object
my $my_match_SE = MTToolbox::Match::SE->new({id => $id,
                                  desc => $desc,
                                  fwd_tag => $fwdTag,
                                  fwd_linker => $fwd_linker,
                                  fwd_primer => $fwdPrimer,
                                  fwd_amplicon => $fwdAmplicon,
                                  fwd_tail => $fwdTail,
                                });

### Test attribute getter subroutines
is( $my_match_SE->get_fwd_tail(), $fwdTail, "get_fwd_tail()" );

### Test attribute setter subroutines
# Make a new tail
my $fwdTailQualsStr2 = 'H#';
my $fwdTail2 = BioUtils::FastqSeq->new( {seq => "TC", quals_str => $fwdTailQualsStr2} );

# Set the new tail
is( $my_match_SE->set_fwd_tail($fwdTail2), 1, "set_fwd_tail()" );

# Check that the new tail was set correctly
is( $my_match_SE->get_fwd_tail(), $fwdTail2, "get_fwd_tail()" );

### Test other getter subroutines
is( $my_match_SE->get_fasta_str(), ">P1_1 UNC20:4:000000000-A0T9D:1:1:17685:2032\nTACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACC\n",
   "get_fasta_str()" );
is( $my_match_SE->get_fastq_str(), '@P1_1 UNC20:4:000000000-A0T9D:1:1:17685:2032' . "\n" . 'TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACC' .
   "\n+\n" .
   $fwdAmpliconQualStr .
   "\n",
   "get_fastq_str()" );
is_deeply( $my_match_SE->get_quals_aref(),
          [ (split //, $fwdAmpliconQualStr) ],
          "get_quals_aref()" );

