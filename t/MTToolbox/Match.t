
use strict;
use warnings;

use MTToolbox::Match;
use BioUtils::FastqSeq 1.0.0;
use Test::More tests => 22;

# Remember this is kind of like an abstract class.  Some methods are not yet
# implemented, and the ones that are implemented are very general.

BEGIN { use_ok( 'MTToolbox::Match'); }

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

### Create the MTToolbox::Match::SE object
my $my_match = MTToolbox::Match->new({id => $id,
                           desc => $desc,
                           fwd_tag => $fwdTag,
                           fwd_linker => $fwd_linker,
                           fwd_primer => $fwdPrimer,
                           fwd_amplicon => $fwdAmplicon,
                        });

### Tests the simple getter subroutines
is( $my_match->get_id(), $id, "get_id()" );
is( $my_match->get_desc(), $desc, "get_desc()" );
is( $my_match->get_header(), "$id $desc", "get_header()" );
is( $my_match->get_fwd_tag(), $fwdTag, "get_fwd_tag()" );
is( $my_match->get_fwd_linker(), $fwd_linker, "get_fwd_linker()" );
is( $my_match->get_fwd_primer(), $fwdPrimer, "get_fwd_primer()" );
is( $my_match->get_fwd_amplicon(), $fwdAmplicon, "get_fwd_amplicon()" );

### Test the simple setter routines (I need to make new elements)
my $fwdTagQualsStr2 = 'C@CDFDFFH';
my $fwdTag2 = BioUtils::FastqSeq->new( {seq => "TTTTGAGTG", quals_str => $fwdTagQualsStr2} );

my $fwdLinkerQualsStr2 = "IHHH";
my $fwd_linker2 = BioUtils::FastqSeq->new( {seq => "CAGT", quals_str => $fwdLinkerQualsStr2} );

my $fwdPrimerQualsStr2 = "HJJIJGIIGGHIGGIIHID";
my $fwdPrimer2 = BioUtils::FastqSeq->new( {seq => "GTGCCAGCAGCCGCGGTAA", quals_str => $fwdPrimerQualsStr2} );

my $fwdAmpliconQualStr2 = 'HGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC#';
my $fwdAmplicon2 = BioUtils::FastqSeq->new( {seq => "AACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACC", quals_str => $fwdAmpliconQualStr2} );

my $fwdTailQualsStr2 = 'H#';
my $fwdTail2 = BioUtils::FastqSeq->new( {seq => "TC", quals_str => $fwdTailQualsStr2} );

# Set the new values with the set subroutines
is( $my_match->set_id("P2_1"), 1, "set_id()");
is( $my_match->set_desc(""), 1, "set_desc");
is( $my_match->get_header(), "P2_1", "get_header(P2_1)" );
is( $my_match->set_fwd_tag($fwdTag2), 1, "set_fwd_tag()" );
is( $my_match->set_fwd_linker($fwd_linker2), 1, "set_fwd_linker()" );
is( $my_match->set_fwd_primer($fwdPrimer2), 1, "set_fwd_primer()" );
is( $my_match->set_fwd_amplicon($fwdAmplicon2), 1, "set_fwd_amplicon()" );

# Test that the setters work correctly with the get subroutines
is( $my_match->get_id(), "P2_1", "get_id(P2_1)" );
is( $my_match->get_desc(), "", "get_desc(\"\")" );
is( $my_match->get_header(), "P2_1", "get_header(P2_1)" );
is( $my_match->get_fwd_tag(), $fwdTag2, "get_fwd_tag()" );
is( $my_match->get_fwd_linker(), $fwd_linker2, "get_fwd_linker()" );
is( $my_match->get_fwd_primer(), $fwdPrimer2, "get_fwd_primer()" );
is( $my_match->get_fwd_amplicon(), $fwdAmplicon2, "get_fwd_amplicon()" );

