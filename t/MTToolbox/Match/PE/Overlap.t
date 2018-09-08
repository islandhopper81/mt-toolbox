
use strict;
use warnings;

use MTToolbox::Match::PE::Overlap;
use BioUtils::FastqSeq 1.0.0;
use Test::More tests => 4;


BEGIN { use_ok( 'MTToolbox::Match::PE::Overlap'); }

# Create a the pieces of a Match::SE object to test.
my $id = "P1_1";
my $desc = "UNC20:4:000000000-A0T9D:1:1:17685:2032";

my $fwd_tag_quals_str = '@@CDFDFFH';
my $fwd_tag = BioUtils::FastqSeq->new( {seq => "ATTTGAGTG", quals_str => $fwd_tag_quals_str} );

my $fwd_linker_quals_str = "HHHH";
my $fwd_linker = BioUtils::FastqSeq->new( {seq => "CAGT", quals_str => $fwd_linker_quals_str} );

my $fwd_primer_quals_str = "IJJIJGIIGGHIGGIIHID";
my $fwd_primer = BioUtils::FastqSeq->new( {seq => "GTGCCAGCAGCCGCGGTAA", quals_str => $fwd_primer_quals_str} );

my $fwd_amplicon_quals_str = 'GGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC#';
my $fwd_amplicon = BioUtils::FastqSeq->new( {seq => "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACC", quals_str => $fwd_amplicon_quals_str} );

my $rev_tag_quals_str = '#1:B';
my $rev_tag = BioUtils::FastqSeq->new( {seq => "AGAA", quals_str => $rev_tag_quals_str} );

my $rev_linker_quals_str = "????";
my $rev_linker = BioUtils::FastqSeq->new( {seq => "TGAC", quals_str => $fwd_linker_quals_str} );

my $rev_primer_quals_str = 'DDFDHHHBHIIIHIIIIIII';
my $rev_primer = BioUtils::FastqSeq->new( {seq => "GGACTACCGGGGTTTCTAAT", quals_str => $rev_primer_quals_str} );


### Create the Match::SE object
my $my_match = MTToolbox::Match::PE::Overlap->new({
                                           id => $id,
                                           desc => $desc,
                                           fwd_tag => $fwd_tag,
                                           fwd_linker => $fwd_linker,
                                           fwd_primer => $fwd_primer,
                                           fwd_amplicon => $fwd_amplicon,
                                           rev_tag => $rev_tag,
                                           rev_linker => $rev_linker,
                                           rev_primer => $rev_primer,
                                        });

### Test attribute getter subroutines


### Test attribute setter subroutines


### Test other getter subroutines
is( $my_match->get_fasta_str(),
    ">P1_1 UNC20:4:000000000-A0T9D:1:1:17685:2032\n" .
    $fwd_amplicon->get_seq() .
    "\n",
    "get_fasta_str()"
    );
is( $my_match->get_fastq_str(),
   '@P1_1 UNC20:4:000000000-A0T9D:1:1:17685:2032' .
   "\n" .
   $fwd_amplicon->get_seq() .
   "\n+\n" .
   $fwd_amplicon_quals_str .
   "\n",
   "get_fastq_str()"
   );
is_deeply( $my_match->get_quals_aref(),
          [ 
            (split //, $fwd_amplicon_quals_str),
          ],
          "get_quals_aref()"
          );


