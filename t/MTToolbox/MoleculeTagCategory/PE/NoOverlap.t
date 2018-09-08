
use strict;
use warnings;

use MTToolbox::MoleculeTagCategory::PE::NoOverlap;
use MTToolbox::Match::PE::NoOverlap;
use BioUtils::FastqSeq 1.0.0;
use Test::More tests => 23;
use Test::Exception;
use IPC::Cmd qw(can_run);

BEGIN { use_ok( 'MTToolbox::MoleculeTagCategory::PE::NoOverlap'); }

## Create a MoleculeTagCategory object to test.  This object will test PEs
# I will add the MTToolbox::Match::PE::NoOverlap objects after I make their
# pieces below
my $mt = "GTCGTGCAG_AGAA";
my $fmt = MTToolbox::MoleculeTagCategory::PE::NoOverlap->new({tag => $mt});

is( $fmt->get_tag(), $mt, "get_tag()" );

### Now I need to make a MTToolbox::Match::PE::NoOverlap so I can add it to the
# MTToolbox::MoleculeTagCategory::PE::NoOverlap object
# first I need to make the BioUtils::FastqSeq pieces
my $fwdTagQualsStr = '@@CDFDFFH';
my $fwdTag = BioUtils::FastqSeq->new( {seq => "ATTTGAGTG", quals_str => $fwdTagQualsStr} );

my $fwdlinkerQualsStr = "HHHH";
my $fwd_linker = BioUtils::FastqSeq->new( {seq => "CAGT", quals_str => $fwdlinkerQualsStr} );

my $fwdPrimerQualsStr = "IJJIJGIIGGHIGGIIHID";
my $fwdPrimer = BioUtils::FastqSeq->new( {seq => "GTGCCAGCAGCCGCGGTAA", quals_str => $fwdPrimerQualsStr} );

my $fwdAmpliconQualStr = 'GGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC#';
my $fwdAmplicon = BioUtils::FastqSeq->new( {seq => "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTA", quals_str => $fwdAmpliconQualStr} );

my $fwdTailQualsStr = '##';
my $fwdTail = BioUtils::FastqSeq->new( {seq => "CC", quals_str => $fwdTailQualsStr} );

my $revTagQualsStr = '#1:B';
my $revTag = BioUtils::FastqSeq->new( {seq => "AGAA", quals_str => $revTagQualsStr} );

my $rev_linker_quals_str = '????';
my $rev_linker = BioUtils::FastqSeq->new( {seq => "TGCA", quals_str => $rev_linker_quals_str} );

my $revPrimerQualsStr = 'DDFDHHHBHIIIHIIIIIII';
my $revPrimer = BioUtils::FastqSeq->new( {seq => "GGACTACCGGGGTTTCTAAT", quals_str => $revPrimerQualsStr} );

my $revAmpliconQualsStr = 'IIIEHIIIGIIIIIIIIIIIIIIIIIIIHHHHGHHDEFFFFCCBBBBBCCB<ACCCECCCCBBBBCBB?@BBCECCCCBBBBBCEECBB?@BCDAA(8>B><@CA>B><B:3>:@ACCCC';
my $revAmplicon = BioUtils::FastqSeq->new( {seq => "CCCATTCGCTCCCCTAGCTTTCGTCTCTCAGTGTCAGTGTCGGCCCAGCAGAGTGCTTTCGCCGTTGGTGTTCTTTCCGATCTCTACGCATTTCACCGCTCCACCGGAAATTCCCTCTGC", quals_str => $revAmpliconQualsStr} );

my $revTailQualsStr = 'A?';
my $revTail = BioUtils::FastqSeq->new( {seq => "CC", quals_str => $revTailQualsStr} );

# Now that all the pieces have been created I can create a
# MTToolbox::Match::PE::NoOverlap object
my $id = 'P1_1';
my $desc = 'UNC20:2:000000000-A0UEL:1:1:13352:2073';
my $mfqs = MTToolbox::Match::PE::NoOverlap->new({id => $id,
                                      desc => $desc,
                                      fwd_tag => $fwdTag,
                                      fwd_linker => $fwd_linker,
                                      fwd_primer => $fwdPrimer,
                                      fwd_amplicon => $fwdAmplicon,
                                      fwd_tail => $fwdTail,
                                      rev_tag => $revTag,
                                      rev_linker => $rev_linker,
                                      rev_primer => $revPrimer,
                                      rev_amplicon => $revAmplicon,
                                      rev_tail => $revTail,
                                      });

# Now that I have a MTToolbox::Match::PE::NoOverlap object I can test the
# add_match method in MoleculeTagCategory.pm
is( $fmt->add_match($mfqs), 1, "add_match()" );
is( $fmt->get_seq_count(), 1, "get_seq_count()" );

# Test the getter methods
my $fwdClustalwStr = '>' . $mfqs->get_header() . "\n";
$fwdClustalwStr .= $mfqs->get_fwd_amplicon()->get_seq() . "\n";
is( $fmt->get_unaligned_seqs_fasta_str('fwd'), $fwdClustalwStr,
   "get_unaligned_seqs_clustalw_str(fwd)" );

my $revClustalwStr = '>' . $mfqs->get_header() . "\n";
$revClustalwStr .= $mfqs->get_rev_amplicon()->get_seq() . "\n";
is( $fmt->get_unaligned_seqs_fasta_str('rev'), $revClustalwStr,
   "get_unaligned_seqs_clustalw_str(rev)" );

my $fwdFastqStr = '@' . $mfqs->get_header() . "\n";
$fwdFastqStr .= $mfqs->get_fwd_amplicon()->get_seq() . "\n";
$fwdFastqStr .= "+\n";
$fwdFastqStr .= $mfqs->get_fwd_amplicon()->get_quals_str() . "\n";
is($fmt->get_unaligned_seqs_fastq_str('fwd'), $fwdFastqStr,
   "get_unaligned_seqs_fastq_str(fwd)");



### Take two MTToolbox::Match::PE::NoOverlap sequences through all major processes of
# MTToolbox::MoleculeTagCategory::PE::NoOverlap.pm
# (adding them, building and MSA with clustalw, building a consensus)

# First create a MTToolbox::MoleculeTagCategory::PE::NoOverlap object to test PE sequences
# I am using the same molecule tag string as I did above
my $fmt2 = MTToolbox::MoleculeTagCategory::PE::NoOverlap->new({tag => $mt});

# Next I need to make the BioUtils::FastqSeq pieces
$fwdTagQualsStr = '@@CDFDFFH';
$fwdTag = BioUtils::FastqSeq->new( {seq => "ATTTGAGTG", quals_str => $fwdTagQualsStr} );

$fwdlinkerQualsStr = "HHHH";
$fwd_linker = BioUtils::FastqSeq->new( {seq => "CAGT", quals_str => $fwdlinkerQualsStr} );

$fwdPrimerQualsStr = "IJJIJGIIGGHIGGIIHID";
$fwdPrimer = BioUtils::FastqSeq->new( {seq => "GTGCCAGCAGCCGCGGTAA", quals_str => $fwdPrimerQualsStr} );

$fwdAmpliconQualStr = 'GGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC#';
$fwdAmplicon = BioUtils::FastqSeq->new( {seq => "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTA", quals_str => $fwdAmpliconQualStr} );

$fwdTailQualsStr = '##';
$fwdTail = BioUtils::FastqSeq->new( {seq => "CC", quals_str => $fwdTailQualsStr} );

$revTagQualsStr = '#1:B';
$revTag = BioUtils::FastqSeq->new( {seq => "AGAA", quals_str => $revTagQualsStr} );

$rev_linker_quals_str = '????';
$rev_linker = BioUtils::FastqSeq->new( {seq => "TGCA", quals_str => $rev_linker_quals_str} );

$revPrimerQualsStr = 'DDFDHHHBHIIIHIIIIIII';
$revPrimer = BioUtils::FastqSeq->new( {seq => "GGACTACCGGGGTTTCTAAT", quals_str => $revPrimerQualsStr} );

$revAmpliconQualsStr = 'IIIEHIIIGIIIIIIIIIIIIIIIIIIIHHHHGHHDEFFFFCCBBBBBCCB<ACCCECCCCBBBBCBB?@BBCECCCCBBBBBCEECBB?@BCDAA(8>B><@CA>B><B:3>:@ACCCC';
$revAmplicon = BioUtils::FastqSeq->new( {seq => "CCCATTCGCTCCCCTAGCTTTCGTCTCTCAGTGTCAGTGTCGGCCCAGCAGAGTGCTTTCGCCGTTGGTGTTCTTTCCGATCTCTACGCATTTCACCGCTCCACCGGAAATTCCCTCTGC", quals_str => $revAmpliconQualsStr} );

$revTailQualsStr = 'A?';
$revTail = BioUtils::FastqSeq->new( {seq => "CC", quals_str => $revTailQualsStr} );

# Now that all the pieces have been created I can create a
# MTToolbox::Match::PE::NoOverlap object
$id = 'P1_2';
$desc = 'UNC20:2:000000000-A0UEL:1:1:13352:2073';
$mfqs = MTToolbox::Match::PE::NoOverlap->new({id => $id,
                                   desc => $desc,
                                   fwd_tag => $fwdTag,
                                   fwd_linker => $fwd_linker,
                                   fwd_primer => $fwdPrimer,
                                   fwd_amplicon => $fwdAmplicon,
                                   fwd_tail => $fwdTail,
                                   rev_tag => $revTag,
                                   rev_linker => $rev_linker,
                                   rev_primer => $revPrimer,
                                   rev_amplicon => $revAmplicon,
                                   rev_tail => $revTail,
                                   });


# Make a second MTToolbox::Match::PE::NoOverlap with different ID
# I don't think I need to test this with different sequences because all that is tested
# pretty thoroughly in other modules such as FastqConseusns.pm and FastqColumn.pm.
$id = 'P1_3';
$desc = 'UNC20:2:000000000-A0UEL:1:1:13352:2074';  # I changed the last number from 3 to 4.  All other variables are the same
my $mfqs2 = MTToolbox::Match::PE::NoOverlap->new({id => $id,
                                       desc => $desc,
                                       fwd_tag => $fwdTag,
                                       fwd_linker => $fwd_linker,
                                       fwd_primer => $fwdPrimer,
                                       fwd_amplicon => $fwdAmplicon,
                                       fwd_tail => $fwdTail,
                                       rev_tag => $revTag,
                                       rev_linker => $rev_linker,
                                       rev_primer => $revPrimer,
                                       rev_amplicon => $revAmplicon,
                                       rev_tail => $revTail,
                                    });

# Add the MTToolbox::Match::PE::NoOverlap objects
is( $fmt2->add_match($mfqs), 1, "add_match()" );
is( $fmt2->add_match($mfqs2), 1, "add_match()" );
is( $fmt2->get_seq_count(), 2, "get_seq_count()" );

### Now build an MSA using clustalw.
### There are 3 methods that I can use to do this:
# run_command_line_clustalw OR run_bio_perl_clustalw OR run_command_line_muscle

# Build the output files directory structure
my $fasta_files_dir = "fasta_files";
my $aln_files_dir = "aln_files";
mkdir $fasta_files_dir;
mkdir $aln_files_dir;

# Option 1: run_command_line_clustalw
SKIP: {
    skip "Clustalw command not found", 2 if ( ! can_run("clustalw") );
    
    # Run clustalw
    $fmt2->run_command_line_clustalw($fasta_files_dir, $aln_files_dir);
    
    # Build both the FWD and REV consensus
    my ($fwdFastqCon, $revFastqCon) = $fmt2->build_con_from_clustalw_file("$aln_files_dir", 'fwd');
    
    # Test the FWD consensus
    is( $fwdFastqCon->get_seq(), $fwdAmplicon->get_seq(),
       "FWD consensus from build_con_from_clustalw_file -- fwdFastqCon->get_seq()" );
    is( $fwdFastqCon->get_quals_str(), $fwdAmplicon->get_quals_str(),
       "FWD consensus from build_con_from_clustalw_file -- fwdFastqCon->get_quals_str()" );
    
    # Test the REV consensus
    is( $revFastqCon->get_seq(), $revAmplicon->get_seq(),
       "REV consensus from build_con_from_clustalw_file -- revFastqCon->get_seq()" );
    is( $revFastqCon->get_quals_str(), $revAmplicon->get_quals_str(),
       "REV consensus from build_con_from_clustalw_file -- revFastqCon->get_quals_str()" );
}

# Option: run_command_line_muscle
SKIP: {
    skip "Muscle command not found", 2 if ( ! can_run("muscle") );
    
    # Run muscle
    $fmt2->run_command_line_muscle($fasta_files_dir, $aln_files_dir);
    
    # Build both the FWD and REV consensus
    my ($fwdFastqCon, $revFastqCon) = $fmt2->build_con_from_fasta_file("$aln_files_dir", 'fwd');
    
    # Test the FWD consensus
    is( $fwdFastqCon->get_seq(), $fwdAmplicon->get_seq(),
       "FWD consensus from build_con_from_fasta_file -- fwdFastqCon->get_seq()" );
    is( $fwdFastqCon->get_quals_str(), $fwdAmplicon->get_quals_str(),
       "FWD consensus from build_con_from_fasta_file -- fwdFastqCon->get_quals_str()" );
    
    # Test the REV consensus
    is( $revFastqCon->get_seq(), $revAmplicon->get_seq(),
       "REV consensus from build_con_from_fasta_file -- revFastqCon->get_seq()" );
    is( $revFastqCon->get_quals_str(), $revAmplicon->get_quals_str(),
       "REV consensus from build_con_from_fasta_file -- revFastqCon->get_quals_str()" );
}

# Delete the directories that I made to store the output files
system("rm -rf $fasta_files_dir");
system("rm -rf $aln_files_dir");



###################
# I added a new approach to building the consensus sequence.  With this new
# approach I wrote methods in the MoleculeTagCategory.pm objects called
# build_con.  This method calls the
# BioUtils::ConsensusBuilder::ConsensusBuilder::build_consensus method to build
# the consensus sequence from all the amplicon regions of the match objects
# stored in each MoleculeTagCategory object.  The following test cases are
# associated with this new method

# test build_con
{
    # I am using the two match objects that are already stored in $category_se2
    my $fwd_fastq_con;
    my $rev_fastq_con;
    lives_ok( sub{ ($fwd_fastq_con, $rev_fastq_con) = $fmt2->build_con() },
             "build_con() -- lives" );
    is ( $fwd_fastq_con->get_seq(),
        $fwdAmplicon->get_seq(),
        "build_con() - checking seq accuracy" );
    is_deeply( $fwd_fastq_con->get_quals_aref,
           [
            ( split //, $fwdAmplicon->get_quals_str() ),
           ],
           "build_con -- check qual accuracy" );
    is ( $rev_fastq_con->get_seq(),
        $revAmplicon->get_seq(),
        "build_con() - checking seq accuracy" );
    is_deeply( $rev_fastq_con->get_quals_aref,
           [
            ( split //, $revAmplicon->get_quals_str() ),
           ],
           "build_con -- check qual accuracy" );
}
