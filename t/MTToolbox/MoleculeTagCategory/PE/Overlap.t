
use strict;
use warnings;

use MTToolbox::MoleculeTagCategory::PE::Overlap;
use MTToolbox::Match::PE::Overlap;
use BioUtils::FastqSeq 1.0.0;
use Test::More tests => 16;
use Test::Exception;
use IPC::Cmd qw(can_run);

BEGIN { use_ok( 'MTToolbox::MoleculeTagCategory::PE::Overlap'); }

### Create a MoleculeTagCategory object to test
#   MTToolbox::MoleculeTagCategory::PE::Overlap sequences
my $mt = "GTCGTGCAG_AGAA";
my $my_mt_category = MTToolbox::MoleculeTagCategory::PE::Overlap->new({tag => $mt});

is( $my_mt_category->get_tag(), $mt, "get_tag()" );

### Now I need to make a MTToolbox::Match::PE::Overlap object so I can add it to the
# MoleculeTagCategory object.
# First I need to make the BioUtils::FastqSeq pieces
my $fwdTagQualsStr = '@@CDFDFFH';
my $fwdTag = BioUtils::FastqSeq->new( {seq => "ATTTGAGTG", quals_str => $fwdTagQualsStr} );

my $fwdLinkerQualsStr = "HHHH";
my $fwd_linker = BioUtils::FastqSeq->new( {seq => "CAGT", quals_str => $fwdLinkerQualsStr} );

my $fwdPrimerQualsStr = "IJJIJGIIGGHIGGIIHID";
my $fwdPrimer = BioUtils::FastqSeq->new( {seq => "GTGCCAGCAGCCGCGGTAA",
                                quals_str => $fwdPrimerQualsStr} );

my $fwdAmpliconQualStr = 'GGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC#';
my $fwdAmplicon = BioUtils::FastqSeq->new( {seq => "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACC",
                                  quals_str => $fwdAmpliconQualStr} );

my $revTagQualsStr = '#1:B';
my $revTag = BioUtils::FastqSeq->new( {seq => "AGAA", quals_str => $revTagQualsStr} );

my $revLinkerQualsStr = "????";
my $rev_linker = BioUtils::FastqSeq->new( {seq => "TGAC", quals_str => $revLinkerQualsStr} );

my $revPrimerQualsStr = 'DDFDHHHBHIIIHIIIIIII';
my $revPrimer = BioUtils::FastqSeq->new( {seq => "GGACTACCGGGGTTTCTAAT",
                                quals_str => $revPrimerQualsStr} );


# Now that all the pieces are created I can make a
# MTToolbox::Match::PE::Overlap object
my $id = 'P1_1';
my $desc = 'UNC20:2:000000000-A0UEL:1:1:13352:2073';
my $my_match = MTToolbox::Match::PE::Overlap->new({id => $id,
                                        desc => $desc,
                                        fwd_tag => $fwdTag,
                                        fwd_linker => $fwd_linker,
                                        fwd_primer => $fwdPrimer,
                                        fwd_amplicon => $fwdAmplicon,
                                        rev_primer => $revPrimer,
                                        rev_linker => $rev_linker,
                                        rev_tag => $revTag,
                                    });

# Now that I have a MTToolbox::Match::PE::Overlap object I can test the
# add_match method in MoleculeTagCategory.pm
is( $my_mt_category->add_match($my_match), 1, "add_match()" );
is( $my_mt_category->get_seq_count(), 1, "get_seq_count()" );
my $fastaStr = '>' . $my_match->get_header() .
               "\n" .
               $my_match->get_fwd_amplicon()->get_seq() .
               "\n"
               ;
is( $my_mt_category->get_unaligned_seqs_fasta_str(), $fastaStr,
   "get_unaligned_seqs_clustalw_str()" );

my $fastqStr = '@' . $my_match->get_header() .
               "\n" . 
               $my_match->get_fwd_amplicon()->get_seq() .
               "\n" .
               "+\n" .
               $my_match->get_fwd_amplicon()->get_quals_str() .
               "\n"
               ;
is( $my_mt_category->get_unaligned_seqs_fastq_str(), $fastqStr,
   "get_unaligned_seqs_fastq_str()" );



### Take two MTToolbox::Match::PE::Overlap sequences through all major processes
# of MoleculeTagCategory.pm
# (adding them, building and MSA with clustalw, building a consensus)

# First create a MoleculeTagCategory object to test PE::Overlap sequences
my $my_mt_category2 = MTToolbox::MoleculeTagCategory::PE::Overlap->new({tag => $mt});

# Next I need to make the BioUtils::FastqSeq pieces
$fwdTagQualsStr = '@@CDFDFFH';
$fwdTag = BioUtils::FastqSeq->new( {seq => "ATTTGAGTG", quals_str => $fwdTagQualsStr} );

$fwdLinkerQualsStr = "HHHH";
$fwd_linker = BioUtils::FastqSeq->new( {seq => "CAGT", quals_str => $fwdLinkerQualsStr} );

$fwdPrimerQualsStr = "IJJIJGIIGGHIGGIIHID";
$fwdPrimer = BioUtils::FastqSeq->new( {seq => "GTGCCAGCAGCCGCGGTAA",
                             quals_str => $fwdPrimerQualsStr} );

$fwdAmpliconQualStr = 'GGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC#';
$fwdAmplicon = BioUtils::FastqSeq->new( {seq => "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTA",
                               quals_str => $fwdAmpliconQualStr} );

$revTagQualsStr = '#1:B';
$revTag = BioUtils::FastqSeq->new( {seq => "AGAA", quals_str => $revTagQualsStr} );

$revLinkerQualsStr = "????";
$rev_linker = BioUtils::FastqSeq->new( {seq => "TGAC", quals_str => $revLinkerQualsStr} );

$revPrimerQualsStr = 'DDFDHHHBHIIIHIIIIIII';
$revPrimer = BioUtils::FastqSeq->new( {seq => "GGACTACCGGGGTTTCTAAT", quals_str => $revPrimerQualsStr} );

# Now that all the pieces have been created I can create a
# MTToolbox::Match::PE::Overlap object
$id = 'P1_2';
$desc = 'UNC20:2:000000000-A0UEL:1:1:13352:2073';
$my_match = MTToolbox::Match::PE::Overlap->new({
                                     id => $id,
                                     desc => $desc,
                                     fwd_tag => $fwdTag,
                                     fwd_linker => $fwd_linker,
                                     fwd_primer => $fwdPrimer,
                                     fwd_amplicon => $fwdAmplicon,
                                     rev_primer => $revPrimer,
                                     rev_linker => $rev_linker,
                                     rev_tag => $revTag,
                                    });


# Now make the second sequence which will be exactly the same but with different ID.
# I changed the last number from 3 to 4.  All other variables are the same
$id = 'P1_3';
$desc = 'UNC20:2:000000000-A0UEL:1:1:13352:2074';  
my $my_match2 = MTToolbox::Match::PE::Overlap->new({
                                         id => $id,
                                         desc => $desc,
                                         fwd_tag => $fwdTag,
                                         fwd_linker => $fwd_linker,
                                         fwd_primer => $fwdPrimer,
                                         fwd_amplicon => $fwdAmplicon,
                                         rev_primer => $revPrimer,
                                         rev_linker => $rev_linker,
                                         rev_tag => $revTag,
                                        });


# Now that I have a MTToolbox::Match::PE::Overlap object I can test the
# add_match method in MoleculeTagCategory.pm
is( $my_mt_category2->add_match($my_match), 1, "add_match()" );
is( $my_mt_category2->add_match($my_match2), 1, "add_match()" );
is( $my_mt_category2->get_seq_count(), 2, "get_seq_count()" );

# Now build an MSA using clustalw.  There are 3 methods that I can use to do this:
# run_command_line_clustalw OR run_bio_perl_clustalw OR run_command_line_muscle
# Option: run_command_line_clustalw
my $fasta_files_dir = "fasta_files/";
my $aln_files_dir = "aln_files/";
my $fastqCon;
mkdir $fasta_files_dir;
mkdir $aln_files_dir;
SKIP: {
    skip "Clustalw command not found", 2 if ( ! can_run("clustalw") );
    
    $my_mt_category2->run_command_line_clustalw($fasta_files_dir, $aln_files_dir);
    $fastqCon = $my_mt_category2->build_con_from_clustalw_file($aln_files_dir);
    
    is(
        $fastqCon->get_seq(),
        $fwdAmplicon->get_seq(),
        "run_command_line_clustalw -- fastqCon->get_seq()"
        );
    
    is_deeply( $fastqCon->get_quals_aref,
               [
                ( split //, $fwdAmplicon->get_quals_str() ),
               ],
               "run_command_line_clustalw -- fastqCon->get_quals_aref"
            );
}

# Option: run_command_line_muscle
SKIP: {
    skip "Muscle command not found", 2 if ( ! can_run("muscle") );
    
    $my_mt_category2->run_command_line_muscle($fasta_files_dir, $aln_files_dir);
    $fastqCon = $my_mt_category2->build_con_from_fasta_file($aln_files_dir);
    
    is(
        $fastqCon->get_seq(),
        $fwdAmplicon->get_seq(),
        "run_command_line_muscle -- fastqCon->get_seq()"
        );
    
    is_deeply( $fastqCon->get_quals_aref,
               [
                ( split //, $fwdAmplicon->get_quals_str() ),
               ],
               "run_command_line_muscle -- fastqCon->get_quals_aref"
            );
}



# Remove these temp directories and their files
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
    # I am using the two match objects that are already stored in $my_mt_category2
    my $fastq_con;
    lives_ok( sub{ $fastq_con = $my_mt_category2->build_con() },
             "build_con() -- lives" );
    is ( $fastq_con->get_seq(),
        $fwdAmplicon->get_seq(),
        "build_con() - checking seq accuracy" );
    is_deeply( $fastq_con->get_quals_aref,
           [
            ( split //, $fwdAmplicon->get_quals_str() ),
           ],
           "build_con -- check qual accuracy" );
}
