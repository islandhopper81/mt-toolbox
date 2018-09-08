
use strict;
use warnings;

use MTToolbox::MoleculeTagCategory::SE;
use MTToolbox::Match::SE;
use BioUtils::FastqSeq 1.0.0;
use Test::More tests => 16;
use Test::Exception;

BEGIN { use_ok( 'MTToolbox::MoleculeTagCategory::SE'); }

### Create a MoleculeTagCategory object to test SE sequences
my $mt = "GTCGTGCAG";
my $category_se = MTToolbox::MoleculeTagCategory::SE->new({tag => $mt});

is( $category_se->get_tag(), $mt, "get_tag()" );

### Now I need to make a MTToolbox::Match::SE so I can add it to the
# MTToolbox::MoleculeTagCategory object first I need to make the
# BioUtils::FastqSeq pieces
my $fwdTagQualsStr = '@@CDFDFFH';
my $fwdTag = BioUtils::FastqSeq->new( {seq => "ATTTGAGTG", quals_str => $fwdTagQualsStr} );

my $fwdlinkerQualsStr = "HHHH";
my $fwd_linker = BioUtils::FastqSeq->new( {seq => "CAGT", quals_str => $fwdlinkerQualsStr} );

my $fwdPrimerQualsStr = "IJJIJGIIGGHIGGIIHID";
my $fwdPrimer = BioUtils::FastqSeq->new( {seq => "GTGCCAGCAGCCGCGGTAA", quals_str => $fwdPrimerQualsStr} );

my $fwdAmpliconQualStr = 'GGHFHEHIGIIJJGGIJJI>HEHDF?CCE@CD@ABBBDDDDDDDDDDDDEDD@CDDDDDDD?ACDDCDDDDDBBDAC@CD>CBDDDDCDDDDDCBABCBDBDD>BB8::AC#';
my $fwdAmplicon = BioUtils::FastqSeq->new( {seq => "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACC", quals_str => $fwdAmpliconQualStr} );

my $fwdTailQualsStr = '##';
my $fwdTail = BioUtils::FastqSeq->new( {seq => "CC", quals_str => $fwdTailQualsStr} );

# Now that all the pieces have been created I can create a
# MTToolbox::Match::SE object
my $id = 'P1_1';
my $desc = 'UNC20:2:000000000-A0UEL:1:1:13352:2073';
my $match_se = MTToolbox::Match::SE->new({id => $id,
                           desc => $desc,
                           fwd_tag => $fwdTag,
                           fwd_linker => $fwd_linker,
                           fwd_primer => $fwdPrimer,
                           fwd_amplicon => $fwdAmplicon,
                           fwd_tail => $fwdTail,
                        });

# Now that I have a MTToolbox::Match::SE object I can test the add_match method
# in MTToolbox::MoleculeTagCategory.pm
is( $category_se->add_match($match_se), 1, "add_match()" );
is( $category_se->get_seq_count(), 1, "get_seq_count()" );
my $fastaStr = '>' . $match_se->get_header() . "\n";
$fastaStr .= $match_se->get_fwd_amplicon()->get_seq() . "\n";
is( $category_se->get_unaligned_seqs_fasta_str(), $fastaStr, "get_unaligned_seqs_clustalw_str()" );

my $fastqStr = '@' . $match_se->get_header() . "\n";
$fastqStr .= $match_se->get_fwd_amplicon()->get_seq() . "\n";
$fastqStr .= "+\n";
$fastqStr .= $match_se->get_fwd_amplicon()->get_quals_str() . "\n";
is( $category_se->get_unaligned_seqs_fastq_str(), $fastqStr, "get_unaligned_seqs_fastq_str()" );





### Take two MTToolbox::Match::SE sequences through all major processes of
# MTToolbox::MoleculeTagCategory.pm (adding them, building and MSA with
# clustalw, building a consensus)

# First create a MoleculeTagCategory object to test SE sequences
my $category_se2 = MTToolbox::MoleculeTagCategory::SE->new({tag => $mt});

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

# Now that all the pieces have been created I can create a
# MTToolbox::Match::SE object
$id = 'P1_1';
$desc = 'UNC20:2:000000000-A0UEL:1:1:13352:2073';
$match_se = MTToolbox::Match::SE->new({id => $id,
                            desc => $desc,
                            fwd_tag => $fwdTag,
                            fwd_linker => $fwd_linker,
                            fwd_primer => $fwdPrimer,
                            fwd_amplicon => $fwdAmplicon,
                            fwd_tail => $fwdTail,
                        });


# Now make the second sequence which will be exactly the same but with different ID.
# I don't think I need to test this with different sequences because all that is tested
# pretty thoroughly in other modules such as FastqConseusns.pm and FastqColumn.pm.
$id = 'P1_2';
$desc = 'UNC20:2:000000000-A0UEL:1:1:13352:2074';  # I changed the last number from 3 to 4.  All other variables are the same
my $match_se2 = MTToolbox::Match::SE->new({id => $id,
                                desc => $desc,
                                fwd_tag => $fwdTag,
                                fwd_linker => $fwd_linker,
                                fwd_primer => $fwdPrimer,
                                fwd_amplicon => $fwdAmplicon,
                                fwd_tail => $fwdTail,
                            });

# Now begin add the MTToolbox::Match::SE objects
# Now that I have a MTToolbox::Match::SE object I can test the add_match method
# in MoleculeTagCategory.pm
is( $category_se2->add_match($match_se), 1, "add_match()" );
is( $category_se2->add_match($match_se2), 1, "add_match()" );
is( $category_se2->get_seq_count(), 2, "get_seq_count()" );

# Now build an MSA using clustalw.  
my $fasta_files_dir = "fasta_files/";
my $aln_files_dir = "aln_files/";
mkdir $fasta_files_dir;
mkdir $aln_files_dir;
$category_se2->run_command_line_clustalw($fasta_files_dir, $aln_files_dir);
my $fastqCon = $category_se2->build_con_from_clustalw_file($aln_files_dir);
is( $fastqCon->get_seq(), $fwdAmplicon->get_seq(),
   "run_command_line_clustalw -- fastqCon->get_seq()" );
is_deeply( $fastqCon->get_quals_aref,
          [ (split //, $fwdAmplicon->get_quals_str()) ],
          "run_command_line_clustalw -- fastqCon->getQualsArr" );


# Now build an MSA using clustalw.  
$category_se2->run_command_line_muscle($fasta_files_dir, $aln_files_dir);
$fastqCon = $category_se2->build_con_from_fasta_file($aln_files_dir);
is( $fastqCon->get_seq(), $fwdAmplicon->get_seq(),
   "run_command_line_muscle -- fastqCon->get_seq()" );
is_deeply( $fastqCon->get_quals_aref,
          [ (split //, $fwdAmplicon->get_quals_str()) ],
          "run_command_line_muscle -- fastqCon->getQualsArr" );

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
    my $fastq_con;
    lives_ok( sub{ $fastq_con = $category_se2->build_con() },
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

