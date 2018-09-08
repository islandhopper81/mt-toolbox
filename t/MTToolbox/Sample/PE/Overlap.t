
use strict;
use warnings;

use MTToolbox::Sample::PE::Overlap;
use BioUtils::FastqSeq 1.0.0;
use Test::More tests => 72;
use Test::Exception;
use Test::Warn;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;

BEGIN { use_ok( 'MTToolbox::Sample::PE::Overlap'); }

# Create a FastqSample to test the methods
my $id = "P1";
my $name = "sample1";
my $barcode = "ATCAC";
my $fwd_file = "";
my $fwdLinkerSeq = "CAGT";
my $revLinkerSeq = "TGAC";
my $fwd_primerSeq = "GTGCCAGC[AC]GCCGCGGTAA";
my $revPrimerSeq = "GGACTAC[ACT][ACG]GGGT[AT]TCTAAT";
my $fwd_mer_len = 9;
my $rev_mer_len = 4;
my $fwd_max_shifts = 2;
my $rev_max_shifts = 2;
my $min_con_depth = 2;
my $diginorm_max = "NA";
my $outputDir = "";
my $flash_m = 10;
my $flash_M = 70;
my $flash_x = ".20";
my $flash_p = 33;
my $flash_f = 100;
my $flash_r = 100;
my $flash_s = 20;
my $PE_overlap_obj = MTToolbox::Sample::PE::Overlap->new({
                                                id => $id,
                                                name => $name,
                                                barcode => $barcode,
                                                fwd_linker => $fwdLinkerSeq,
                                                rev_linker => $revLinkerSeq,
                                                fwd_primer => $fwd_primerSeq,
                                                rev_primer => $revPrimerSeq,
                                                fwd_mer_len => $fwd_mer_len,
                                                rev_mer_len => $rev_mer_len,
                                                fwd_max_shifts => $fwd_max_shifts,
                                                rev_max_shifts => $rev_max_shifts,
                                                min_con_depth => $min_con_depth,
                                                diginorm_max => $diginorm_max,
                                                output_dir => $outputDir,
                                                reads_file => $fwd_file,
                                                flash_m => $flash_m,
                                                flash_M => $flash_M,
                                                flash_x => $flash_x,
                                                flash_p => $flash_p,
                                                flash_r => $flash_r,
                                                flash_s => $flash_s,
                                                flash_f => $flash_f,
                                                });

# Getters and Setters are tested in Sample.t (the parent class test file)

# Test the rev_mer_len getter and setter methods
is( $PE_overlap_obj->get_rev_mer_len(), 4, "get_rev_mer_len()" );
is( $PE_overlap_obj->set_rev_mer_len(3), 1, "set_rev_mer_len()" );
is( $PE_overlap_obj->get_rev_mer_len(), 3, "get_rev_mer_len()" );
$PE_overlap_obj->set_rev_mer_len(4);  # reset the rev_mer_len attribute

# Test the _reverse_complement helper method
is( MTToolbox::Sample::PE::Overlap::_reverse_complement("ATCG"), "CGAT",
   '_reverse_complement(ATCG)' );
is( MTToolbox::Sample::PE::Overlap::_reverse_complement("[ATCG]"), "[CGAT]",
   '_reverse_complement([ATCG])' );
is( MTToolbox::Sample::PE::Overlap::_reverse_complement("A[TC]G"), "C[GA]T",
   '_reverse_complement(A[TC]G)' );

### The biggest thing that I need to test is categorize_by_MTs().  However, I
# don't directly test that.  These test cases are done on _parse_mreged_read()

# I am going to make a few sequences and their matching
my $fwdSeq1 =      'AACGGTGGCCAGTGTGCCAGCCGCCGCGGTAATACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACCAAGCTTGAGAACGGGGGGGGGAGGGGGAAATTCCGGGGGGGGGGGGAAAAGGGTAGAGAACGGAAAAAACACCAACGGCGAAAGCACTCTGCTGGGCCGACACTGACACTGAGAGACGAAAGCTAGGGGAGCGAATGGGATTAGATACCCCAGTAGTCCGTCATCTGA';
my $fwdSeq1Quals = '??AA????DDDDDDDDGGGGGGHHHHHHH9CFGHHHGFFHFHIIIHIFHHHHHHHHHCCF##H#III#H##EF##H#E#F##FG#FG#FGGGGGDD#DD#@CE#GGGC#GGG#EG8#?*??CG8#CC?#1C?#44#4C#?E#::CEC?*#:C##########################################################################AAA-##5#,5A++CE>*>ECD#BCHC?-GACAA,HFEDFHFHHFHHHHFC@=HIIIIHHHHIFFFEFFDDDDBDDBB?B?????????';

my ($is_match1, $fwd_tag1, $fwd_linker1, $fwd_primer1, $fwd_amplicon1,
    $rev_primer1, $rev_linker1, $rev_tag1) =
    $PE_overlap_obj->_parse_merged_read($fwdSeq1, $fwdSeq1Quals);

# Now I can start writing the tests for the above information (_parse_merged_read)
# These are testing the sequences returns
is( $is_match1, 1,
   "_parse_merged_read fwdSeq1 isFwdMatch" );
is( $fwd_tag1->get_seq(), "AACGGTGGC",
   "_parse_merged_read fwd_tag1" );
is( $fwd_linker1->get_seq(), $fwdLinkerSeq,
   "_parse_merged_read fwd_linker1" );
is( $fwd_primer1->get_seq(), "GTGCCAGCCGCCGCGGTAA",
   "_parse_merged_read fwd_primer1" );
is( $fwd_amplicon1->get_seq(),
   "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACCAAGCTTGAGAACGGGGGGGGGAGGGGGAAATTCCGGGGGGGGGGGGAAAAGGGTAGAGAACGGAAAAAACACCAACGGCGAAAGCACTCTGCTGGGCCGACACTGACACTGAGAGACGAAAGCTAGGGGAGCGAATGGG",
   "_parse_merged_read fwd_amplicon1"
   );
is( $rev_primer1->get_seq(), "ATTAGATACCCCAGTAGTCC",
   "_parse_merged_read rev_primer1" );
is( $rev_linker1->get_seq(), MTToolbox::Sample::PE::Overlap::_reverse_complement($revLinkerSeq),
   "_parse_merged_read rev_linker1" );
is( $rev_tag1->get_seq(), "TCTGA",
   "_parse_merged_read rev_tag1" );

# These are testing the quality values returned
is ( $fwd_tag1->get_quals_str(), '??AA????D',
    "_parse_merged_read fwd_tag_quals_str1");
is( $fwd_linker1->get_quals_str(), 'DDDD',
   "_parse_merged_read fwdSeq1 fwd_linker_quals1" );
is( $fwd_primer1->get_quals_str(), 'DDDGGGGGGHHHHHHH9CF',
   "_parse_merged_read fwdSeq1 fwd_primer_quals_str1" );
is( $fwd_amplicon1->get_quals_str(),
   'GHHHGFFHFHIIIHIFHHHHHHHHHCCF##H#III#H##EF##H#E#F##FG#FG#FGGGGGDD#DD#@CE#GGGC#GGG#EG8#?*??CG8#CC?#1C?#44#4C#?E#::CEC?*#:C##########################################################################AAA-##5#,5A++CE>*>ECD#BCHC?-GACAA,HFEDFHFHHFHHHHFC@=HIIIIHH',
   "_parse_merged_read fwdSeq1 fwd_amplicon_quals_str1"
   );
is( $rev_linker1->get_quals_str(), "????",
   "_parse_merged_read fwdSeq1 rev_linker_quals1");


# Test what happens when the fwd_linker is set to nothing
{
    $PE_overlap_obj->set_fwd_linker(undef);
    
    my $fwdSeq1 =      'AACGGTGGCGTGCCAGCCGCCGCGGTAATACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACCAAGCTTGAGAACGGGGGGGGGAGGGGGAAATTCCGGGGGGGGGGGGAAAAGGGTAGAGAACGGAAAAAACACCAACGGCGAAAGCACTCTGCTGGGCCGACACTGACACTGAGAGACGAAAGCTAGGGGAGCGAATGGGATTAGATACCCCAGTAGTCCGTCATCTGA';
    my $fwdSeq1Quals = '??AA????DDDDGGGGGGHHHHHHH9CFGHHHGFFHFHIIIHIFHHHHHHHHHCCF##H#III#H##EF##H#E#F##FG#FG#FGGGGGDD#DD#@CE#GGGC#GGG#EG8#?*??CG8#CC?#1C?#44#4C#?E#::CEC?*#:C##########################################################################AAA-##5#,5A++CE>*>ECD#BCHC?-GACAA,HFEDFHFHHFHHHHFC@=HIIIIHHHHIFFFEFFDDDDBDDBB?B?????????';
    
    my ($is_match1, $fwd_tag1, $fwd_linker1, $fwd_primer1, $fwd_amplicon1, $rev_primer1, $rev_linker1, $rev_tag1) =
        $PE_overlap_obj->_parse_merged_read($fwdSeq1, $fwdSeq1Quals);
    
    # Now I can start writing the tests for the above information (_parse_merged_read)
    # These are testing the sequences returns
    is( $is_match1, 1, "_parse_merged_read fwdSeq1 isFwdMatch" );
    is( $fwd_tag1->get_seq(), "AACGGTGGC", "_parse_merged_read fwd_tag1" );
    is( $fwd_linker1->get_seq(), "", "_parse_merged_read fwd_linker1" );
    is( $fwd_primer1->get_seq(), "GTGCCAGCCGCCGCGGTAA", "_parse_merged_read fwd_primer1" );
    is( $fwd_amplicon1->get_seq(),
       "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACCAAGCTTGAGAACGGGGGGGGGAGGGGGAAATTCCGGGGGGGGGGGGAAAAGGGTAGAGAACGGAAAAAACACCAACGGCGAAAGCACTCTGCTGGGCCGACACTGACACTGAGAGACGAAAGCTAGGGGAGCGAATGGG",
       "_parse_merged_read fwd_amplicon1"
       );
    is( $rev_primer1->get_seq(), "ATTAGATACCCCAGTAGTCC", "_parse_merged_read rev_primer1" );
    is( $rev_linker1->get_seq(), MTToolbox::Sample::PE::Overlap::_reverse_complement($revLinkerSeq), "_parse_merged_read rev_linker1" );
    is( $rev_tag1->get_seq(), "TCTGA", "_parse_merged_read rev_tag1" );
    
    # These are testing the quality values returned
    is ( $fwd_tag1->get_quals_str(), '??AA????D', "_parse_merged_read fwd_tag_quals_str1");
    is( $fwd_linker1->get_quals_str(), '', "_parse_merged_read fwdSeq1 fwd_linker_quals1" );
    is( $fwd_primer1->get_quals_str(), 'DDDGGGGGGHHHHHHH9CF', "_parse_merged_read fwdSeq1 fwd_primer_quals_str1" );
    is( $fwd_amplicon1->get_quals_str(),
       'GHHHGFFHFHIIIHIFHHHHHHHHHCCF##H#III#H##EF##H#E#F##FG#FG#FGGGGGDD#DD#@CE#GGGC#GGG#EG8#?*??CG8#CC?#1C?#44#4C#?E#::CEC?*#:C##########################################################################AAA-##5#,5A++CE>*>ECD#BCHC?-GACAA,HFEDFHFHHFHHHHFC@=HIIIIHH',
       "_parse_merged_read fwdSeq1 fwd_amplicon_quals_str1"
    );
    is( $rev_linker1->get_quals_str(), "????", "_parse_merged_read fwdSeq1 rev_linker_quals1");
    
    # reset the fwd_linker from the above sequence tests
    $PE_overlap_obj->set_fwd_linker($fwdLinkerSeq);
}

# Test what happens when the rev_linker is set to nothing
{
    $PE_overlap_obj->set_rev_linker(undef);
    
    my $fwdSeq1 =      'AACGGTGGCCAGTGTGCCAGCCGCCGCGGTAATACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACCAAGCTTGAGAACGGGGGGGGGAGGGGGAAATTCCGGGGGGGGGGGGAAAAGGGTAGAGAACGGAAAAAACACCAACGGCGAAAGCACTCTGCTGGGCCGACACTGACACTGAGAGACGAAAGCTAGGGGAGCGAATGGGATTAGATACCCCAGTAGTCCTCTGA';
    my $fwdSeq1Quals = '??AA????DDDDDDDDGGGGGGHHHHHHH9CFGHHHGFFHFHIIIHIFHHHHHHHHHCCF##H#III#H##EF##H#E#F##FG#FG#FGGGGGDD#DD#@CE#GGGC#GGG#EG8#?*??CG8#CC?#1C?#44#4C#?E#::CEC?*#:C##########################################################################AAA-##5#,5A++CE>*>ECD#BCHC?-GACAA,HFEDFHFHHFHHHHFC@=HIIIIHHHHIFFFEFFDDDDBDDBB?B?????';
    
    my ($is_match1, $fwd_tag1, $fwd_linker1, $fwd_primer1, $fwd_amplicon1, $rev_primer1, $rev_linker1, $rev_tag1) =
        $PE_overlap_obj->_parse_merged_read($fwdSeq1, $fwdSeq1Quals);
    
    # Now I can start writing the tests for the above information (_parse_merged_read)
    # These are testing the sequences returns
    is( $is_match1, 1, "_parse_merged_read fwdSeq1 isFwdMatch" );
    is( $fwd_tag1->get_seq(), "AACGGTGGC", "_parse_merged_read fwd_tag1" );
    is( $fwd_linker1->get_seq(), "CAGT", "_parse_merged_read fwd_linker1" );
    is( $fwd_primer1->get_seq(), "GTGCCAGCCGCCGCGGTAA", "_parse_merged_read fwd_primer1" );
    is( $fwd_amplicon1->get_seq(),
       "TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACCAAGCTTGAGAACGGGGGGGGGAGGGGGAAATTCCGGGGGGGGGGGGAAAAGGGTAGAGAACGGAAAAAACACCAACGGCGAAAGCACTCTGCTGGGCCGACACTGACACTGAGAGACGAAAGCTAGGGGAGCGAATGGG",
       "_parse_merged_read fwd_amplicon1"
       );
    is( $rev_primer1->get_seq(), "ATTAGATACCCCAGTAGTCC", "_parse_merged_read rev_primer1" );
    is( $rev_linker1->get_seq(), "", "_parse_merged_read rev_linker1" );
    is( $rev_tag1->get_seq(), "TCTGA", "_parse_merged_read rev_tag1" );
    
    # These are testing the quality values returned
    is ( $fwd_tag1->get_quals_str(), '??AA????D', "_parse_merged_read fwd_tag_quals_str1");
    is( $fwd_linker1->get_quals_str(), 'DDDD', "_parse_merged_read fwdSeq1 fwd_linker_quals1" );
    is( $fwd_primer1->get_quals_str(), 'DDDGGGGGGHHHHHHH9CF', "_parse_merged_read fwdSeq1 fwd_primer_quals_str1" );
    is( $fwd_amplicon1->get_quals_str(),
       'GHHHGFFHFHIIIHIFHHHHHHHHHCCF##H#III#H##EF##H#E#F##FG#FG#FGGGGGDD#DD#@CE#GGGC#GGG#EG8#?*??CG8#CC?#1C?#44#4C#?E#::CEC?*#:C##########################################################################AAA-##5#,5A++CE>*>ECD#BCHC?-GACAA,HFEDFHFHHFHHHHFC@=HIIIIHH',
       "_parse_merged_read fwdSeq1 fwd_amplicon_quals_str1"
    );
    is( $rev_linker1->get_quals_str(), '', "_parse_merged_read fwdSeq1 rev_linker_quals1");
    
    # reset the rev_linker from the above sequence tests
    $PE_overlap_obj->set_rev_linker($revLinkerSeq);
}


## Seqeunce 2
{
    $PE_overlap_obj->set_fwd_linker("GA");
    $PE_overlap_obj->set_rev_linker("AC");
    
    my $fwdSeq =       'GGTGGTTAGAGAGTGCCAGCCGCCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATTACTGGGCGTAAAGGGAGCGTAGGCGGACATTTAAGTCAGGGGTGAAATCCCAGAGCTCAACTCTGGAACTGCCTTTGATACTGGGTGTCTTGAGTGTGAGAGAGGTATGTGGAACTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGACATACTGGCTCATTACGGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCGAGTAGTCCGTGCGT';
    my $fwdSeqQuals = q{????????DDD<BBBDCFFCCFCEF>CCCHHIHEGF8AFE/AFDFHIHHCE*55CEHEEEHGHIIIH<#DFGHDHH,??C+#DEBCE8EDEFFFEEDEDEEFBEFEFFEFFFFFEBCBEEEEE=5#?EDD??'8?:ACC*?ACEE<AEEEEEEFF;?EFFFC?AD<8:A?:EEEEEFC#EEE@@@@@EEEEEEEEEEEDCEEDEFFCE=?EEDEEFFEEEEDEE=EDCECACCCCAAEEA??C##:AEE#@EC8FFEEA>7,7.-.8+++-+@A9/.@GDECCEEAAEA>AA6E@@-@@@<@><<????5};
    
    my ($is_match, $fwd_tag, $fwd_linker, $fwd_primer, $fwd_amplicon, $rev_primer, $rev_linker, $rev_tag) =
        $PE_overlap_obj->_parse_merged_read($fwdSeq, $fwdSeqQuals);
    
    is( $is_match, 1, "_parse_merged_read - is_match" );
    is( $fwd_tag->get_seq(), "GGTGGTTAGA", "_parse_merged_read - fwd_tag" );
    is( $fwd_linker->get_seq(), "GA", "_parse_merged_read - fwd_linker" );
    is( $fwd_primer->get_seq(), "GTGCCAGCCGCCGCGGTAA", "_parse_merged_read - fwd_primer" );
    is( $fwd_amplicon->get_seq(), "TACGAAGGGGGCTAGCGTTGCTCGGAATTACTGGGCGTAAAGGGAGCGTAGGCGGACATTTAAGTCAGGGGTGAAATCCCAGAGCTCAACTCTGGAACTGCCTTTGATACTGGGTGTCTTGAGTGTGAGAGAGGTATGTGGAACTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGACATACTGGCTCATTACGGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAGG",
       "_parse_merged_read - fwd_amplicon" );
    is( $rev_primer->get_seq(), "ATTAGATACCCGAGTAGTCC", "_parse_merged_read rev_primer" );
    is( $rev_linker->get_seq(), "GT", "_parse_merged_read - rev_linker" );
    is( $rev_tag->get_seq(), "GCGT", "_parse_merged_read rev_tag" );
    
    # These are testing the array returns
    is( $fwd_tag->get_quals_str(), '????????DD', "_parse_merged_read fwdSeq fwd_tag_quals_str");
    is( $fwd_linker->get_quals_str(), 'D<', "_parse_merged_read - fwd_linker quals" );
    is( $fwd_primer->get_quals_str(), 'BBBDCFFCCFCEF>CCCHH', "_parse_merged_read - fwd primer quals" );
    is( $fwd_amplicon->get_quals_str(), q{IHEGF8AFE/AFDFHIHHCE*55CEHEEEHGHIIIH<#DFGHDHH,??C+#DEBCE8EDEFFFEEDEDEEFBEFEFFEFFFFFEBCBEEEEE=5#?EDD??'8?:ACC*?ACEE<AEEEEEEFF;?EFFFC?AD<8:A?:EEEEEFC#EEE@@@@@EEEEEEEEEEEDCEEDEFFCE=?EEDEEFFEEEEDEE=EDCECACCCCAAEEA??C##:AEE#@EC8FFEEA>7,7.-.8+++-+@A9/.@GDECCE},
       "_parse_merged_read - fwd amplicon quals" );
    is( $rev_primer->get_quals_str(), 'EAAEA>AA6E@@-@@@<@><', "_parse_merged_read - rev primer quals");
    is( $rev_linker->get_quals_str(), "<?", "_parse_merged_read - rev_linker quals" );
    is( $rev_tag->get_quals_str(), "???5", "_parse_merged_read - reg tag quals" );
}

# I found a bug that caused this sequence to pass with the incorrect rev tag.
# I am adding test cases for this situation.
{
    # first test that the sequence doesn't pass
    my $fwdSeq =       'TGGACCCGCTAGAGTGCCAGCCGCCGCGGTAATACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGGTGGCTTTTTAAGTCCGCCGTCAAATCCCAGGGCTCAACCCTGGACAGGCGGTGGAAACTACCAAGCTTGAGTACGGTAGGGGCAGAGGGAATTTCCGGTGGAGCGGTGAAATGCGTAGAGATCGGAAAGAACACCAACGGCTAAAGCACTCTGCTGGGCCGACACTGNCACTGAGAGACGAAAGCTAGGGGAGCGAATGTGATTAGATACCCCTGTAGTCCGTTAAGTT';
    my $fwdSeqQuals = q{????ABABDDDDDDDDFGGGGGCHHHHHHHHIGHHIHHHIIHGHIIIIHHHHIIIHHHHIIH#HH#HHHHH#FHH#HDEGGGG#G#GGGGGGGGGGGGGGGGEGEGGGGGGG?:8?C??C:0000:*:C:00?>>24?C:C*:CCCEE#<=<E=#</'///6E?66;6<EEE;66#-';EEE@888::#@E##EEE@EEEEEE@DEEEEEEEEEEEEEEEDDC6+#+DD#BD<AECEEE#EC7EFFEEA55#ECDFDCA+EEEFFFFFDEBEDCDC@EFFGE88EA.CA88+ECEE@@@-=55+=>>=====9};
    
    my ($is_match, $fwd_tag, $fwd_linker, $fwd_primer, $fwd_amplicon, $rev_primer, $rev_linker, $rev_tag) =
        $PE_overlap_obj->_parse_merged_read($fwdSeq, $fwdSeqQuals);
    
    is( $is_match, 1, "_parse_merged_read - is_match" );
    is( $rev_tag->get_seq(), "TAAGTT", "_parse_merged_read - rev_tag seq" );
}

# test _get_flash_params_str -- could use futher testing
{
    my $str = "-m 10 -M 70 -x .20 -p 33 -r 100 -f 100 -s 20 ";
    is ($PE_overlap_obj->_get_flash_params_str(), $str,
        "_get_flash_params_str");
}

# test categorize_by_MTs
{
    # make a MTToolbox::Sample::PE::Overlap object with the correct variables
    my $in_file = dirname($0) . "/test_files/PE_overlap_seqs.fastq";
    my $dir = tempdir();
    my $PE = MTToolbox::Sample::PE::Overlap->new({
                                        id => $id,
                                        name => $name,
                                        barcode => $barcode,
                                        reads_file => $in_file,
                                        fwd_linker => "GA",
                                        rev_linker => "AC",
                                        fwd_primer => "GTGCCAGC[AC]GCCGCGGTAA",
                                        rev_primer => "GGACTAC[ACT][ACG]GGGT[AT]TCTAAT",
                                        fwd_mer_len => 8,
                                        fwd_max_shifts => 5,
                                        rev_mer_len => 5,
                                        rev_max_shifts => 5,
                                        min_con_depth => 6,
                                        diginorm_max => 15,
                                        output_dir => $dir,
                                        });
    
    # run categorize_by_MTs
    lives_ok( sub{ $PE->categorize_by_MTs() }, "lives - categorize_by_MTs()" );
    
    # get the resulting MT objects
    my $MT_objs_href;
    lives_ok( sub{ $MT_objs_href = $PE->get_MT_objs_href() },
             "lives - get_MT_objs_href()" );
    
    # get the MoleculeTagCategory::SE object in the MT_objs_href
    my $MT_obj;
    lives_ok( sub{ $MT_obj = $MT_objs_href->{"GGAATGACTAGGT-CAGAGTA"} },
             "lives- get the MT_obj" );
    
    # make sure the MT obj only has 15 seqs in it (because of the diginorm step)
    is( $MT_obj->get_seq_count(), 15, "check diginorm" );
    
    # make sure that min_con_depth works
    $MT_obj = $MT_objs_href->{"GGAAACTGCAG-CTAAGGG"};
    is( $MT_obj->get_seq_count(), 5, "check min_con_depth" );
    lives_ok( sub{ $PE->build_MSAs("Muscle") },
             "lives - build_MSAs(Muscle)" );
    lives_ok( sub{ $PE->build_consensi("Muscle") },
             "lives - build_consensi(Muscle)" );
    warnings_exist{ $MT_obj->get_con_fastq_str() }
              [qr/Consensus is not defined/],
              "throws - Consensus is not defined";
}

