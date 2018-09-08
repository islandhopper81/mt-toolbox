
use strict;
use warnings;

use MTToolbox::Sample;
use BioUtils::FastqSeq 1.0.0;
use Test::More tests => 35;

BEGIN { use_ok( 'MTToolbox::Sample'); }

# Create a FastqMTToolbox::Sample to test the methods
my $id = "P1";
my $name = "sample1";
my $barcode = "ATCAC";
my $fwdlinkerSeq = "CAGT";
my $fwdPrimerSeq = "GTGCCAGC[AC]GCCGCGGTAA";
my $revPrimerSeq = "GGACTAC[ACT][ACG]GGGT[AT]TCTAAT";
my $fwd_mer_len = 9;
my $fwd_max_shifts = 6;
my $min_con_depth = 2;
my $diginorm_max = "NA";
my $DIGI_MAX = 1000000;
my $outputDir = "";
my $sample = MTToolbox::Sample->new({id => $id,
                          name => $name,
                          barcode => $barcode,
                          fwd_linker => $fwdlinkerSeq,
                          fwd_primer => $fwdPrimerSeq,
                          rev_primer => $revPrimerSeq,
                          fwd_mer_len => $fwd_mer_len,
                          fwd_max_shifts => $fwd_max_shifts,
                          min_con_depth => $min_con_depth,
                          diginorm_max => $diginorm_max,
                          output_dir => $outputDir,
                          });

# Test getters
is( $sample->get_id(), $id, "get_id()" );
is( $sample->get_name(), $name, "get_name()" );
is( $sample->get_barcode(), $barcode, "get_barcode()" );
is( $sample->get_fwd_linker(), $fwdlinkerSeq, "get_fwd_linker()" );
is( $sample->get_fwd_primer(), $fwdPrimerSeq, "get_fwd_primer()" );
is( $sample->get_rev_primer(), $revPrimerSeq, "get_rev_primer()" );
is( $sample->get_fwd_mer_len(), $fwd_mer_len, "get_fwd_mer_len()" );
is( $sample->get_fwd_max_shifts(), $fwd_max_shifts, "get_fwd_max_shifts()" );
is( $sample->get_min_con_depth(), $min_con_depth, "get_min_con_depth()" );
is( $sample->get_diginorm_max(), $DIGI_MAX, "get_diginorm_max()" );
is( $sample->get_output_dir(), $outputDir, "get_output_dir()" );

# Test setters
is( $sample->set_id("P2"), 1, "set_id(P2)" );
is( $sample->set_name("sample2"), 1, "set_name(sample2)" );
is( $sample->set_barcode("ATCAG"), 1, "set_barcode()" );
is( $sample->set_fwd_linker("CAGG"), 1, "set_fwd_linker()" );
is( $sample->set_fwd_primer("GTGCCAGC[AC]GCCGCGGTAT"), 1, "set_fwd_primer()" );
is( $sample->set_rev_primer("GGACTAC[AT][ACG]GGGT[AT]TCTAAT"), 1, "set_rev_primer()" );
is( $sample->set_fwd_mer_len(8), 1, "set_fwd_mer_len(8)" );
is( $sample->set_fwd_max_shifts(5), 1, "set_fwd_max_shifts(5)" );
is( $sample->set_min_con_depth(3), 1, "set_min_con_depth(3)" );
is( $sample->set_diginorm_max(5), 1, "set_diginorm_max(5)" );
is( $sample->set_output_dir("/home/scott/output"), 1, "set_output_dir()" );

# Test the setter accuracy
is( $sample->get_id(), "P2", "get_id()" );
is( $sample->get_name(), "sample2", "get_name()" );
is( $sample->get_barcode(), "ATCAG", "get_barcode()" );
is( $sample->get_fwd_linker(), "CAGG", "get_fwd_linker()" );
is( $sample->get_fwd_primer(), "GTGCCAGC[AC]GCCGCGGTAT", "get_fwd_primer()" );
is( $sample->get_rev_primer(), "GGACTAC[AT][ACG]GGGT[AT]TCTAAT", "get_rev_primer()" );
is( $sample->get_fwd_mer_len(), 8, "get_fwd_mer_len()" );
is( $sample->get_fwd_max_shifts(), 5, "get_fwd_max_shifts()" );
is( $sample->get_min_con_depth(), 3, "get_min_con_depth()" );
is( $sample->get_diginorm_max(), 5, "get_diginorm_max()" );
is( $sample->get_output_dir(), "/home/scott/output", "get_output_dir()" );

# Test the forward linker
# Sometimes I will want to leave out the forward linker.  I need to make sure when
# that happens it returns an empty string on a get_fwd_linker call
$sample->set_fwd_linker(undef);
is( $sample->get_fwd_linker(), "", "get_fwd_linker() on undef" );
