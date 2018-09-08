
use strict;
use warnings;

use MTToolbox::MoleculeTagCategory;
use MTToolbox::Match::SE;
use BioUtils::FastqSeq 1.0.0;
use Test::More tests => 6;

BEGIN { use_ok( 'MTToolbox::MoleculeTagCategory'); }

### Create a MTToolbox::MoleculeTagCategory object to test SE sequences
my $mt = "GTCGTGCAG_AGAA";
my $fmt = MTToolbox::MoleculeTagCategory->new({tag => $mt});

is( $fmt->get_tag(), $mt, "get_tag()" );
is( $fmt->get_seq_count(), 0, "get_seq_count()" );
is( $fmt->is_con_defined(), 0, "is_con_defined()" );

is( $fmt->set_is_con_defined(1), 1, "set_is_con_defined(1)" );
is( $fmt->is_con_defined(), 1, "is_con_defined() - yes" );

# NOTE: I test the add_match methods in the child class test files.


