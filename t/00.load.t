
use strict;
use warnings;

use Test::More tests => 22;

BEGIN {
    use_ok( 'MTToolbox' );
    
    use_ok( 'MTToolbox::MTParams' );
    use_ok( 'MTToolbox::MTDepthDist' );
    use_ok( 'MTToolbox::RunSummary');
    
    use_ok( 'MTToolbox::Match' );
    use_ok( 'MTToolbox::Match::SE' );
    use_ok( 'MTToolbox::Match::PE' );
    use_ok( 'MTToolbox::Match::PE::NoOverlap' );
    use_ok( 'MTToolbox::Match::PE::Overlap' );
    
    use_ok( 'MTToolbox::Mismatch' );
    
    use_ok( 'MTToolbox::MoleculeTagCategory' );
    use_ok( 'MTToolbox::MoleculeTagCategory::SE' );
    use_ok( 'MTToolbox::MoleculeTagCategory::PE::NoOverlap' );
    use_ok( 'MTToolbox::MoleculeTagCategory::PE::Overlap' ) ;
    
    use_ok( 'MTToolbox::Sample' );
    use_ok( 'MTToolbox::Sample::SE' );
    use_ok( 'MTToolbox::Sample::PE::NoOverlap' );
    use_ok( 'MTToolbox::Sample::PE::Overlap' );
    
    use_ok( 'MTToolbox::SampleSummary' );
    use_ok( 'MTToolbox::SampleSummary::SE' );
    use_ok( 'MTToolbox::SampleSummary::PE::Overlap' );
    use_ok( 'MTToolbox::SampleSummary::PE::NoOverlap' );
}

# For some reason when I uncomment this it says I fail one of the five above
# tests but it doesn't give any diagnostics.  ???
#diag( "Testing MTToolbox $MTToolbox::VERSION" );
