
use strict;
use warnings;

use MTToolbox::Mismatch;
use Test::More tests => 7;

BEGIN { use_ok( 'MTToolbox::Mismatch'); }

### Create a SampleSummaryInfo object to test
my $quals_str = '?@$AS(.'; 
my $seq = "ATGTGA";
my $mismatch = MTToolbox::Mismatch->new({id => "seq1",
                              desc => "my sequence for testing",
                              seq => $seq,
                              quals_str => $quals_str,
                            });

is( $mismatch->get_id(), "seq1", "get_id()" );
is( $mismatch->get_desc(), "my sequence for testing", "get_desc()" );
is( $mismatch->get_header(), "seq1 my sequence for testing", "get_header()" );
is( $mismatch->get_seq(), $seq, "get_seq()" );
is( $mismatch->get_quals_str(), $quals_str, "get_quals_str()" );

is( $mismatch->get_fastq_str(),
('@seq1 my sequence for testing' . "\n" .
"$seq\n" .
"+\n" .
"$quals_str\n")
, "get_fastq_str()" );
