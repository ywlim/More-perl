#!/usr/bin/perl
# density_to_wig.pl by YWL on 2/25/13
# convert density format to wig format
# density format looks like this:
# chr1    0       -1
# chr1    25      -1
# chr1    50      -1
# chr1    75      -1
# chr1    100     -1
# chr1    125     -1
# chr1    150     -1

# wig format looks like this:
# variableStep chrom=chr2 span=0
# 300701 12.5
# 300702 12.5
# 300703 12.5
# 300704 12.5

# note that this script was specially written for GSM327666_hES.RING1B.densities.txt
# and assume that the span is 25

use strict; use warnings;

die "usage: density_to_wig.pl <density file>\n" unless @ARGV == 1;
open (IN, $ARGV[0]) or die "error opening $ARGV[0]\n";

my $oldchr;

while (my $line = <IN>)
{
        chomp $line;
        my ($chr, $coor, $value) = split("\t", $line);
        print "variableStep chrom=$chr span=25\n" unless ($chr eq $oldchr);
        print "$coor    $value\n" unless (($value == -1) or ($value == 0));
        $oldchr = $chr;
}

close IN;
