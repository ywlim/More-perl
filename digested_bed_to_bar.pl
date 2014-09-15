#!/usr/bin/perl
# digested_bed_to_bar.pl by YWL 0n 2/14/13
# input compiled digest bed file made by insilicodigest.pl
# convert it into bars that separate each bed entry
# for visualization purpose on UCSC browser

use strict; use warnings;

die "usage: digested_bed_to_bar.pl <compiled bed> \n" unless @ARGV == 1;

open(IN, $ARGV[0]) or die "error opening $ARGV[0]\n";

while (my $line = <IN>)
{
	chomp $line;
	my ($chr, $start, $end) = split("\t", $line);
	print "$chr	$end	", $end+1, "\n";
}

close IN;