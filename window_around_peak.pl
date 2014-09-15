#!/usr/bin/perl
# window_around_peak.pl by YWL
# make +-5000 windows around peaks in a bed file

use strict; use warnings;
open (IN, $ARGV[0]) or die "error opening file\n";
while (my $line = <IN>)
{
    chomp $line;
    my @stuff = split("\t", $line);
    my $chr = $stuff[0];
    my $start = $stuff[1];
    my $end = $stuff[2];
    my $mid = int((($end - $start) / 2) + $start);
    my $left = $mid - 5000;
    my $right = $mid + 5000;
    print "$chr\t$left\t$right\n";
}

close IN;