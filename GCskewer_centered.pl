#!/usr/bin/perl
# GCskewer_centered.pl by Yoong Wearn Lim
# modified from GCskewer_flexi.pl on 6/16/14
# calculate average GC skew of a list of genes (fasta file) around a user-defined range
# window size = 100

use strict; use warnings;
use FAlite;

die "usage: GCskewer_centered.pl <fasta>\n" unless @ARGV == 1;

my ($filename) = $ARGV[0];
my ($window) = 100;

open(IN, $filename) or die "error reading $filename";
my $fasta = new FAlite(\*IN);

my @skew;
my $w = 0;
my $length;

while (my $entry = $fasta->nextEntry)
{
	my $seq = $entry->seq;
	$length = length($seq);
	#die "length of your seq is $length and your window is $window\n";

	for (my $k = 0; $k < $length - ($window + 1); $k++)					# loop each nucleotide (sliding window)
	{
		my $chunk = substr($seq, $k, $window);					# extract 10bp chunk at a time
		my $G = $chunk =~ tr/Gg/Gg/;							# count number of G
		my $C = $chunk =~ tr/Cc/Cc/;							# count number of C
		my $GCsum = $G + $C;
		if ($GCsum == 0)	{$skew[$w][$k] = 0}				# avoid divided by zero
		else			{$skew[$w][$k] = ($G - $C)/($G + $C)}
	}
	$w++;
}

my $sumskew = 0;

for (my $k = 0; $k < $length - ($window + 1); $k++)						# loop each chunk
{
	for (my $j = 0; $j < $w; $j++)					# loop each gene
	{
		if (!exists $skew[$j][$k])	{$skew[$j][$k] = 0}
		$sumskew += $skew[$j][$k];							# calculate sum of GC skew at chunk k for each gene
	}

	my $aveskew = $sumskew / $w;						# calculate average skew at each chunk
	$sumskew = 0;											# reset sum of skew for the next chunk
	print 0-($length/2) + ($window / 2) + $k,"	$aveskew\n";
}
