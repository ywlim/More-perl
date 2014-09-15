#!/usr/bin/perl
# modified from Ian Korf's de-barcoder.pl by Yoong Wearn Lim on 2/24/12
# this script reads the barcode automatically (requires no knowledge of barcodes)
# print out a file for each barcode it encounters
# prints out just + for the third line to minimize file size
# trim off barcodes and T overhang (7 bases) and quality

use strict; use warnings;
use Getopt::Std;
use FileHandle;

die "usage: de-barcoder_YWL.pl <fastq>\n" unless @ARGV == 1;
my ($file) = $ARGV[0];

open(my $in, $file) or die "error reading $file\n";
my $report = "debarcoder_report.txt";
open (OUT, ">$report") or die "error writing $report\n";

my %count;
my %FH;
my $barcode;

while (my $line1 = <$in>) 
{
	chomp $line1;
	
	my $line2 = <$in>;
	chomp $line2;
	$barcode = substr($line2, 0, 6);	# read what the barcode is
	$line2 = substr($line2, 7);	# trim of 7 bases (barcode + T overhang)
	
	<$in>; # line 3 is not important, just repeating ID
	my $line3 = "+";
	
	my $line4 = <$in>;
	chomp $line4;
	$line4 = substr($line4, 7);	# trim of 7 quality scores

	if (!exists $FH{$barcode})	# see barcode for the first time, open a new file
	{
		$FH{$barcode} = new FileHandle;
		$FH{$barcode}->open(">$barcode\.fastq");
		$FH{$barcode}->print("$line1\n$line2\n$line3\n$line4\n");
	}
		
	else
	{
		$FH{$barcode}->print("$line1\n$line2\n$line3\n$line4\n");
	}
	
	$count{$barcode}++;
}

foreach $barcode (sort {$count{$b} <=> $count{$a}} keys %count)
{
	print OUT "$barcode	$count{$barcode}\n";
}
