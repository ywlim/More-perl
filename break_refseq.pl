#!/usr/bin/perl
# break_refseq.pl by Yoong Wearn Lim
# break a raw refseq table into different bed files

use strict; use warnings;

open(IN, "/data/wearn/Perl/refseq_table.txt") or die "error opening refseq_table.txt\n";

open(PROMOTER, ">tss.bed") or die "can't open tss.bed\n";
open(FIRSTEXON, ">firstexon.bed") or die "can't open firstexon.bed\n";
open(FIVEUTR, ">5utr.bed") or die "can't open 5utr.bed\n";
open(THREEUTR, ">3utr.bed") or die "can't open 3utr.bed\n";
open(EXON, ">exon.bed") or die "can't open exon.bed\n";
open(INTRON, ">intron.bed") or die "can't open intron.bed\n";
open(GENE, ">gene.bed") or die "can't open gene.bed\n";
open (TTS, ">tts.bed") or die "can't open tts.bed\n";

my %seen;

while (my $line = <IN>)
{
	chomp $line;
	my ($chr, $strand) = $line =~ m/\S+\t\S+\t(\S+)\t(\S)/;	
	
	# note that txnend is actually txnstart if it's minus strand
	my ($txnstart, $txnend, $cdsstart, $cdsend, $exonstart, $exonend) = $line =~ m/\S+\t\S+\t\S+\t\S\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t\S+\t(\S+)\t(\S+)/;

	# filter for unique start sites only
	if ((($strand eq "+") && (!$seen{$txnstart})) or (($strand eq "-") && (!$seen{$txnend})))
	{

		print GENE "$chr	$txnstart	$txnend\n";
		my @estart = split(",", $exonstart);
		my @eend= split(",", $exonend);

		for (my $i = 0; $i < @estart; $i++)
		{			
			if (($strand eq "+") && ($i == 0))	{print FIRSTEXON "$chr	$estart[$i]	$eend[$i]\n"}	
			elsif(($strand eq "-") && ($i == @estart - 1))	{print FIRSTEXON "$chr	$estart[$i]	$eend[$i]\n"}
			else	{print EXON "$chr	$estart[$i]	$eend[$i]\n"}
			
			# print intron
			if ($i < @estart - 1)	# because there is no intron after the last exon
			{
				if (($estart[$i + 1] - 1) > ($eend[$i] + 1))	# make sure it's not negative
				{
					print INTRON "$chr	", $eend[$i] + 1, "	", $estart[$i + 1] - 1 , "\n";
				}
			}
		}

		if ($strand eq "+")
		{
		
			print PROMOTER "$chr	", $txnstart - 2000, "	", $txnstart + 2000, "\n";
			print TTS "$chr	", $txnend - 2000, "	", $txnend + 2000, "\n";

			# make sure 5 and 3'UTR exists before printing them
			if (($cdsstart - 1) > $txnstart)
			{		
				print FIVEUTR "$chr	$txnstart	", $cdsstart - 1, "\n";
			}

			if ($txnend > ($cdsend + 1))
			{
				print THREEUTR "$chr	", $cdsend + 1, "	$txnend\n";
			}
		}

		elsif ($strand eq "-")
		{
			print PROMOTER "$chr	", $txnend - 2000, "	", $txnend + 2000, "\n";
			print TTS "$chr	", $txnstart - 2000, "	", $txnstart + 2000, "\n";
			if ($txnend > ($cdsend + 1))
			{
				print FIVEUTR "$chr	", $cdsend + 1, "	$txnend\n";
			}
			if (($cdsstart - 1) > $txnstart)
			{
				print THREEUTR "$chr	$txnstart	", $cdsstart - 1, "\n";
			}
		}
	}

	$seen{$txnstart} = 1;
	$seen{$txnend} = 1;
}

close IN;
close PROMOTER;
close FIVEUTR;
close THREEUTR;
close EXON;
close INTRON;
close FIRSTEXON;
