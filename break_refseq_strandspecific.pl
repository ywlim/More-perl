#!/usr/bin/perl
# break_refseq_strandspecific.pl by Yoong Wearn Lim
# break a raw refseq table into different bed files
# this version takes strand into consideration

use strict; use warnings;

open(IN, "/data/wearn/Perl/refseq_table.txt") or die "error opening refseq_table.txt\n";

open(PROMOTER_P, ">tss_plus.bed") or die "can't open tss_plus.bed\n";
open(FIRSTEXON_P, ">firstexon_plus.bed") or die "can't open firstexon_plus.bed\n";
open(FIVEUTR_P, ">5utr_plus.bed") or die "can't open 5utr_plus.bed\n";
open(THREEUTR_P, ">3utr_plus.bed") or die "can't open 3utr_plus.bed\n";
open(EXON_P, ">exon_plus.bed") or die "can't open exon_plus.bed\n";
open(INTRON_P, ">intron_plus.bed") or die "can't open intron_plus.bed\n";
open(GENE_P, ">gene_plus.bed") or die "can't open gene_plus.bed\n";
open(PROMOTER_M, ">tss_minus.bed") or die "can't open tss_minus.bed\n";
open(FIRSTEXON_M, ">firstexon_minus.bed") or die "can't open firstexon_minus.bed\n";
open(FIVEUTR_M, ">5utr_minus.bed") or die "can't open 5utr_minus.bed\n";
open(THREEUTR_M, ">3utr_minus.bed") or die "can't open 3utr_minus.bed\n";
open(EXON_M, ">exon_minus.bed") or die "can't open exon_minus.bed\n";
open(INTRON_M, ">intron_minus.bed") or die "can't open intron_minus.bed\n";
open(GENE_M, ">gene_minus.bed") or die "can't open gene_minus.bed\n";

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
		if ($strand eq "+")
		{
			print GENE_P "$chr	$txnstart	$txnend\n";
		}
		if ($strand eq "-")
                {
                        print GENE_M "$chr	$txnstart	$txnend\n";
                }

		my @estart = split(",", $exonstart);
		my @eend= split(",", $exonend);

		for (my $i = 0; $i < @estart; $i++)
		{			
			if (($strand eq "+") && ($i == 0))	{print FIRSTEXON_P "$chr	$estart[$i]	$eend[$i]\n"}	
			elsif(($strand eq "-") && ($i == @estart - 1))	{print FIRSTEXON_M "$chr	$estart[$i]	$eend[$i]\n"}
			else
			{
				if ($strand eq "+")
				{	
					print EXON_P "$chr	$estart[$i]	$eend[$i]\n";
				}
			
				elsif ($strand eq "-")     
                                {
                                        print EXON_M "$chr	$estart[$i]	$eend[$i]\n";
                                }
			}

			# print intron
			if ($i < @estart - 1)	# because there is no intron after the last exon
			{
				if (($estart[$i + 1] - 1) > ($eend[$i] + 1))	# make sure it's not negative
				{
					if ($strand eq "+")     
                                	{
						print INTRON_P "$chr	", $eend[$i] + 1, "	", $estart[$i + 1] - 1 , "\n";
					}

					elsif ($strand eq "-")  
                                        {
                                                print INTRON_M "$chr	", $eend[$i] + 1, "	", $estart[$i + 1] - 1 , "\n";
                                        }
				}
			}
		}

		if ($strand eq "+")
		{
		
			print PROMOTER_P "$chr	", $txnstart - 2000, "	", $txnstart + 2000, "\n";

			# make sure 5 and 3'UTR exists before printing them
			if (($cdsstart - 1) > $txnstart)
			{		
				print FIVEUTR_P "$chr	$txnstart	", $cdsstart - 1, "\n";
			}

			if ($txnend > ($cdsend + 1))
			{
				print THREEUTR_P "$chr	", $cdsend + 1, "	$txnend\n";
			}
		}

		elsif ($strand eq "-")
		{
			print PROMOTER_M "$chr	", $txnend - 2000, "	", $txnend + 2000, "\n";
			if ($txnend > ($cdsend + 1))
			{
				print FIVEUTR_M "$chr	", $cdsend + 1, "	$txnend\n";
			}
			if (($cdsstart - 1) > $txnstart)
			{
				print THREEUTR_M "$chr	$txnstart	", $cdsstart - 1, "\n";
			}
		}
	}

	$seen{$txnstart} = 1;
	$seen{$txnend} = 1;
}

close IN;
close PROMOTER_P;
close PROMOTER_M;
close FIVEUTR_P;
close FIVEUTR_M;
close THREEUTR_P;
close THREEUTR_M;
close EXON_P;
close EXON_M;
close INTRON_P;
close INTRON_M;
close FIRSTEXON_P;
close FIRSTEXON_M;
