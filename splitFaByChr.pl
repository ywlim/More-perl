#!/usr/bin/perl
# splitFaByChr.pl by Yoong Wearn Lim
# split up a fasta file by chromosome
# I wrote this to split up hg19.fa
# so that I can cat them back together
# in the correct order 1-2-3-4...-X-Y-M
# because BS-seeker is stupid and will change
# my chro name to 0001 0002 etc by order of encounter


use strict; use warnings; use FAlite; use FileHandle;

die "usage: splitFaByChr.pl <fa file>\n" unless (@ARGV == 1);
open (IN, $ARGV[0]) or die "error opening $ARGV[0]\n";

# file handles
my @genome = (1..22, 'X', 'Y', 'M');

my %FH;
for (my $i = 0; $i < @genome; $i++)
{
	$FH{"$genome[$i]"} = new FileHandle;
	$FH{"$genome[$i]"}->open(">$genome[$i]") or die "error opening $genome[$i]\n";
}

my $fasta = new FAlite(\*IN);

while (my $entry = $fasta->nextEntry)
{
	my ($head) = $entry->def;
	my ($chro) = $head =~ m/>chr(\w+)/;
	print "chr is $chro\n";
	my $seq = $entry->seq;

	$FH{"$chro"}->print("$head\n");
	$FH{"$chro"}->print("$seq\n");
}

close IN;

__END__


my $mishere = 0;
while (my $line = <IN>)
{
	chomp $line;
	if ($line =~ m/^>chrM/)
	{
		$mishere = 1;
		next;
	}

	if ($mishere == 0)
	{
		print "$line\n";
	}

	elsif ($mishere == 1)
	{
		if ($line =~ m/^>chr/)
		{
			$mishere = 0;	# reset at next chromosome after M
			print "$line\n";
		}
	}
}

close IN;


__END__
my $fasta = new FAlite(\*IN);

while (my $entry = $fasta->nextEntry)
{
	my $def = $entry->def;
	my $seq = $entry->seq;

	if ($def !~ m/^>chrY/)
	{
		print "$def\n";
		print "$seq\n";
	}
}

close IN;