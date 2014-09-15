#!/usr/bin/perl
# decrease_graph_res.pl by YWL
# decrease text file size for plotting

use strict; use warnings;

my $window = 100;
my $count = 1;
my @sum;
open (IN, $ARGV[0]) or die "error opening $ARGV[0]\n";
while (my $line = <IN>)
{
    chomp $line;

   	my @stuff = split("\t", $line);

	for (my $i = 1; $i < @stuff; $i++)
	{
		#print "$i\t$stuff[$i]\n";
		if (!exists $sum[$i])    {$sum[$i] = 0}
		$sum[$i] += $stuff[$i];
    }

	if ($count == $window)
	{
		my $newbp = $stuff[0] - $window/2;
        print "$newbp\t";

		for (my $j = 1; $j < @stuff; $j++)
        {
			my $avg = $sum[$j] / $window;
      		printf("%.2f", $avg);
			print "\t";
			$sum[$j] = 0;
		}
        $count = 1;
		print "\n";
    }
    $count++;
}

close IN;
