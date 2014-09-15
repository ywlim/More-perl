#!/usr/bin/perl
# wig_to_metagene.pl by Yoong Wearn Lim on 9/5/2014
# modified from wig_to_metaplot_strand_specific_lowmem.pl
# input wig file and generate gene metaplot for the bed region
# current script breaks gene into 10000 parts
# and also take gene length up and down stream regions
# so total part is 30000
# temp files are not deleted for debug purpose
# this script only works for strand specific wig files (eg dripc)
# as it will report sense and antisense signals

use strict; use warnings;
use Getopt::Std;	# first time using get opt! excited!
use FileHandle;

our ($opt_h, $opt_w, $opt_b, $opt_n);
getopts('hw:b:');	# w and b take arguments, thus the : after them

my $usage = "usage: wig_to_metaplot.pl -w plus.wig,minus.wig -b bed \n";

if (($opt_h) || (@ARGV == 0))
{
	die $usage;
}

if (!$opt_b)	{die "-b not given\n"}
if (!$opt_w)	{die "-w not given\n"}

my @wig = split(",", $opt_w);
my $bed = $opt_b;

open (OUT, ">metagene.txt") or die "can't open outfile\n";

my %val;
my $span = 10;

#### WIG FILE ##########
# q = 0 is plus wig file
# q = 1 is minus wig file
# breaking wig file up
# my $trigger = 0;
# print "Pre-processing wig files...\n";
# for (my $q = 0; $q < 2; $q++)
# {
# 	if ($wig[$q] =~ /\.gz$/)
# 	{
# 		open (WIG, "gunzip -c $wig[$q] |") || die "can't open pipe to $wig[$q]\n";
# 	}
# 	else
# 	{
# 		open (WIG, $wig[$q]) || die "can't open $wig[$q]\n";
# 	}
#
# 	my $chr;
#
# 	while (my $line = <WIG>)
# 	{
# 		chomp $line;
# 		next if (($line !~ m/^variableStep/) and ($line !~ /^\d+/));
# 		die "sorry, fixedStep wig file not supported\n" if ($line =~ m/^fixedStep/);
#
# 		if ($line =~ m/^variableStep/)
# 		{
# 			close TEMP if ($trigger == 1);
# 			($chr, $span) = $line =~ m/^variableStep\schrom=chr(\w+)\sspan=(\d+)/;
# 			# print "chr is $chr and span is $span\n";
# 			print "Breaking chr$chr...\n";
# 			my $wig_out = $chr . "_" . $q . "_wig.temp";
# 			open (TEMP, ">$wig_out") or die "error writing to $wig_out\n";
# 		}
# 		elsif ($line =~ m/\d+\t\d+/)
# 		{
# 			print TEMP "$line\n";
# 			$trigger = 1;
# 		}
# 	}
# }

# breaking bed files up
# my %FH;
# my @chromosome;
# my %seen;
#
# open (BED, $bed) or die "can't open $bed bed file\n";
#
# print "Pre-processing bed file: $bed\n";
# while (my $line = <BED>)
# {
#     chomp $line;
#     my ($chrom) = $line =~ m/^chr(\w+)/;
#     if (!defined $seen{$chrom})
#     {
#         push (@chromosome, $chrom); # these are the chromosomes with bed values
#         $seen{$chrom} = $chrom;
#     }
#
#     if (!exists $FH{$chrom})	# see chrom for the first time, open a new file
#     {
#
#         $FH{$chrom} = new FileHandle;
#         $FH{$chrom}->open(">$chrom\_bed.temp");
#         $FH{$chrom}->print("$line\n");
#     }
#
#     else
#     {
#         $FH{$chrom}->print("$line\n");
#     }
# }
# close BED;

#my @chromosome = qw(21);
my @chromosome = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);



my %depth; my %count; my %winval; my %wincount;
my $windowsize;
my $k_int;

for (my $c = 0; $c < @chromosome; $c++)
{
	print "Loading chr$chromosome[$c] miniwig into memory...\n";
	# load miniwig into memory
	for (my $q = 0; $q < 2; $q++)
	{
		open (MINIWIG, "$chromosome[$c]\_$q\_wig.temp") or last;    # last because that wig for this chro does't exist;

		while (my $line = <MINIWIG>)
		{
			chomp $line;
			my ($position, $value) = $line =~ m/(\d+)\t(\S+)/;

			for (my $i = $position; $i <= ($position + $span); $i++)
			{
				if ($q == 0)
				{
					if (!defined $val{plus}{$chromosome[$c]}{$i})
					{
						$val{plus}{$chromosome[$c]}{$i} = 0;
					}
					$val{plus}{$chromosome[$c]}{$i} += $value;
				}
				elsif ($q == 1)
				{
					if (!defined $val{minus}{$chromosome[$c]}{$i})
					{
						$val{minus}{$chromosome[$c]}{$i} = 0;
					}
					$val{minus}{$chromosome[$c]}{$i} += $value;
				}
			}
		}
		close MINIWIG;
	}

	# process the same chrom bed file

    # need to close all bed filehandles
#     if (exists $FH{$chromosome[$c]})
#     {
#         $FH{$chromosome[$c]}->close;
#     }

    print "Working on chr$chromosome[$c]...\n";
    open (MINIBED, "$chromosome[$c]\_bed.temp") or last;

    while (my $line2 = <MINIBED>)
    {
        chomp $line2;
        my ($chro, $start, $end, $strand) = $line2 =~ m/^chr(\w+)\t(\d+)\t(\d+)\t\S+\t\S+\t(\S)/;
        #print "$chro	$start	$end\n";

        # make one gene length window up and downstream of gene
        my $genelength = $end - $start;
        my $up = $start - $genelength;
        my $down = $end + $genelength;
        $windowsize = $genelength / 10000; # 10000 parts in gene, 10000 up and 10000 down, total 30000 parts in my graph

        my $w = 1;
        # loop the windows in the 30000 parts
        # non overlapping windows
        for (my $k = $up; $k <= $down - $windowsize; $k += $windowsize)
        {
            # loop within the individual window
            #print "test $k $windowsize\n";
            $k_int = int($k);
            for (my $j = 0; $j <= $windowsize; $j++)
            {
                # get value from plus wig file
                #print "k is $k and j is $j and windowsize is $windowsize\n";
                if (exists $val{plus}{$chromosome[$c]}{$k_int+$j})
                {
                    #print "YA plus: chr$c: j is $j and value is $val{plus}{$chromosome[$c]}{$j}\n";
                    $depth{plus} += $val{plus}{$chromosome[$c]}{$k_int+$j};
                    $count{plus}++;
                }

                # get value from minus wig file
                if (exists $val{minus}{$chromosome[$c]}{$k_int+$j})
                {
                    #print "YA minus: chr$c: j is $j and value is $val{minus}{$chromosome[$c]}{$j}\n";
                    $depth{minus} += $val{minus}{$chromosome[$c]}{$k_int+$j};
                    $count{minus}++;
                }
            }

            # the window gets a value for this gene
            if ($strand eq "+")
            {
                if (defined $count{plus} && defined $depth{plus})
                {
                    $winval{sense}[$w] += $depth{plus} / $count{plus};
                    $wincount{sense}[$w]++;
                }
                if (defined $count{minus} && defined $depth{minus})
                {
                    $winval{antisense}[$w] += $depth{minus} / $count{minus};
                    $wincount{antisense}[$w]++;
                }
            }

            elsif ($strand eq "-")
            {
                # need to assign this to the last window
                if (defined $count{minus} && defined $depth{minus})
                {
                    $winval{sense}[30000 - $w] += $depth{minus} / $count{minus};
                    $wincount{sense}[30000 - $w]++;
                }
                if (defined $count{plus} && defined $depth{plus})
                {
                    $winval{antisense}[30000 - $w] += $depth{plus} / $count{plus};
                    $wincount{antisense}[30000 - $w]++;
                }
            }

            # restart values for individual window
            undef $depth{plus};
            undef $count{plus};
            undef $depth{minus};
            undef $count{minus};
            $w++;   # shift one window in the 1000 parts
        }
    }
    close MINIBED;
}

# remove miniwig from memory
undef %val;

# print the result for each window
print "Now compiling results...\n";

my $avg_sense; my $avg_antisense;
for (my $w = 1; $w <= 30000; $w++)
{
    if (!defined $wincount{sense}[$w])
    {
        $avg_sense = 0;
    }
    else
    {
        $avg_sense = $winval{sense}[$w] / $wincount{sense}[$w];
    }
    if (!defined $wincount{antisense}[$w])
    {
        $avg_antisense = 0;
    }
    else
    {
        $avg_antisense = $winval{antisense}[$w] / $wincount{antisense}[$w] * -1;    # report negative value for antisense
    }

    print OUT "$w\t$avg_sense\t$avg_antisense\n";
}

close OUT;


__END__

#### RESULTS ############
# print result (average depth at each position h)
# assume that the bed coordinates are centered (eg. +- 500 TSS)

# output file header
my @longname;
for (my $w = 0; $w < @name; $w++)
{
	push (@longname, ($name[$w] . "_sense"));
	push (@longname, ($name[$w] . "_antisense"));
}

my $num_longname = @longname;

open (OUT, ">metagene.txt") or die "can't write to metagene.txt\n";
print "Now printing result\n";
my $header = join("\t", @longname);
print OUT "bp\t$header\n";

my %avg_depth;

my $h_adjusted = 0 - ($k / 2);
for (my $h = 1; $h < $k; $h++)
{
	print OUT $h_adjusted + $h, "\t";
	for (my $f = 0; $f < @bed; $f+=2)	# loop each bed set (plus and minus)
	{
		# get average for plus and minus strand, flipping minus strand coordinate with minus array [-$h]
		# sense: plus ($f) gene with plus wig signal; minus ($f+1) gene with minus wig signal

		# sense
		if ((!defined $count{plus}[$f][$h]) and (!defined $count{minus}[$f+1][-$h]))
		{
			$avg_depth{sense} = "NA";
		}
		elsif (!defined $count{plus}[$f][$h])
		{
			$avg_depth{sense} = $depth{minus}[$f+1][-$h] / $count{minus}[$f+1][-$h];
		}
		elsif (!defined $count{minus}[$f+1][-$h])
		{
			$avg_depth{sense} = $depth{plus}[$f][$h] / $count{plus}[$f][$h];
		}
		else
		{
			$avg_depth{sense} = ($depth{plus}[$f][$h] + $depth{minus}[$f+1][-$h]) / ($count{plus}[$f][$h] + $count{minus}[$f+1][-$h]);
		}

		# antisense (* -1 to get negative values)
		if ((!defined $count{minus}[$f][$h]) and (!defined $count{plus}[$f+1][-$h]))
		{
			$avg_depth{antisense} = "NA";
		}
		elsif (!defined $count{plus}[$f+1][-$h])
		{
			$avg_depth{antisense} = ($depth{minus}[$f][$h] / $count{minus}[$f][$h]) * -1;
		}
		elsif (!defined $count{minus}[$f][$h])
		{
			$avg_depth{antisense} = ($depth{plus}[$f+1][-$h] / $count{plus}[$f+1][-$h]) * -1;
		}
		else
		{
			$avg_depth{antisense} = (($depth{minus}[$f][$h] + $depth{plus}[$f+1][-$h]) / ($count{minus}[$f][$h] + $count{plus}[$f+1][-$h])) * -1;
		}

		print OUT "$avg_depth{sense}\t$avg_depth{antisense}\t";

	}
	print OUT "\n";
}

close OUT;

__END__

############ make R script for graphing result ##############
### R script needs fixing: need to combine bedname with sense or antisense



open (R, ">metaplot.R") or die "can't open metaplot.R\n";
print R "library(ggplot2)\n";
print R "library(reshape)\n";
print R "pdf(file=\"metaplot.pdf\", family=\"Helvetica\", width=12, height=8)\n";
print R "plot<-read.table(\"metaplot.txt\", header=T)\n";
print R "plot.melt <- melt(plot[,c('bp', ";

for (my $w = 0; $w < @longname; $w++)
{
	print R "'$longname[$w]'";
	print R ", " unless ($w == $num_longname - 1);
}

print R ")], id.vars=1)\n";
print R "ggplot(plot.melt, aes(x=bp, y=value, colour=variable, group=variable)) + geom_smooth() + theme_bw() + opts(panel.grid.minor=theme_blank()) + scale_colour_brewer(palette=\"Set1\", name=\"Bed\")\n";

close R;

################ run that R script! ##############

`R --vanilla < metaplot.R`;
`rm *.temp`;
