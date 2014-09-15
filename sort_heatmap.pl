#!/usr/bin/perl
# sort_heatmap.pl by Yoong Wearn Lim
# sort the output of calculate_heatmap.pl
# basically this is a sub function of create_heatmap.pl
# first file will be used to order
# the rest of the files will follow the first file's order

use strict; use warnings; use mitochy;
die "usage: sort_heatmap.pl <TSV file 1> <TSV file 2> <TSV file 3> ..\n" unless (@ARGV >= 1);

my $in; my $out;
open ($in, $ARGV[0]) or die "error opening $ARGV[0]\n";
my $outname = mitochy::getFilename($ARGV[0]);
open ($out, ">", "$outname\.new") or last "error opening $outname\.new\n";
my %row; my %sum; my @data;

while (my $line = <$in>)
{
    chomp $line;
    @data = split("\t", $line);
    my ($gene) = $data[0];
    $row{$gene} = $line;
    $sum{$gene} = 0;
    for (my $i = 1; $i < @data; $i++)
    {
        next if ($data[$i] eq "NA");
        $sum{$gene} += $data[$i];
    }
}

close $in;

my @ordered;
foreach my $gene (sort {$sum{$b} <=> $sum{$a}} keys %sum)
{
    print $out "$row{$gene}\n";
    push(@ordered, $gene);
}
close $out;

# now sort the rest of the files following the first file's order
for (my $j = 1; $j < @ARGV; $j++)
{
    open ($in, $ARGV[$j]) or last "error opening $ARGV[$j]\n";
    my $outname = mitochy::getFilename($ARGV[$j]);
    open ($out, ">", "$outname\.new") or last "error opening $outname\.new\n";
    while (my $line = <$in>)
    {
        chomp $line;
        my ($gene) = $line =~ m/^(\w+)\t/;
        $row{$gene} = $line;
    }
    close $in;

    for (my $k = 0; $k < @ordered; $k++)
    {
        print $out "$row{$ordered[$k]}\n";
    }

    close $out;
}