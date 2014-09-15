#!/usr/bin/perl
# check_fastq.pl by YWL on 5/15/14
# check fastq file to see if it has correct format

use strict; use warnings;

die "usage: <check_fastq.pl> <fastq file>\n" unless @ARGV == 1;
my ($file) = $ARGV[0];

if ($file =~ /\.gz$/)
{
    open(IN, "gunzip -c $file |") || die "can't open pipe to $file";
}

else
{
    open(IN, $file) || die "can't open $file";
}

my $count = 0;
while (my $line = <IN>)
{
    chomp $line;
    my $read = <IN>;
    chomp $read;
    if ($read =~ m/[^ATGCN]/)    # not A T C or G
    {
        die "error at line $count: $read\n";
    }
    my $id = <IN>;
    chomp $id;
    my $quality = <IN>;
    chomp $quality;
    if (length($quality) != length($read))
    {
        die "error at line $count: $quality\n";
    }
    $count++;
}

close IN;

print "Fastq is good!\n";
