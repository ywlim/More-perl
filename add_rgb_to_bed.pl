#!/usr/bin/perl
# add_rgb_to_bed.pl by YWL
# add color code to bed with different values
# for visualization on genome browser

use strict; use warnings;

# color scheme from dark to light
my $color0 = "160,160,160"; # grey color for no data point
my $color1 = "25,0,51";     # shades of purple
my $color2 = "51,0,102";
my $color3 = "76,0,153";
my $color4 = "102,0,204";
my $color5 = "127,0,255";
my $color6 = "153,51,255";
my $color7 = "178,102,255";
my $color8 = "204,153,255";
my $color9 = "229,204,255";

my $rgb;

print "track name=\"NT2 methylation\" description=\"NT2 methylation\" itemRgb=\"On\"\n";

open (IN, $ARGV[0]) or die "error opening $ARGV[0]\n";
while (my $line = <IN>)
{
    chomp $line;
    my @stuff = split("\t", $line);

    if ($stuff[3] eq "NA")
    {
        $rgb = $color0;
        $stuff[3] = "-1";   # need to assign a value because genome browser doesn't like NA
    }
    elsif ((0 <= $stuff[3]) && ($stuff[3] <= 10))    {$rgb = $color1}
    elsif ((10 < $stuff[3]) && ($stuff[3] <= 20))    {$rgb = $color2}
    elsif ((20 < $stuff[3]) && ($stuff[3] <= 30))    {$rgb = $color3}
    elsif ((30 < $stuff[3]) && ($stuff[3] <= 40))    {$rgb = $color4}
    elsif ((40 < $stuff[3]) && ($stuff[3] <= 50))    {$rgb = $color4}
    elsif ((50 < $stuff[3]) && ($stuff[3] <= 60))    {$rgb = $color5}
    elsif ((60 < $stuff[3]) && ($stuff[3] <= 70))    {$rgb = $color6}
    elsif ((70 < $stuff[3]) && ($stuff[3] <= 80))    {$rgb = $color7}
    elsif ((80 < $stuff[3]) && ($stuff[3] <= 90))    {$rgb = $color8}
    elsif ((90 < $stuff[3]) && ($stuff[3] <= 100))   {$rgb = $color9}
    else                                           {die "WHAT?\n"}

    print "$stuff[0]\t$stuff[1]\t$stuff[2]\t.\t$stuff[3]\t.\t0\t0\t$rgb\n";
}

close IN;