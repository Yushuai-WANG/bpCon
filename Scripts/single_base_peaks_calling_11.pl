#!/usr/bin/perl
use strict;
use warnings;

die "Usage: perl Program_script <genome> <out>\n", if (@ARGV != 2); #Usage information
my ($list,$out) = @ARGV;

open IN,$list or die "ERROR in IN: $!\n";
open (OUT,">$out") or die "Can not open the output_file:$!\n";

while(my $line = <IN>){
    chomp $line;
    my @array = split /\t+/,$line;
    $array[3] =~ m/(.+)_(\d+)_(\d+)_(.+)/;
    my $gene = $4;
    $array[3] = $gene;
    print OUT join("\t",@array),"\n";
}

close IN;
close OUT;

