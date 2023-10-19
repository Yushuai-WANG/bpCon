#!/usr/bin/perl
use strict;
use warnings;

die "Usage: perl Program_script <genome> <out>\n", if (@ARGV != 3); #Usage information
my ($list,$out,$flag) = @ARGV;

open IN,$list or die "ERROR in IN: $!\n";
open (OUT,">$out") or die "Can not open the output_file:$!\n";

while(my $line = <IN>){
    chomp $line;
    my @array = split /\t+/,$line;
    
    # $array[3] =~ m/(.+)_(\d+)_(\d+)_(.+)/;
    # my $gene = $4;
    # $array[3] = $gene;

    my @len = split /,/,$array[10];
    my @start = split /,/,$array[11];

    my $i = 0;
    while($i < $array[9]){
        print OUT $array[0],"\t",$array[1] + $start[$i],"\t",$array[1] + $start[$i] + $len[$i],"\t",$array[3],"\t",$len[$i],"\t",$array[5],"\t",$flag,"\n";
        $i++;
    }
}

close IN;
close OUT;

