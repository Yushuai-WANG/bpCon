#!/usr/bin/perl
use strict;
use warnings;
 
die "Usage: perl Program_script <input_file> <output_file>\n", if (@ARGV != 3); #Usage information
my ($list,$out,$flag) = @ARGV;

open IN,$list or die "ERROR in IN: $!\n";
open OUT, ">$out" or die "ERROR in OUT: $!\n";

my $hash;
while(my $line = <IN>){
    chomp $line;
    my @array = split /\s+/,$line;
    if(($array[5] >= $flag && $array[6] >= $flag) || ($array[7] >= $flag && $array[8] >= $flag)){
        my $name = $array[0]."_".$array[1]."_".$array[2]."_".$array[3]."_".$array[4];
        if(exists $hash->{$name}){print $line,"\n";}
        $hash->{$name}=$array[5]."\t".$array[6]."\t".$array[7]."\t".$array[8];
    }
}

close IN;

foreach my $name(keys %{$hash}){
    print OUT $name,"\t",$hash->{$name},"\n";
}

close OUT;
