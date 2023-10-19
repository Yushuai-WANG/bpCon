#!/usr/bin/perl
use warnings;
use strict;

die "Usage: perl program <in> <out>" if (@ARGV != 4); # usage infromation
my ($list1,$list2,$list3,$out) = @ARGV;

open (IN,$list1) or die "Can not open the $list1:$!\n";
open (OUT,">$out") or die "Can not open the output_file:$!\n";

my $hash;
while(my $line = <IN>){
    chomp $line;
    if($line =~ m/^#/){next;}
    my @array = split /\t+/,$line;
    $hash->{$array[0]}->{"FC"}=$array[2];
    $hash->{$array[0]}->{"q"}=$array[6];
    
}
close IN;

open (IN,$list2) or die "Can not open the $list1:$!\n";
while(my $line = <IN>){
    chomp $line;
    my @array = split /\t+/,$line;
    my $exp = $array[1]."\t".$array[2]."\t".$array[3]."\t".$array[4];
    $hash->{$array[0]}->{"exp"}=$exp;
}

close IN;

open (IN,$list3) or die "Can not open the $list1:$!\n";
while(my $line = <IN>){
    chomp $line;
    my @array = split /\t+/,$line;
    if(exists $hash->{$array[3]}){
        print OUT $line,"\t",$hash->{$array[3]}->{"exp"},"\t",$hash->{$array[3]}->{"FC"},"\t",$hash->{$array[3]}->{"q"},"\n";
    }
}

close IN;
close OUT;

