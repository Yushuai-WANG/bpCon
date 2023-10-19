#!/usr/bin/perl
use warnings;
use strict;

die "Usage: perl program <in> <out>" if (@ARGV != 3); # usage infromation
my ($list1,$list2,$out) = @ARGV;

open (IN,$list1) or die "Can not open the $list1:$!\n";
open (OUT,">$out") or die "Can not open the output_file:$!\n";

my $hash;
while(my $line = <IN>){
    chomp $line;
    my @array = split /\t+/,$line;
    $array[1]=sprintf "%.2f",$array[1];
    $array[2]=sprintf "%.2f",$array[2];
    $array[3]=sprintf "%.2f",$array[3];
    $array[4]=sprintf "%.2f",$array[4];
    $hash->{$array[0]}=$array[1]."\t".$array[2]."\t".$array[3]."\t".$array[4];
}
close IN;

open (IN,$list2) or die "Can not open the $list1:$!\n";
while(my $line = <IN>){
    chomp $line;
    my @array = split /\t+/,$line;
    if($array[6] eq "NA"){$array[6] = 1;}
    if($array[6] <= 1){
        if(exists $hash->{$array[0]}){
            my $name = $array[0];
            $name =~ s/(.+)_(\d+)_(\d+)_(.+?)_(.+)/$1\t$2\t$3\t$4\t$5/g;
            $array[2]=sprintf "%.2f",$array[2];
            print OUT $name,"\t",$array[2],"\t",$array[6],"\t",$hash->{$array[0]},"\n";
        }
    }
    
}

close IN;
close OUT;


