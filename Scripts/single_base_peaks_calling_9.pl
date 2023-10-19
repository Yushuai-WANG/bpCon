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
    if($line =~ m/^#/){next;}
    my @array = split /\t+/,$line;
    my ($st,$ed);
    if($array[1] =~ m/,/){
        my @start = split /,/,$array[1];
        $st = $start[0];
    }
    else{
        $st = $array[1];
    }

    if($array[2] =~ m/,/){
        my @end  = split /,/,$array[2];
        $ed = $end[-1];
    }
    else{
        $ed = $array[2]-1;
    }
    
    my $name = $array[0]."\t".$st."\t".$ed;
    $hash->{$name}=$array[3]."\t".$array[4]."\t".$array[5]."\t".$array[6]."\t".$array[7]."\t".$array[8]."\t".$array[9]."\t".$array[10];
    
}
close IN;

open (IN,$list2) or die "Can not open the $list1:$!\n";
while(my $line = <IN>){
    chomp $line;
    my @array = split /\t+/,$line;
    my $name = $array[0]."\t".$array[1]."\t".$array[2];
    my $name1 = $array[0]."\t".$array[1]."\t".($array[2]+1);
    if(exists $hash->{$name}){
        print OUT $array[3],"\t",$hash->{$name},"\n";
    }
    elsif(exists $hash->{$name1}){
        print OUT $array[3],"\t",$hash->{$name1},"\n";
    }
    else{
        print $line,"\n";
    }
    
}

close IN;
close OUT;

