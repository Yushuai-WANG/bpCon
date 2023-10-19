#!/usr/bin/perl
use warnings;
use strict;

die "Usage: perl program <in> <out>" if (@ARGV != 10); # usage infromation
my ($list1,$list2,$list3,$list4,$list5,$out,$flag1,$flag2,$flag3,$flag4) = @ARGV;

open (OUT,">$out") or die "Can not open the output_file:$!\n";

my $hash1;
my $hash2;

open (IN,$list2) or die "Can not open the $list1:$!\n";
while(my $line = <IN>){
    chomp $line;
    my @array = split /\t+/,$line;
    my $exp = (($array[1] + $array[2])/2)."\t".(($array[3] + $array[4])/2);
    $hash1->{$array[0]}->{"exp"}=$exp;
    if((($array[1] >= $flag1 && $array[2] >= $flag1) || ($array[3] >= $flag1 && $array[4] >= $flag1)) && ((($array[1] + $array[2])/2 >= $flag2) || (($array[3] + $array[4])/2 >= $flag2))){
        $hash1->{$array[0]}->{"flag"} = 1;
    }
    else{
        $hash1->{$array[0]}->{"flag"} = 0;
    }
}

close IN;

open (IN,$list4) or die "Can not open the $list1:$!\n";
while(my $line = <IN>){
    chomp $line;
    my @array = split /\t+/,$line;
    my $exp = (($array[1] + $array[2])/2)."\t".(($array[3] + $array[4])/2);
    $hash2->{$array[0]}->{"exp"}=$exp;
    if((($array[1] >= $flag3 && $array[2] >= $flag3) || ($array[3] >= $flag3 && $array[4] >= $flag3)) && ((($array[1] + $array[2])/2 >= $flag4) || (($array[3] + $array[4])/2 >= $flag4))){
        $hash2->{$array[0]}->{"flag"} = 1;
    }
    else{
        $hash2->{$array[0]}->{"flag"} = 0;
    }
}

open (IN,$list1) or die "Can not open the $list1:$!\n";
while(my $line = <IN>){
    chomp $line;
    if($line =~ m/^#/){next;}
    my @array = split /\t+/,$line;
    if($array[6] eq "NA"){$array[6] = 1;}
    if($hash1->{$array[0]}->{"flag"} == 1){
        $hash1->{$array[0]}->{"FC"}=$array[2];
        $hash1->{$array[0]}->{"q"}=$array[6];
    }
    else{
        $hash1->{$array[0]}->{"FC"} = 0;
        $hash1->{$array[0]}->{"q"} = 1;
    }
}
close IN;


open (IN,$list3) or die "Can not open the $list1:$!\n";
while(my $line = <IN>){
    chomp $line;
    if($line =~ m/^#/){next;}
    my @array = split /\t+/,$line;
    if($array[6] eq "NA"){$array[6] = 1;}
    $hash2->{$array[0]}->{"FC"}=$array[2];
    $hash2->{$array[0]}->{"q"}=$array[6];
    if($hash2->{$array[0]}->{"flag"} == 1){
        $hash2->{$array[0]}->{"FC"}=$array[2];
        $hash2->{$array[0]}->{"q"}=$array[6];
    }
    else{
        $hash2->{$array[0]}->{"FC"} = 0;
        $hash2->{$array[0]}->{"q"} = 1;
    }
}
close IN;


close IN;

open (IN,$list5) or die "Can not open the $list1:$!\n";
while(my $line = <IN>){
    chomp $line;
    my @array = split /\t+/,$line;
    my $gene =  $array[0]."_".$array[1]."_".$array[2]."_".$array[3];

    if(exists $hash1->{$gene}->{"exp"}){
        print OUT $line,"\t",$hash1->{$gene}->{"exp"},"\t";
    }
    else{
        print OUT $line,"\t",0,"\t",0,"\t";
    }

    if(exists $hash2->{$gene}->{"exp"}){
        print OUT $hash2->{$gene}->{"exp"},"\t";
    }
    else{
        print OUT 0,"\t",0,"\t";
    }

    if(exists $hash1->{$gene}->{"FC"}){
        print OUT $hash1->{$gene}->{"FC"},"\t",$hash1->{$gene}->{"q"},"\t";
    }
    else{
        print OUT 0,"\t",1,"\t";
    }

    if(exists $hash2->{$gene}->{"FC"}){
        print OUT $hash2->{$gene}->{"FC"},"\t",$hash2->{$gene}->{"q"},"\n";
    }
    else{
        print OUT 0,"\t",1,"\n";
    }
}

close IN;
close OUT;

