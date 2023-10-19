#!/usr/bin/perl
use warnings;
use strict;

die "Usage: perl <program> <in>",if (@ARGV != 10); # usage infromation
my ($bam1,$bam2,$bam3,$bam4,$bam5,$bam6,$bam7,$bam8,$in,$out) = @ARGV;

open IN,$in or die "Can not open IN1: $!\n";
open (OUT,">>$out") or die "Can not open the output_file:$!\n";

my $i = 1;
while(my $line = <IN>){
    chomp $line;
    my @array = split /\s+/,$line;
    my $chr = $array[0];
    my $start = $array[1];
    my $end = $array[2];
    my $gene = $array[3];
    my $strand = $array[5];

    # generate reads count of each site across exons
    my $region = $chr.":".$start.":".$end;
    my $now_out = $out.".".$i;
    if($strand eq "+"){
        `multiBamSummary bins --binSize 1 --bamfiles $bam1 $bam2 $bam3 $bam4 --labels input1 input2 IP1 IP2 --region $region --outRawCounts $now_out`;
    }
    else{
        `multiBamSummary bins --binSize 1 --bamfiles $bam5 $bam6 $bam7 $bam8 --labels input1 input2 IP1 IP2 --region $region --outRawCounts $now_out`;
    }
    $i++;

    open IN2,$now_out or die "Can not open IN2: $!\n";
    while(my $line2 = <IN2>){
        chomp $line2;
        if($line2 =~ m/^#/){next;}
        my @array2 = split /\s+/,$line2;
        print OUT $array2[0],"\t",$array2[1],"\t",$array2[2],"\t",$strand,"\t",$gene,"\t",
                  $array2[3],"\t",$array2[4],"\t",$array2[5],"\t",$array2[6],"\n";
    }
    close IN2;

    `rm $now_out`;
}

close IN;
close OUT;
