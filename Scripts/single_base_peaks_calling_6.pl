#!/usr/bin/perl
use strict;
use warnings;

die "Usage: perl Program_script <genome> <out>\n", if (@ARGV != 3); #Usage information
my ($list,$out,$flag) = @ARGV;

open IN,$list or die "ERROR in IN: $!\n";
open (OUT,">$out") or die "Can not open the output_file:$!\n";

my $line = <IN>;
chomp $line;
my @array = split /\t+/,$line;

my $chr_temp = $array[0];
my $start_temp = $array[1];
my $end_temp = $array[2];
my $strand_temp = $array[3];
my $gene_temp = $array[4];

while(my $line = <IN>){
    chomp $line;
    my @array = split /\t+/,$line;

    my $chr = $array[0];
    my $start = $array[1];
    my $end = $array[2];
    my $strand = $array[3];
    my $gene = $array[4];
    
    if(($gene eq $gene_temp) && ($chr eq $chr_temp) && ($strand eq $strand_temp)){
        if($start - $end_temp <= $flag){
            $chr_temp = $chr_temp;
            $start_temp = $start_temp;
            $end_temp = $end;
            $strand_temp = $strand_temp;
            $gene_temp = $gene_temp;
        }
        else{
            print OUT $chr_temp,"\t",$start_temp,"\t",$end_temp,"\t",$strand_temp,"\t",$gene_temp,"\n";
            $chr_temp = $chr;
            $start_temp = $start;
            $end_temp = $end;
            $strand_temp = $strand;
            $gene_temp = $gene;
        }
    }
    else{
        print OUT $chr_temp,"\t",$start_temp,"\t",$end_temp,"\t",$strand_temp,"\t",$gene_temp,"\n";
        $chr_temp = $chr;
        $start_temp = $start;
        $end_temp = $end;
        $strand_temp = $strand;
        $gene_temp = $gene;
    }
}
print OUT $chr_temp,"\t",$start_temp,"\t",$end_temp,"\t",$strand_temp,"\t",$gene_temp,"\n";

close IN;
close OUT;

