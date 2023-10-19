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
my $strand_temp = $array[5];
my $len_temp = $array[4];
my $line_temp = $line;

while(my $line = <IN>){
    chomp $line;
    my @array = split /\t+/,$line;

    my $chr = $array[0];
    my $start = $array[1];
    my $end = $array[2];
    my $strand = $array[5];
    my $len = $array[4];
    my $line = $line;
    
    if(($chr eq $chr_temp) && ($strand eq $strand_temp)){
        if(($end_temp - $start)/($end_temp - $start_temp) >= $flag || ($end_temp - $start)/($end - $start) >= $flag){
            if($len_temp >= $len){
                $chr_temp = $chr_temp;
                $start_temp = $start_temp;
                $end_temp = $end_temp;
                $strand_temp = $strand_temp;
                $len_temp = $len_temp;
                $line_temp = $line_temp;
            }
            else{
                $chr_temp = $chr;
                $start_temp = $start;
                $end_temp = $end;
                $strand_temp = $strand;
                $len_temp = $len;
                $line_temp = $line;
            }
        }
        else{
            print OUT $line_temp,"\n";
            $chr_temp = $chr;
            $start_temp = $start;
            $end_temp = $end;
            $strand_temp = $strand;
            $len_temp = $len;
            $line_temp = $line;
        }
    }
    else{
        print OUT $line_temp,"\n";
        $chr_temp = $chr;
        $start_temp = $start;
        $end_temp = $end;
        $strand_temp = $strand;
        $len_temp = $len;
        $line_temp = $line;
    }
}
print OUT $line_temp,"\n";

close IN;
close OUT;

