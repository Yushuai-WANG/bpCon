#!/usr/bin/perl
use strict;
use warnings;

die "Usage: perl Program_script <genome> <out>\n", if (@ARGV != 4); #Usage information
my ($list,$out,$flag,$cutoff) = @ARGV;

open IN,$list or die "ERROR in IN: $!\n";
open (OUT,">$out") or die "Can not open the output_file:$!\n";

my $line = <IN>;
chomp $line;
my @array = split /\t+/,$line;

my $chr_temp = $array[0];
my $start_temp = $array[1];
my $end_temp = $array[2];
my $strand_temp = $array[5];
$array[3] =~ m/.+_\d+_\d+_(.+)/;
my $gene_temp = $1;
my $dif_temp = $array[20];
my @print_temp = ($array[0],$array[1],$array[2],$gene_temp,$array[4],$array[5],$array[6],$array[7],$array[8],$array[9],$array[10],$array[11]);

while(my $line = <IN>){
    chomp $line;
    my @array = split /\t+/,$line;

    my $chr = $array[0];
    my $start = $array[1];
    my $end = $array[2];
    my $strand = $array[5];
    $array[3] =~ m/.+_\d+_\d+_(.+)/;
    my $gene = $1;
    my $dif = $array[20];
    
    if(($gene eq $gene_temp) && ($chr eq $chr_temp) && ($strand eq $strand_temp) && 
        (($dif_temp > $cutoff && $dif > $cutoff) || 
         ($dif_temp < -$cutoff && $dif < -$cutoff) || 
         (($dif_temp >= -$cutoff && $dif_temp < $cutoff) && ($dif >= -$cutoff && $dif < $cutoff))
        )){
        if($start - $end_temp <= $flag){
            $print_temp[0] = $chr_temp;
            $print_temp[1] = $start_temp;
            $print_temp[2] = $end;
            $print_temp[5] = $strand_temp;
            $print_temp[3] = $gene_temp;

            $print_temp[6] = $print_temp[1];
            $print_temp[7] = $print_temp[2];

            $print_temp[4] = $print_temp[4] + $array[4];
            $print_temp[9] = $print_temp[9] + $array[9];
            $print_temp[10] = $print_temp[10].",".$array[10];

            my @st = split /,/,$array[11];
            my $i = 0;
            while ($i < $array[9]){
                my $now = $start + $st[$i] - $print_temp[1];
                $print_temp[11] = $print_temp[11].",".$now;
                $i++;
            }

            $chr_temp = $chr_temp;
            $start_temp = $start_temp;
            $end_temp = $end;
            $strand_temp = $strand_temp;
            $gene_temp = $gene_temp;
            $dif_temp = $dif;
        }
        else{
            print OUT join ("\t", @print_temp),"\n";
            @print_temp = ($array[0],$array[1],$array[2],$gene,$array[4],$array[5],$array[6],$array[7],$array[8],$array[9],$array[10],$array[11]);
            $chr_temp = $chr;
            $start_temp = $start;
            $end_temp = $end;
            $strand_temp = $strand;
            $gene_temp = $gene;
            $dif_temp = $dif;
        }
    }
    else{
        print OUT join ("\t", @print_temp),"\n";
        @print_temp = ($array[0],$array[1],$array[2],$gene,$array[4],$array[5],$array[6],$array[7],$array[8],$array[9],$array[10],$array[11]);
        $chr_temp = $chr;
        $start_temp = $start;
        $end_temp = $end;
        $strand_temp = $strand;
        $gene_temp = $gene;
        $dif_temp = $dif;
    }
}
print OUT join ("\t", @print_temp),"\n";

close IN;
close OUT;

