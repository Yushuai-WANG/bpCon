#!/usr/bin/perl
use strict;
use warnings;

die "Usage: perl Program_script <genome> <out>\n", if (@ARGV != 2); #Usage information
my ($list,$out) = @ARGV;

open IN,$list or die "ERROR in IN: $!\n";
open (OUT,">$out") or die "Can not open the output_file:$!\n";

my $line = <IN>;
chomp $line;
my @array = split /\t+/,$line;

my $chr_temp = $array[0];
my $start_temp = $array[1];
my $end_temp = $array[2];
my $gene_temp = $array[3];
my $strand_temp = $array[5];
my $type_temp = $array[6];

while(my $line = <IN>){
    chomp $line;
    my @array = split /\t+/,$line;

    my $chr = $array[0];
    my $start = $array[1];
    my $end = $array[2];
    my $gene = $array[3];
    my $strand = $array[5];
    my $type = $array[6];
    
    if(($gene eq $gene_temp) && ($chr eq $chr_temp) && ($strand eq $strand_temp) && ($type_temp ne $type)){
        if($start >= $start_temp && $end <= $end_temp){
            if((($start - 1) - $start_temp + 1) != 0){
                print OUT $chr_temp,"\t",$start_temp,"\t",$start-1,"\t",$gene_temp,"\t",($start - 1) - $start_temp + 1,"\t",$strand_temp,"\t",$type_temp,"\n";
            }
            if(($end - $start + 1) != 0){
                print OUT $chr,"\t",$start,"\t",$end,"\t",$gene_temp,"\t",$end - $start + 1,"\t",$strand,"\t",$type_temp."_".$type,"\n";
            }
            if(($end_temp - ($end + 1) + 1) != 0){
                $chr_temp = $chr_temp;
                $start_temp = $end + 1;
                $end_temp = $end_temp;
                $gene_temp = $gene_temp;
                $strand_temp = $strand_temp;
                $type_temp = $type_temp;
            }
            else{
                $line = <IN>;
                if(defined $line){
                    chomp $line;
                    @array = split /\t+/,$line;
                    
                    $chr_temp = $array[0];
                    $start_temp = $array[1];
                    $end_temp = $array[2];
                    $gene_temp = $array[3];
                    $strand_temp = $array[5];
                    $type_temp = $array[6];
                }
            }
        }
        elsif($start >= $start_temp && $start <= $end_temp && $end > $end_temp){
            if((($start - 1) - $start_temp + 1) != 0){
                print OUT $chr_temp,"\t",$start_temp,"\t",$start - 1,"\t",$gene_temp,"\t",($start - 1) - $start_temp + 1,"\t",$strand_temp,"\t",$type_temp,"\n";
            }
            if(($end_temp - $start + 1) != 0){
                print OUT $chr,"\t",$start,"\t",$end_temp,"\t",$gene,"\t",$end_temp - $start + 1,"\t",$strand,"\t",$type_temp."_".$type,"\n";
            }
            if(($end - ($end_temp + 1) + 1) != 0){
                $chr_temp = $chr;
                $start_temp = $end_temp + 1;
                $end_temp = $end;
                $gene_temp = $gene;
                $strand_temp = $strand;
                $type_temp = $type;
            }
        }
        else{
            print OUT $chr_temp,"\t",$start_temp,"\t",$end_temp,"\t",$gene_temp,"\t",$end_temp - $start_temp + 1,"\t",$strand_temp,"\t",$type_temp,"\n";
            $chr_temp = $chr;
            $start_temp = $start;
            $end_temp = $end;
            $gene_temp = $gene;
            $strand_temp = $strand;
            $type_temp = $type;
        }
    }
    else{
        print OUT $chr_temp,"\t",$start_temp,"\t",$end_temp,"\t",$gene_temp,"\t",$end_temp - $start_temp + 1,"\t",$strand_temp,"\t",$type_temp,"\n";
        $chr_temp = $chr;
        $start_temp = $start;
        $end_temp = $end;
        $gene_temp = $gene;
        $strand_temp = $strand;
        $type_temp = $type;
    }
}
print OUT $chr_temp,"\t",$start_temp,"\t",$end_temp,"\t",$gene_temp,"\t",$end_temp - $start_temp + 1,"\t",$strand_temp,"\t",$type_temp,"\n";

close IN;
close OUT;

