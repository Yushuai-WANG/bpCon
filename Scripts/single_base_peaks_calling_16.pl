#!/usr/bin/perl
use warnings;
use strict;

die "Usage: perl program <in> <out>" if (@ARGV != 5); # usage infromation
my ($list1,$list2,$out,$flag,$flag2) = @ARGV;

open (IN,$list1) or die "Can not open the $list1:$!\n";
open (OUT,">$out") or die "Can not open the output_file:$!\n";

my $hash;
while(my $line = <IN>){
    chomp $line;
    my @array = split /\t+/,$line;

    foreach my $site (($array[1] - $flag)..($array[1] + $flag)){
        $hash->{$array[0]}->{$site}->{$array[3]}->{"start"}->{$array[4]} = $array[7];
    }
    foreach my $site (($array[2] - $flag)..($array[2] + $flag)){
        $hash->{$array[0]}->{$site}->{$array[3]}->{"end"}->{$array[4]} = $array[7];
    }

}
close IN;

open (IN,$list2) or die "Can not open the $list1:$!\n";
my $line = <IN>;
chomp $line;
my @array = split /\t+/,$line;

my $chr_temp = $array[0];
my $start_temp = $array[1];
my $end_temp = $array[2];
my $gene_temp = $array[3];
my $strand_temp = $array[5];
my $count = 1;
my $len = $end_temp - $start_temp + 1;
my $bed = 0;
my $filter = $len;

while(my $line = <IN>){
    chomp $line;
    my @array = split /\t+/,$line;

    my $chr = $array[0];
    my $start = $array[1];
    my $end = $array[2];
    my $gene = $array[3];
    my $strand = $array[5];

    if($gene eq $gene_temp){
        my $near = 0;

        if(exists $hash->{$chr}->{$end_temp}->{$gene}->{"end"} && exists $hash->{$chr}->{$start}->{$gene}->{"start"}){
            foreach my $trans (keys %{$hash->{$chr}->{$end_temp}->{$gene}->{"end"}}){
                my $Fir = $hash->{$chr}->{$end_temp}->{$gene}->{"end"}->{$trans};
                
                if(exists $hash->{$chr}->{$start}->{$gene}->{"start"}->{$trans}){
                    my $Sec = $hash->{$chr}->{$start}->{$gene}->{"start"}->{$trans};

                    if(abs($Sec - $Fir) == 1){
                        $near = 1;
                    }
                }
            }
        }

        if($near == 1){
            $count += 1;
            my $tmp_len = $end - $start + 1;
            $len = $len.",".$tmp_len;
            $filter += $tmp_len;
            my $tmp_start = $start - $start_temp - 1;
            $bed = $bed.",".$tmp_start;

            $chr_temp = $chr;
            $start_temp = $start_temp;
            $end_temp = $end;
            $strand_temp = $strand;
            $gene_temp = $gene;

            $near = 0;
        }
        else{
            if($filter >= $flag2){
                print OUT $chr_temp,"\t",$start_temp,"\t",$end_temp,"\t",$chr_temp."_".$start_temp."_".$end_temp."_".$gene_temp,"\t",$filter,"\t",
                          $strand_temp,"\t",$start_temp,"\t",$end_temp,"\t",0,"\t",$count,"\t",$len,"\t",$bed,"\n";
            }
            $count = 1;
            $len = $end - $start + 1;
            $bed = 0;
            $filter = $len;
            $chr_temp = $chr;
            $start_temp = $start;
            $end_temp = $end;
            $strand_temp = $strand;
            $gene_temp = $gene;
        }
    }
    else{
        if($filter >= $flag2){
            print OUT $chr_temp,"\t",$start_temp,"\t",$end_temp,"\t",$chr_temp."_".$start_temp."_".$end_temp."_".$gene_temp,"\t",$filter,"\t",
                      $strand_temp,"\t",$start_temp,"\t",$end_temp,"\t",0,"\t",$count,"\t",$len,"\t",$bed,"\n";
        }
        $count = 1;
        $len = $end - $start + 1;
        $bed = 0;
        $filter = $len;
        $chr_temp = $chr;
        $start_temp = $start;
        $end_temp = $end;
        $strand_temp = $strand;
        $gene_temp = $gene;
    }
}

if($filter >= $flag2){
    print OUT $chr_temp,"\t",$start_temp,"\t",$end_temp,"\t",$chr_temp."_".$start_temp."_".$end_temp."_".$gene_temp,"\t",$filter,"\t",
              $strand_temp,"\t",$start_temp,"\t",$end_temp,"\t",0,"\t",$count,"\t",$len,"\t",$bed,"\n";
}

close IN;
close OUT;

# chr2    69507468    69527266    peak_15 1193    -   69507468    69527266    0   7   111,230,279,287,155,80,51   0,1762,7002,11485,13365,17564,19747 

