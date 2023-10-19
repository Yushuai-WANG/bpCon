#!/usr/bin/perl
use warnings;
use strict;

die "Usage: perl program <in> <out>" if (@ARGV != 4); # usage infromation
my ($list1,$list2,$out,$flag) = @ARGV;

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
my $len_temp = $array[4];
my $count_temp = $array[9];
my $lenarray_temp = $array[10];
my $startarray_temp = $array[11];

while(my $line = <IN>){
    chomp $line;
    my @array = split /\t+/,$line;

    my $chr = $array[0];
    my $start = $array[1];
    my $end = $array[2];
    my $gene = $array[3];
    my $strand = $array[5];
    my $len = $array[4];
    my $count = $array[9];
    my $lenarray = $array[10];
    my $startarray = $array[11];

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

            #print $line,"\n";

            $chr_temp = $chr;
            $start_temp = $start_temp;
            $end_temp = $end;
            $gene_temp = $gene;
            $strand_temp = $strand;
            $len_temp = $len_temp + $len;
            $count_temp = $count_temp + $count;
            $lenarray_temp = $lenarray_temp.",".$lenarray;

            my @st = split /,/,$startarray;
            my $i = 0;
            while ($i < $count){
                my $now = $start + $st[$i] - $start_temp;
                $startarray_temp = $startarray_temp.",".$now;
                $i++;
            }
        }
        else{
            print OUT $chr_temp,"\t",$start_temp,"\t",$end_temp,"\t",$chr_temp."_".$start_temp."_".$end_temp."_".$gene_temp,"\t",$len_temp,"\t",
                          $strand_temp,"\t",$start_temp,"\t",$end_temp,"\t",0,"\t",$count_temp,"\t",$lenarray_temp,"\t",$startarray_temp,"\n";
            
            $chr_temp = $chr;
            $start_temp = $start;
            $end_temp = $end;
            $gene_temp = $gene;
            $strand_temp = $strand;
            $len_temp = $len;
            $count_temp = $count;
            $lenarray_temp = $lenarray;
            $startarray_temp = $startarray;
        }
    }
    else{
        print OUT $chr_temp,"\t",$start_temp,"\t",$end_temp,"\t",$chr_temp."_".$start_temp."_".$end_temp."_".$gene_temp,"\t",$len_temp,"\t",
                  $strand_temp,"\t",$start_temp,"\t",$end_temp,"\t",0,"\t",$count_temp,"\t",$lenarray_temp,"\t",$startarray_temp,"\n";
            
        $chr_temp = $chr;
        $start_temp = $start;
        $end_temp = $end;
        $gene_temp = $gene;
        $strand_temp = $strand;
        $len_temp = $len;
        $count_temp = $count;
        $lenarray_temp = $lenarray;
        $startarray_temp = $startarray;
    }
}

print OUT $chr_temp,"\t",$start_temp,"\t",$end_temp,"\t",$chr_temp."_".$start_temp."_".$end_temp."_".$gene_temp,"\t",$len_temp,"\t",
          $strand_temp,"\t",$start_temp,"\t",$end_temp,"\t",0,"\t",$count_temp,"\t",$lenarray_temp,"\t",$startarray_temp,"\n";
            

close IN;
close OUT;

# chr2    69507468    69527266    peak_15 1193    -   69507468    69527266    0   7   111,230,279,287,155,80,51   0,1762,7002,11485,13365,17564,19747 

