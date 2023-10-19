#!/usr/bin/perl
use strict;
use warnings;

die "Usage: perl Program_script <genome> <out>\n", if (@ARGV != 3); #Usage information
my ($gtf,$out1,$out2) = @ARGV;

open IN,$gtf or die "ERROR in IN: $!\n";
open (OUT1,">$out1") or die "Can not open the output_file:$!\n";
open (OUT2,">$out2") or die "Can not open the output_file:$!\n";

while(my $line = <IN>){
    chomp $line;
    if($line =~ m/^#/){next;}
    my @array = split /\t+/,$line;
    
    my $chr = $array[0];
    my $region = $array[2];
    my $start = $array[3];
    my $end = $array[4];
    my $annotation = $array[8];

    if($region eq "exon"){
        if($annotation =~ m/gene_id "(.+)"; transcript_id "(.+)"; gene_type .+ gene_name "(.+?)"; transcript_type .+; transcript_name .+; exon_number (.+?); exon_id .+; level .+; transcript_support_level "(.+)"; .+_id .+; .+/){
            my $gene_name = $1;
            my $transcript_id = $2;
            my $gene_symbol = $3;
            my $exon_number = $4;
            my $transcript_support_level = $5;
            if($transcript_support_level ne "" && ($transcript_support_level eq "NA" || $transcript_support_level == 1)){
                print OUT1 $chr,"\t",$start,"\t",$end,"\t",$gene_name,"\t",$transcript_id,"\t",$array[6],"\t",$gene_symbol,"\t",$exon_number,"\n";
            }
        }
        # else{
        #     print $line,"\n";
        # }
        
    }

    if($region eq "gene"){
        $annotation =~ m/gene_id "(.+)"; gene_type .+; gene_name "(.+?)"; level.*/;
        my $gene_name = $1;
        my $gene_symbol = $2;
        print OUT2 $chr,"\t",$start,"\t",$end,"\t",$gene_name,"\t",".","\t",$array[6],"\t",$gene_symbol,"\n";
    }

}

close IN;
close OUT1;
close OUT2;
