#!/usr/bin/perl
use strict;
use warnings;
 
die "Usage: perl Program_script <input_file> <output_file>\n", if (@ARGV != 5); #Usage information
my ($list,$fq1,$fq2,$out1,$out2) = @ARGV;

open IN,$list or die "ERROR in IN: $!\n";

my $hash;
while(my $line = <IN>){
    chomp $line;
	my @array = split /\s+/,$line;
	my $name = "@".$array[0];
	$hash->{$name}=1;
}

close IN;

open IN1,$fq1 or die "ERROR in IN: $!\n";
open IN2,$fq2 or die "ERROR in IN: $!\n";
open OUT1, ">$out1" or die "ERROR in OUT: $!\n";
open OUT2, ">$out2" or die "ERROR in OUT: $!\n";

my $count = 0;
while(my $line1 = <IN1>){
    chomp $line1;
    my @array1 = split /\s+/,$line1;
    if(exists $hash->{$array1[0]}){
    	my $line2 = <IN1>;chomp $line2;
		my $line3 = <IN1>;chomp $line3;
		my $line4 = <IN1>;chomp $line4;

		my $line21 = <IN2>;chomp $line21;
		my @array2 = split /\s+/,$line21;

		my $line22 = <IN2>;chomp $line22;
		my $line23 = <IN2>;chomp $line23;
		my $line24 = <IN2>;chomp $line24;

		#print $line1,"\t",$line2,"\t",$line3,"\t",$line4,"\n"; 
		#print $line21,"\t",$line22,"\t",$line23,"\t",$line24,"\n"; 

		if($array1[0] eq $array2[0]){
			print OUT1 $line1,"\n",$line2,"\n",$line3,"\n",$line4,"\n"; 
			print OUT2 $line21,"\n",$line22,"\n",$line23,"\n",$line24,"\n"; 

			$count++;
		}
		else{
			print $list,"\t",$line1,"\t",$line2,"\n";
		}
    }
    else{
    	<IN1>;<IN1>;<IN1>;
    	<IN2>;<IN2>;<IN2>;<IN2>;
    }
}

# print $list,"\t",$count,"\n";

close IN1;
close IN2;
close OUT1;
close OUT2;
