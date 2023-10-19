#!/usr/bin/perl
use strict;
use warnings;
 
die "Usage: perl Program_script <input_file> <output_file>\n", if (@ARGV != 2); #Usage information
my ($list,$out) = @ARGV;

open IN,$list or die "ERROR in IN: $!\n";
open OUT, ">$out" or die "ERROR in OUT: $!\n";

my $hash1;
my $hash2;
while(my $line = <IN>){
	if($line =~ m/^@/){next;}
	
    chomp $line;
	my @array = split /\t+/,$line;
	if($array[2] eq "*"){
		$hash1->{$array[0]}=1;
	}
	elsif($array[2] ne "*"){
		$hash2->{$array[0]}=1;
	}
}

close IN;

foreach my $read (keys %{$hash1}){
	print OUT $read,"\n";
}

close OUT;

print $list,"\t","rRNA:\t",scalar keys %{$hash2},"\t",
				 "non_rRNA:\t",scalar keys %{$hash1},"\t",
				 "total:\t",(scalar keys %{$hash1}) + (scalar keys %{$hash2}),"\t",
				 (scalar keys %{$hash2})/((scalar keys %{$hash1}) + (scalar keys %{$hash2}))*100,"\n";

