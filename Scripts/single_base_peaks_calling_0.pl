#!/usr/bin/perl
use warnings;
use strict;

die "Usage: perl program <in> <output>",if (@ARGV != 1); # usage infromation
my ($in) = @ARGV;

open FILE,$in or die "IN: $!\n";

my $flag = 1000;
my $i = 1;
my $j = 1;

while(my $line = <FILE>){
	chomp $line;
	if($i > 1000){
		$j++;
		$i = 1;
		close OUT;
	}
	my $out = $in.".".$j;
	open (OUT,">>$out") or die "Can not open the output_file:$!\n";

	if($i <= 1000){
		print OUT $line,"\n";
		$i++;
	}
}
close FILE;
close OUT;

