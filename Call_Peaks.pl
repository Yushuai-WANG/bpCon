#! /usr/bin/perl -w
#===============================================================================
#  DESCRIPTION:  A pipeline for Call peaks base the enrichment of reads in IP than Input samples
#       AUTHOR:  Yushuai Wang
#      VERSION:  1.0
#      CREATED:  10/17/2023
#===============================================================================
my $version = 1.00;

use strict;
use Getopt::Long;
use Parallel::ForkManager;

#======================================= Usage =============================================
sub usage{
	print <<"USAGE";
Version $version
Usage:
	$0 
	   -f_IN <Forward bam file for input1 and input2>
	   -f_IP <Forward bam file for IP1 and IP2>
	   -r_IN <Backward bam file for input1 and input2>
	   -r_IP <Backward bam file for IP1 and IP2>
	   -c <Mapped reads count for input1, input2, IP1 and IP2>
	   -bed <bed file>
	   -o <output file>
options:
	-cores CORES; Number of CPU cores to use. Use 0 to auto-detect. Default: 1
	-m|max_processes; Max parallel process. Default: 8
	-f_IN; Forward bam file for input1 and input2
	-f_IP; Forward bam file for IP1 and IP2
	-r_IN; Backward bam file for input1 and input2
	-r_IP; Backward bam file for IP1 and IP2
	-c; Mapped reads count for input1, input2, IP1 and IP2
	-bed; bed file
	-adj_p; Adjusted p value. Default: 0.05
	-fc; log2 fold change. Default: 0.58
	-CPM; CPM cutoff for the peaks
	-r|resolution [normal | high]; Reresolution for the peak calling, high mode generated more peaks, but also sacrificed true positive rate. Default: normal
	-o; output file
USAGE
	exit(1);
}
#======================================== END ==============================================

## option variables with default value
my $cores = 1;
my $max_processes = 8;
my $f_IN = "";
my $f_IP = "";
my $r_IN = "";
my $r_IP = "";
my $count = "";
my $bed = "";
my $output = "";
my $adj_p = 0.05;
my $fc = 0.58;
my $CPM = 0;
my $resolution = "normal";

# necessary arguments
GetOptions('f_IN=s' => \$f_IN, 
		   'f_IP=s' => \$f_IP, 
		   'r_IN=s' => \$r_IN, 
		   'r_IP=s' => \$r_IP, 
		   'c|count=s' => \$count, 
		   'bed=s' => \$bed, 
		   'adj_p=f' => \$adj_p, 
		   'fc=f' => \$fc, 
		   'CPM=f' => \$CPM, 
		   'o|output=s' => \$output,
		   'cores=i' => \$cores,
		   'r|resolution=s' => \$resolution,
		   'm|max_processes=i' => \$max_processes);

if ($f_IN eq "" || $f_IP eq "" || $r_IN  eq "" || $r_IP eq "" || $output eq "" || $count eq "") {&usage;}


my @f_IN_files = ();
my @f_IP_files = ();
my @r_IN_files = ();
my @r_IP_files = ();
if (defined $f_IN) {@f_IN_files = split(',', $f_IN);}
if (defined $f_IP) {@f_IP_files = split(',', $f_IP);}
if (defined $r_IN) {@r_IN_files = split(',', $r_IN);}
if (defined $r_IP) {@r_IP_files = split(',', $r_IP);}

my @all_count = ();
if (defined $count) {@all_count = split(',', $count);}

# max_processes
my $pm = Parallel::ForkManager->new($max_processes);
# file_num
my $file_num = scalar(@f_IN_files);

########################
###### Call Peaks ######
########################
# sort bed file
system("msort -k f1 -k f4 -k n2 $bed > $bed.std");
system("perl ./Scripts/single_base_peaks_calling_P3.pl $bed.std $bed.merge 0");
system("rm $bed.std");

# get single base reads coverage
system("mkdir $output.base");
system("perl ./Scripts/single_base_peaks_calling_0.pl $bed.merge");
system("mv $bed.merge.* $output.base");

# count file number
my @files = glob("$output.base/*");
my $num = scalar @files;

my @commands = ();
for my $i (1..$num){
	my $com = "perl ./Scripts/single_base_peaks_calling_1.pl $f_IN_files[0] $f_IN_files[1] $f_IP_files[0] $f_IP_files[1] $r_IN_files[0] $r_IN_files[1] $r_IP_files[0] $r_IP_files[1] ./$output.base/$bed.merge.$i ./$output.base/$bed.merge.$i.base";
	push @commands, $com;
}

# running log
print " ---------  Now running on:  ---------\n";
print join("\n",@commands),"\n";

for my $cmd (@commands) {
    $pm->start and next; system($cmd); $pm->finish;
}

$pm->wait_all_children;

# merge file
system("cat ./$output.base/$bed.merge.*.base > $output.base.merge");

# filter site lower than 10 reads
system("perl ./Scripts/single_base_peaks_calling_2.pl $output.base.merge $output.base.filt 10");

# call difference with GLM statistical testing
system("Rscript ./Scripts/peak_difference_2_rep2.R $output.base.filt $output.DEGs $output.count $all_count[0] $all_count[1] $all_count[2] $all_count[3]");

# format
system("perl ./Scripts/single_base_peaks_calling_3.pl $output.count $output.DEGs $output.site");

# filter
if ($resolution eq "normal"){
	system(`awk 'BEGIN{ FS=" ";OFS="\t" }{if((\$7 < 0.001) && ((\$8 >= 0.5 && \$9 >= 0.5) || (\$10 >= 0.5 && \$11 >= 0.5))) print \$0}' $output.site > $output.site2`);
}
elsif($resolution eq "high"){
	system(`awk 'BEGIN{ FS=" ";OFS="\t" }{if(\$7 < 0.00001) print \$0}' $output.site > $output.site2`);
}

# cleaning
system("rm $output.base.merge $output.DEGs $output.count $output.site");

## merge to region
# sort
system("sort -S 40% --parallel=$cores -k1,1 -k5,5 -k2,2n -k3,3n $output.site2 > $output.site.std");
system("rm $output.site2");

# merge sites within 100bp into regions
if ($resolution eq "normal"){
	system("perl ./Scripts/single_base_peaks_calling_6.pl $output.site.std $output.site.std.merge 100");
}
elsif($resolution eq "high"){
	system("perl ./Scripts/single_base_peaks_calling_6.pl $output.site.std $output.site.std.merge 50");
}

# fileter small regions (less than 5bp)
system(`awk 'BEGIN{ FS=" ";OFS="\t" }{if(\$3 - \$2 >= 5) print \$0}' $output.site.std.merge > $output.site.std.merge2`);

# sort again
system("msort -k f1 -k f5 -k n2 -k n3 -k f4 $output.site.std.merge2 > $output.site.std.merge.std");
system("rm $output.site.std $output.site.std.merge $output.site.std.merge2");

# merge region into peaks (for exon), tolerant 25bp distance，filter peaks less than 50bp
if ($resolution eq "normal"){
	system("perl ./Scripts/single_base_peaks_calling_4.pl $bed $output.site.std.merge.std $output.site.std.merge.peak 25 50");
}
elsif($resolution eq "high"){
	system("perl ./Scripts/single_base_peaks_calling_4.pl $bed $output.site.std.merge.std $output.site.std.merge.peak 10 25");
}

# if peaks overlapped more than 30%, retain the longer one
system("msort -k f1 -k n2 -k n3 $output.site.std.merge.peak > $output.site.std.merge.peak.std");
system("perl ./Scripts/single_base_peaks_calling_20.pl $output.site.std.merge.peak.std $output.site.std.merge.peak 0.3");

## enrichment of IP among peaks
system(`awk '\$6 == \"+\"' $output.site.std.merge.peak > $output.forward`);
system(`awk '\$6 == \"-\"' $output.site.std.merge.peak > $output.reverse`);

# count reads
system("multiBamSummary BED-file --numberOfProcessors $cores --metagene --BED $output.forward --bamfiles $f_IN_files[0] $f_IN_files[1] $f_IP_files[0] $f_IP_files[1] --outRawCounts $output.forward.tab");
system("multiBamSummary BED-file --numberOfProcessors $cores --metagene --BED $output.reverse --bamfiles $r_IN_files[0] $r_IN_files[1] $r_IP_files[0] $r_IP_files[1] --outRawCounts $output.reverse.tab");
system("cat $output.forward.tab $output.reverse.tab > $output.total.tab");
system("rm $output.forward $output.reverse $output.forward.tab $output.reverse.tab");

# format for DEseq2
system("perl ./Scripts/single_base_peaks_calling_5.pl $output.total.tab $output.site.std.merge.peak $output.total.tab.test");
system("rm $output.total.tab");

# difference
system("Rscript ./Scripts/peak_difference_2_rep2.R $output.total.tab.test $output.DEGs2 $output.count2 $all_count[0] $all_count[1] $all_count[2] $all_count[3]");

# format peaks with difference
system("perl ./Scripts/single_base_peaks_calling_7.pl $output.DEGs2 $output.count2 $output.site.std.merge.peak $output.site.std.merge.peak.dif");
system("rm *.test *.DEGs2 *.count2");

#### cutoff for peaks
system("perl ./Scripts/single_base_peaks_calling_11.pl $output.site.std.merge.peak.dif $output.site.std.merge.peak.raw");
system(`awk '\$18 < $adj_p && \$17 > $fc && ((\$13 >= $CPM && \$14 >= $CPM) || (\$15 >= $CPM && \$16 >= $CPM))' $output.site.std.merge.peak.raw | msort -k f1 -k f4 -k n2 -k n3 > $output.site.std.merge.peak.dif.filter`);
system("rm $output.site.std.merge.std $output.site.std.merge.peak $output.site.std.merge.peak.std $output.site.std.merge.peak.dif");
