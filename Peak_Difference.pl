#! /usr/bin/perl -w
#===============================================================================
#  DESCRIPTION:  A pipeline for Call peaks base the enrichment of reads in IP than Input samples
#       AUTHOR:  Yushuai Wang
#      VERSION:  1.0
#      CREATED:  10/18/2023
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
	   -peak_C <Peaks for Control>
	   -peak_T <Peaks for Treatment>
	   -f1_IN <Forward bam file for input1 and input2>
	   -f1_IP <Forward bam file for IP1 and IP2>
	   -f2_IN <Forward bam file for input3 and input4>
	   -f2_IP <Forward bam file for IP3 and IP4>
	   -r1_IN <Backward bam file for input1 and inpu
	   -r1_IP <Backward bam file for IP1 and IP2>
	   -r2_IN <Backward bam file for input3 and inpu
	   -r2_IP <Backward bam file for IP3 and IP4>
	   -c <Mapped reads count for input1, input2, IP1 and IP2>
	   -bed <bed file>
	   -o <output file>
options:
	-peak_C <Peaks for Control>
	-peak_T <Peaks for Treatment>
	-f1_IN <Forward bam file for input1 and input2>
	-f1_IP <Forward bam file for IP1 and IP2>
	-f2_IN <Forward bam file for input3 and input4>
	-f2_IP <Forward bam file for IP3 and IP4>
	-r1_IN <Backward bam file for input1 and inpu
	-r1_IP <Backward bam file for IP1 and IP2>
	-r2_IN <Backward bam file for input3 and inpu
	-r2_IP <Backward bam file for IP3 and IP4>
	-cores CORES; Number of CPU cores to use. Use 0 to auto-detect. Default: 1
	-m|max_processes; Max parallel process. Default: 8
	-c; Mapped reads count for input1, input2, IP1 and IP2
	-bed; bed file
	-adj_p; Adjusted p value. Default: 0.001
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
my $peak_C = "";
my $peak_T = "";
my $f1_IN = "";
my $f1_IP = "";
my $f2_IN = "";
my $f2_IP = "";
my $r1_IN = "";
my $r1_IP = "";
my $r2_IN = "";
my $r2_IP = "";
my $count = "";
my $bed = "";
my $output = "";
my $adj_p = 0.001;
my $fc = 0.58;
my $CPM = 0;
my $resolution = "normal";

# necessary arguments
GetOptions('peak_C=s' => \$peak_C, 
		   'peak_T=s' => \$peak_T,
		   'bed=s' => \$bed,  
		   'r|resolution=s' => \$resolution,
		   'cores=i' => \$cores,
		   'f1_IN=s' => \$f1_IN, 
		   'f1_IP=s' => \$f1_IP, 
		   'f2_IN=s' => \$f2_IN,
		   'f2_IP=s' => \$f2_IP,
		   'r1_IN=s' => \$r1_IN, 
		   'r1_IP=s' => \$r1_IP, 
		   'r2_IN=s' => \$r2_IN,
		   'r2_IP=s' => \$r2_IP,
		   'c|count=s' => \$count, 
		   'adj_p=f' => \$adj_p, 
		   'fc=f' => \$fc, 
		   'CPM=f' => \$CPM, 
		   'o|output=s' => \$output
		   );

if ($f1_IN eq "" || $f1_IP eq "" || $f2_IN  eq "" || $f2_IP eq "" || $peak_C  eq "" || $peak_T eq "" ||
	$r1_IN eq "" || $r1_IP eq "" || $r2_IN  eq "" || $r2_IP eq "" || $output eq "" || $count eq "") {&usage;}


my @f1_IN_files = ();
my @f1_IP_files = ();
my @f2_IN_files = ();
my @f2_IP_files = ();

my @r1_IN_files = ();
my @r1_IP_files = ();
my @r2_IN_files = ();
my @r2_IP_files = ();
if (defined $f1_IN) {@f1_IN_files = split(',', $f1_IN);}
if (defined $f1_IP) {@f1_IP_files = split(',', $f1_IP);}
if (defined $f2_IN) {@f2_IN_files = split(',', $f2_IN);}
if (defined $f2_IP) {@f2_IP_files = split(',', $f2_IP);}

if (defined $r1_IN) {@r1_IN_files = split(',', $r1_IN);}
if (defined $r1_IP) {@r1_IP_files = split(',', $r1_IP);}
if (defined $r2_IN) {@r2_IN_files = split(',', $r2_IN);}
if (defined $r2_IP) {@r2_IP_files = split(',', $r2_IP);}

my @all_count = ();
if (defined $count) {@all_count = split(',', $count);}

####################################  
####### Call peak difference ####### 
####################################
## intersect peaks
# bed12 to bed6
system("perl ./Scripts/single_base_peaks_calling_14.pl $peak_C $peak_C.format IN");
system("perl ./Scripts/single_base_peaks_calling_14.pl $peak_T $peak_T.format IP");
system("cat $peak_C.format $peak_T.format | msort -k f1 -k f4 -k n2 -k n3 -k f7 > $output.all");
system("perl ./Scripts/single_base_peaks_calling_15.pl $output.all $output.all.intersect");

system(`awk '\$5 >= 10 && \$7 == "IN"' $output.all.intersect > $output.all.intersect.IN`);
system(`awk '\$5 >= 10 && \$7 == "IP"' $output.all.intersect > $output.all.intersect.IP`);
system(`awk '\$5 >= 10 && \$7 != "IN" && \$7 != "IP"' $output.all.intersect > $output.all.intersect.mix`);

if ($resolution eq "normal"){
	system("perl ./Scripts/single_base_peaks_calling_16.pl $bed $output.all.intersect.IN  $output.all.intersect.IN.merge  25 50");
	system("perl ./Scripts/single_base_peaks_calling_16.pl $bed $output.all.intersect.IP  $output.all.intersect.IP.merge  25 50");
	system("perl ./Scripts/single_base_peaks_calling_16.pl $bed $output.all.intersect.mix $output.all.intersect.mix.merge 25 50");
}
elsif($resolution eq "high"){
	system("perl ./Scripts/single_base_peaks_calling_16.pl $bed $output.all.intersect.IN  $output.all.intersect.IN.merge  10 25");
	system("perl ./Scripts/single_base_peaks_calling_16.pl $bed $output.all.intersect.IP  $output.all.intersect.IP.merge  10 25");
	system("perl ./Scripts/single_base_peaks_calling_16.pl $bed $output.all.intersect.mix $output.all.intersect.mix.merge 10 25");
}

system("cat $output.all.intersect.IN.merge $output.all.intersect.IP.merge $output.all.intersect.mix.merge > $output.all.merge");
system("rm *.format $output.all $output.all.intersect $output.all.intersect.IN $output.all.intersect.IP $output.all.intersect.mix $output.all.intersect.IN.merge $output.all.intersect.IP.merge $output.all.intersect.mix.merge");

## calculate difference
# format
system(`awk '\$6 == "+"' $output.all.merge > $output.forward`);
system(`awk '\$6 == "-"' $output.all.merge > $output.reverse`);

# count
system(`multiBamSummary BED-file --numberOfProcessors $cores --metagene --BED $output.forward --bamfiles $f1_IN_files[0] $f1_IN_files[1] $f1_IP_files[0] $f1_IP_files[1] $f2_IN_files[0] $f2_IN_files[1] $f2_IP_files[0] $f2_IP_files[1] --outRawCounts $output.forward.tab`);
system(`multiBamSummary BED-file --numberOfProcessors $cores --metagene --BED $output.reverse --bamfiles $r1_IN_files[0] $r1_IN_files[1] $r1_IP_files[0] $r1_IP_files[1] $r2_IN_files[0] $r2_IN_files[1] $r2_IP_files[0] $r2_IP_files[1] --outRawCounts $output.reverse.tab`);
system("cat $output.forward.tab $output.reverse.tab > $output.total.tab");
system("rm $output.forward.tab $output.reverse.tab $output.forward $output.reverse");

# format for DEseq2
system("perl ./Scripts/single_base_peaks_calling_9.pl $output.total.tab $output.all.merge $output.total.tab.test");
system("rm $output.total.tab");

system("Rscript ./Scripts/peak_difference_4_rep2.R $output.total.tab.test $output.IN.DEGs $output.IP.DEGs $output.IN.count $output.IP.count $all_count[0] $all_count[1] $all_count[2] $all_count[3] $all_count[4] $all_count[5] $all_count[6] $all_count[7]");

# format peaks with difference
system("perl ./Scripts/single_base_peaks_calling_17.pl $output.IN.DEGs $output.IN.count $output.IP.DEGs $output.IP.count $output.all.merge $output.all.merge.dif 0 0 0 0");

# cutoff for peaks
system(`awk '((\$13 >= $CPM && \$15 >= $CPM) && ((\$17 > $fc && \$18 < $adj_p) || (\$19 > $fc && \$20 < $adj_p)))' $output.all.merge.dif | awk '\{print \$0,"\t",\$19 - \$17\}' > $output.all.merge.dif.filter`);

######################################################
##### Connect adjacent peaks with the same trend #####
######################################################
if ($resolution eq "normal"){
	# sort
	system("sort -k1,1 -k2,2n -k3,3n $output.all.merge.dif.filter > $output.all.merge.dif.filter.std");
	
	# overlap
	system("perl ./Scripts/single_base_peaks_calling_20.pl $output.all.merge.dif.filter.std $output.all.merge.dif.filter.std.filt 0.3");
	
	# merge
	system("perl ./Scripts/single_base_peaks_calling_21.pl $output.all.merge.dif.filter.std.filt $output.all.merge.dif.filter.std.merge 50 0.15");
	
	# filter peak less than 100bp
	system("awk '\$5 >= 100' $output.all.merge.dif.filter.std.merge > $output.all.merge.dif.filter.std.merge.filt");
	
	## call difference
	system(`awk '\$6 == "+"' $output.all.merge.dif.filter.std.merge.filt > $output.forward`);
	system(`awk '\$6 == "-"' $output.all.merge.dif.filter.std.merge.filt > $output.reverse`);
	system(`multiBamSummary BED-file --numberOfProcessors $cores --metagene --BED $output.forward --bamfiles $f1_IN_files[0] $f1_IN_files[1] $f1_IP_files[0] $f1_IP_files[1] $f2_IN_files[0] $f2_IN_files[1] $f2_IP_files[0] $f2_IP_files[1] --outRawCounts $output.forward.tab`);
	system(`multiBamSummary BED-file --numberOfProcessors $cores --metagene --BED $output.reverse --bamfiles $r1_IN_files[0] $r1_IN_files[1] $r1_IP_files[0] $r1_IP_files[1] $r2_IN_files[0] $r2_IN_files[1] $r2_IP_files[0] $r2_IP_files[1] --outRawCounts $output.reverse.tab`);
	system("cat $output.forward.tab $output.reverse.tab > $output.total.tab");
	system("rm $output.forward.tab $output.reverse.tab $output.forward $output.reverse");
	
	# format for DEseq2
	system("perl ./Scripts/single_base_peaks_calling_9.2.pl $output.total.tab $output.all.merge.dif.filter.std.merge.filt $output.total.tab.test");
	system("rm $output.total.tab");
	
	system("Rscript ./Scripts/peak_difference_4_rep2.R $output.total.tab.test $output.IN.DEGs $output.IP.DEGs $output.IN.count $output.IP.count $all_count[0] $all_count[1] $all_count[2] $all_count[3] $all_count[4] $all_count[5] $all_count[6] $all_count[7]");
	system("perl ./Scripts/single_base_peaks_calling_17.2.pl $output.IN.DEGs $output.IN.count $output.IP.DEGs $output.IP.count $output.all.merge.dif.filter.std.merge.filt $output.all.merge.dif.filter.std.merge.filt.dif 0 0 0 0");
	
	# cutoff for peaks
	system(`awk '((\$13 >= $CPM && \$15 >= $CPM) && ((\$17 > $fc && \$18 < $adj_p) || (\$19 > $fc && \$20 < $adj_p)))' $output.all.merge.dif.filter.std.merge.filt.dif | awk '\{print \$0,"\t",\$19 - \$17\}' > $output.all.merge.dif.trend`);
}
elsif($resolution eq "high"){
	system("cp $output.all.merge.dif.filter $output.all.merge.dif.trend");
}

# Cleaning
system("rm $output.total.tab.test $output.IN.DEGs $output.IP.DEGs $output.IN.count $output.IP.count");
system("rm $output.all.merge $output.all.merge.dif $output.all.merge.dif.filter $output.all.merge.dif.filter.std $output.all.merge.dif.filter.std.filt $output.all.merge.dif.filter.std.filt.merge $output.all.merge.dif.filter.std.filt.merge.filt $output.all.merge.dif.filter.std.merge $output.all.merge.dif.filter.std.merge.filt $output.all.merge.dif.filter.std.merge.filt.dif");
