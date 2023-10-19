#! /usr/bin/perl -w
#===============================================================================
#  DESCRIPTION:  A pipeline for remove of cell Contamination such as Mycoplasma, Rickettsia and Acholeplasma and rRNA
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
	   -f_1 <Fastq file list that separated by "," for Reads 1>
	   -f_2 <Fastq file list that separated by "," for Reads 2>
	   -i_C <index>
	   -i_r <index>
	   -label <label of the fastq files>
	   -o <summary file>
options:
	-cores CORES; Number of CPU cores to use. Use 0 to auto-detect. Default: 1
	-m|max_processes; Max parallel process. Default: 8
	-f_1|fastq_1; Fastq files for Reads 1
	-f_2|fastq_2; Fastq files for Reads 2
	-s|strandness; Specify strand-specific information (Default: NONE), other opthions: RF or FR
	-i_C|index_C; Index of the Contamination that mapped to
	-i_r|index_r; Index of the rRNA that mapped to
	-label; Label of the fastq files
	-o|output; Output summary file
USAGE
	exit(1);
}
#======================================== END ==============================================

## option variables with default value
my $cores = 1;
my $max_processes = 8;
my $fastq_1 = "";
my $fastq_2 = "";
my $strandness = "NONE";
my $index_C = "";
my $index_r = "";
my $output = "";
my $label = "";

# necessary arguments
GetOptions('f_1|fastq_1=s' => \$fastq_1, 
		   'f_2|fastq_2=s' => \$fastq_2, 
		   'i_C|index_C=s' => \$index_C, 
		   'i_r|index_r=s' => \$index_r, 
		   'o|output=s' => \$output,
		   'label=s' => \$label,
		   'cores=i' => \$cores,
		   'm|max_processes=i' => \$max_processes,
		   's|strandness=s' => \$strandness);

if ($fastq_1 eq "" || $fastq_2 eq "" || $index_C  eq "" || $index_r eq "" || $output eq "" || $label eq "") {&usage;}

my @fastq_files_1 = ();
my @fastq_files_2 = ();
my @label_file = ();
if (defined $fastq_1) {@fastq_files_1 = split(',', $fastq_1);}
if (defined $fastq_2) {@fastq_files_2 = split(',', $fastq_2);}
if (defined $label) {@label_file = split(',', $label);}

# max_processes
my $pm = Parallel::ForkManager->new($max_processes);
# file_num
my $file_num = scalar(@fastq_files_1);

##################################
###### remove Contamination ######
##################################
# mapping
my @commands = ();

for my $i (0..($file_num - 1)){
	if(exists $fastq_files_1[$i] && $fastq_files_2[$i]){
		my $com = "hisat2 --rna-strandness $strandness -p $cores -a -x $index_C -1 $fastq_files_1[$i] -2 $fastq_files_2[$i] -S $fastq_files_1[$i].Contamination.sam 2>> $fastq_files_1[$i].Contamination.err";
		push @commands, $com;
	}
	else{
		print "Not properly Pair-end reads ?\n";
		last;
	}
}

# running log
print " ---------  Now running on:  ---------\n";
print join("\n",@commands),"\n";

for my $cmd (@commands) {
    $pm->start and next; system($cmd); $pm->finish;
}

$pm->wait_all_children;

# count
@commands = ();

for my $i (0..($file_num - 1)){
	my $com = "perl ./Scripts/filter_rRNA_reads.pl $fastq_files_1[$i].Contamination.sam $fastq_files_1[$i].non_Contamination > $fastq_files_1[$i].non_Contamination.log";
	push @commands, $com;
}

# running log
print " ---------  Now running on:  ---------\n";
print join("\n",@commands),"\n";

for my $cmd (@commands) {
    $pm->start and next; system($cmd); $pm->finish;
}

$pm->wait_all_children;

# merge log file
system("cat *.non_Contamination.log | sort -k1,1 > $output.Contamination.summary");
system("rm *.non_Contamination.log");

# extract clean reads
@commands = ();

for my $i (0..($file_num - 1)){
	my $com = "perl ./Scripts/filter_rRNA_reads2.pl $fastq_files_1[$i].non_Contamination $fastq_files_1[$i] $fastq_files_2[$i] $fastq_files_1[$i].non_C.R1.fastq $fastq_files_1[$i].non_C.R2.fastq";
	push @commands, $com;
}

# running log
print " ---------  Now running on:  ---------\n";
print join("\n",@commands),"\n";

for my $cmd (@commands) {
    $pm->start and next; system($cmd); $pm->finish;
}

$pm->wait_all_children;

# clean temp files
system("rm *.non_Contamination");
system("rm *.Contamination.sam");

#########################
###### remove rRNA ######
#########################
## remove rRNA
# mapping
@commands = ();

for my $i (0..($file_num - 1)){
	if(exists $fastq_files_1[$i] && $fastq_files_2[$i]){
		my $com = "hisat2 --rna-strandness $strandness -p $cores -a -x $index_r -1 $fastq_files_1[$i].non_C.R1.fastq -2 $fastq_files_1[$i].non_C.R2.fastq -S $fastq_files_1[$i].rRNA.sam 2>> $fastq_files_1[$i].rRNA.err";
		push @commands, $com;
	}
}

# running log
print " ---------  Now running on:  ---------\n";
print join("\n",@commands),"\n";

for my $cmd (@commands) {
    $pm->start and next; system($cmd); $pm->finish;
}

$pm->wait_all_children;

# count
@commands = ();

for my $i (0..($file_num - 1)){
	my $com = "perl ./Scripts/filter_rRNA_reads.pl $fastq_files_1[$i].rRNA.sam $fastq_files_1[$i].non_rRNA > $fastq_files_1[$i].non_rRNA.log";
	push @commands, $com;
}

# running log
print " ---------  Now running on:  ---------\n";
print join("\n",@commands),"\n";

for my $cmd (@commands) {
    $pm->start and next; system($cmd); $pm->finish;
}

$pm->wait_all_children;

# merge log file
system("cat *.non_rRNA.log | sort -k1,1 > $output.rRNA.summary");
system("rm *.non_rRNA.log");

# extract clean reads
@commands = ();

for my $i (0..($file_num - 1)){
	my $com = "perl ./Scripts/filter_rRNA_reads2.pl $fastq_files_1[$i].non_rRNA $fastq_files_1[$i].non_C.R1.fastq $fastq_files_1[$i].non_C.R2.fastq $fastq_files_1[$i].non_r.R1.fastq $fastq_files_1[$i].non_r.R2.fastq";
	push @commands, $com;
}

# running log
print " ---------  Now running on:  ---------\n";
print join("\n",@commands),"\n";

for my $cmd (@commands) {
    $pm->start and next; system($cmd); $pm->finish;
}

$pm->wait_all_children;

# rename files
for my $i (0..($file_num - 1)){
	system("mv $fastq_files_1[$i].non_r.R1.fastq $label_file[$i].non_r.R1.fastq");
	system("mv $fastq_files_1[$i].non_r.R2.fastq $label_file[$i].non_r.R2.fastq");
}

# clean temp files
print "Cleaning ...\n";
system("rm *.non_rRNA");
system("rm *.rRNA.sam");
system("rm *.non_C.*.fastq");
system("rm *.Contamination.err *.rRNA.err");
