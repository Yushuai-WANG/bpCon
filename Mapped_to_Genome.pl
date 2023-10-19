#! /usr/bin/perl -w
#===============================================================================
#  DESCRIPTION:  A pipeline for mapping clean reads to the genome and generate related mapping reports
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
	   -i_G <index>
	   -picard <path>
	   -label <label of the fastq files>
	   -o <summary file>
options:
	-cores CORES; Number of CPU cores to use. Use 0 to auto-detect. Default: 1
	-m|max_processes; Max parallel process. Default: 8
	-f_1|fastq_1; Fastq files for Reads 1
	-f_2|fastq_2; Fastq files for Reads 2
	-label; Label of the fastq files
	-s|strandness; Specify strand-specific information (Default: NONE), other opthions: RF or FR
	-i_G|index_G; Index of the Genome that mapped to
	-g|genome; The genome file
	-o|output; Output summary file
	-dedup[yes|no]; remove PCR duplicates or not. Default: yes
	-picard; path to the picard.jar
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
my $index_G = "";
my $genome = "";
my $output = "";
my $dedup = "yes";
my $picard = "";
my $label = "";

# necessary arguments
GetOptions('f_1|fastq_1=s' => \$fastq_1, 
		   'f_2|fastq_2=s' => \$fastq_2, 
		   'i_G|index_G=s' => \$index_G,
		   'label=s' => \$label,
		   'g|genome=s' => \$genome,
		   'o|output=s' => \$output,
		   'cores=i' => \$cores,
		   'm|max_processes=i' => \$max_processes,
		   's|strandness=s' => \$strandness,
		   'picard=s' => \$picard,
		   'dedup=s' => \$dedup);

if ($fastq_1 eq "" || $fastq_2 eq "" || $index_G  eq "" || $output eq "" || $picard eq "") {&usage;}

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

################################
###### Mapping the genome ######
################################
### mapping
my @commands = ();

for my $i (0..($file_num - 1)){
	if(exists $fastq_files_1[$i] && $fastq_files_2[$i]){
		my $com = "hisat2 --rna-strandness $strandness -p $cores --dta -x $index_G -1 $fastq_files_1[$i] -2 $fastq_files_2[$i] -S $label_file[$i].Genome.sam 2>> $label_file[$i].Genome.err";
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

### count
@commands = ();

for my $i (0..($file_num - 1)){
	my $com = "perl ./Scripts/filter_rRNA_reads.pl $label_file[$i].Genome.sam $label_file[$i].non_Genome > $label_file[$i].non_Genome.log";
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
system("cat *.non_Genome.log | sort -k1,1 > $output.Genome.summary");
system("rm *.non_Genome.log");
system("rm *.Genome.err");
system("rm *.non_Genome");

### sam to bam
@commands = ();

for my $i (0..($file_num - 1)){
	my $com = "samtools view -@ $cores -h -b -S -T $genome $label_file[$i].Genome.sam > $label_file[$i].Genome.bam 2>> sam2bam.log";
	push @commands, $com;
}

# running log
print " ---------  Now running on:  ---------\n";
print join("\n",@commands),"\n";

for my $cmd (@commands) {
    $pm->start and next; system($cmd); $pm->finish;
}

$pm->wait_all_children;

### remove PCR duplicates
if($dedup eq "yes"){
	# sort
	@commands = ();

	for my $i (0..($file_num - 1)){
		my $com = "sambamba sort -n -t $cores -o $label_file[$i].Genome.sort.bam $label_file[$i].Genome.bam";
		push @commands, $com;
	}
	
	# running log
	print " ---------  Now running on:  ---------\n";
	print join("\n",@commands),"\n";
	
	for my $cmd (@commands) {
	    $pm->start and next; system($cmd); $pm->finish;
	}
	
	$pm->wait_all_children;

	# REMOVE_DUPLICATES
	@commands = ();

	for my $i (0..($file_num - 1)){
		my $com = "java -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -Xmx8000M -jar $picard MarkDuplicates REMOVE_DUPLICATES=true I=$label_file[$i].Genome.sort.bam O=$label_file[$i].Genome.dedup.bam M=$label_file[$i].dedup.log";
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
	system(`grep "Unknown" *.dedup.log > $output.dedup.summary`);

}
elsif($dedup eq "no"){
	for my $i (0..($file_num - 1)){
		system("cp $label_file[$i].Genome.bam $label_file[$i].Genome.dedup.bam");
	}
}

### retain primary alignment only
@commands = ();

for my $i (0..($file_num - 1)){
	my $com = "samtools view -b -F 4 -F 8 -F 256 $label_file[$i].Genome.dedup.bam > $label_file[$i].Genome.mapped.bam";
	push @commands, $com;
}

# running log
print " ---------  Now running on:  ---------\n";
print join("\n",@commands),"\n";

for my $cmd (@commands) {
    $pm->start and next; system($cmd); $pm->finish;
}

$pm->wait_all_children;

# forward strand, include 81,83 (and 161,163); reverse strand, include 97,99 (and 145,147)
@commands = ();

for my $i (0..($file_num - 1)){
	my $com1 = "samtools view -f 81 -F 32 $label_file[$i].Genome.mapped.bam > $label_file[$i].forward.sam";
	my $com2 = "samtools view -f 97 -F 16 $label_file[$i].Genome.mapped.bam > $label_file[$i].reverse.sam";
	push @commands, $com1;
	push @commands, $com2;
}

# running log
print " ---------  Now running on:  ---------\n";
print join("\n",@commands),"\n";

for my $cmd (@commands) {
    $pm->start and next; system($cmd); $pm->finish;
}

$pm->wait_all_children;

@commands = ();

for my $i (0..($file_num - 1)){
	my $com1 = "samtools view -f 161 -F 16 $label_file[$i].Genome.mapped.bam >> $label_file[$i].forward.sam";
	my $com2 = "samtools view -f 145 -F 32 $label_file[$i].Genome.mapped.bam >> $label_file[$i].reverse.sam";
	push @commands, $com1;
	push @commands, $com2;
}

# running log
print " ---------  Now running on:  ---------\n";
print join("\n",@commands),"\n";

for my $cmd (@commands) {
    $pm->start and next; system($cmd); $pm->finish;
}

$pm->wait_all_children;

# total mapped reads
@commands = ();

for my $i (0..($file_num - 1)){
	my $com = "cat $label_file[$i].forward.sam $label_file[$i].reverse.sam > $label_file[$i].total.sam";
	push @commands, $com;
}

# running log
print " ---------  Now running on:  ---------\n";
print join("\n",@commands),"\n";

for my $cmd (@commands) {
    $pm->start and next; system($cmd); $pm->finish;
}

$pm->wait_all_children;

# count mapped reads
@commands = ();

for my $i (0..($file_num - 1)){
	my $com = "perl ./Scripts/filter_rRNA_reads.pl $label_file[$i].total.sam $label_file[$i].total.sam.non_Genome > $label_file[$i].total.sam.non_Genome.log";
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
system("cat *.non_Genome.log | sort -k1,1 > $output.mapped.summary");
system("rm *.non_Genome.log");

### sam 2 bam
@commands = ();

for my $i (0..($file_num - 1)){
	my $com = "samtools view -@ $cores -h -b -S -T $genome $label_file[$i].total.sam > $label_file[$i].total.bam 2>> sam2bam.log";
	push @commands, $com;
}
for my $i (0..($file_num - 1)){
	my $com = "samtools view -@ $cores -h -b -S -T $genome $label_file[$i].forward.sam > $label_file[$i].forward.bam 2>> sam2bam.log";
	push @commands, $com;
}
for my $i (0..($file_num - 1)){
	my $com = "samtools view -@ $cores -h -b -S -T $genome $label_file[$i].reverse.sam > $label_file[$i].reverse.bam 2>> sam2bam.log";
	push @commands, $com;
}

# running log
print " ---------  Now running on:  ---------\n";
print join("\n",@commands),"\n";

for my $cmd (@commands) {
    $pm->start and next; system($cmd); $pm->finish;
}

$pm->wait_all_children;

# sort
@commands = ();

for my $i (0..($file_num - 1)){
	my $com = "samtools sort -@ $cores -O BAM -o $label_file[$i].total.sort.bam $label_file[$i].total.bam";
	push @commands, $com;
}
for my $i (0..($file_num - 1)){
	my $com = "samtools sort -@ $cores -O BAM -o $label_file[$i].forward.sort.bam $label_file[$i].forward.bam";
	push @commands, $com;
}
for my $i (0..($file_num - 1)){
	my $com = "samtools sort -@ $cores -O BAM -o $label_file[$i].reverse.sort.bam $label_file[$i].reverse.bam";
	push @commands, $com;
}

# running log
print " ---------  Now running on:  ---------\n";
print join("\n",@commands),"\n";

for my $cmd (@commands) {
    $pm->start and next; system($cmd); $pm->finish;
}

$pm->wait_all_children;

# index
@commands = ();

for my $i (0..($file_num - 1)){
	my $com = "samtools index $label_file[$i].total.sort.bam $label_file[$i].total.sort.bam.bai 2>> bam.log";
	push @commands, $com;
}
for my $i (0..($file_num - 1)){
	my $com = "samtools index $label_file[$i].forward.sort.bam $label_file[$i].forward.sort.bam.bai 2>> bam.log";
	push @commands, $com;
}
for my $i (0..($file_num - 1)){
	my $com = "samtools index $label_file[$i].reverse.sort.bam $label_file[$i].reverse.sort.bam.bai 2>> bam.log";
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
print "Cleaning ...\n";
system("rm *.sam");
system("rm *.log");
