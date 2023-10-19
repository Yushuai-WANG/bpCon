# bpCon
bpCon (Base pair connection) is a bioinformatic pipeline for peaks detection suit for kinds of RIP-seq data. Though it currently works, but is still being optimized for a better user experience.

## Table of content
- Background
- Installation and Requirement
- Example and Usage
- Maintainers and Contributing
- License

### Background
Presently, RIP-seq a commonly applied sequencing assay that measures the protein-RNA interactions. Derived from its foundation, sequencing methodologies like MeRIP-seq have been produced to measure the location and abundance of RNA modification sites under specific cellular conditions or specific treatment. In practice, the technique is sensitive to PCR amplification biases commonly found in NGS data, and the efficiency of immunoprecipitation often varies between different IP samples, make it hard to assure high-precision identification of antibody-enriched regions. Here, we developed bpCon (Base pair connection), which provides peaks detection and peak difference estimation for kinds of antibody-based enrichment sequencing assay such as RIP-seq, meRIP-seq, acRIP-seq and so on. In contrast to sliding window-based detection approaches, bpCon performs enrichment analysis on every individual base pair within a given interval and connects consecutive enriched sites to attain high-precision peak regions.
</div>

### Installation and Requirement
bpCon is developed in Perl v5.16.3 and relies on the following software, scripts or packages:
- HISAT2 v2.2.1
- samtools v1.16.1
- sambamba v0.8.2
- Picard MarkDuplicates v2.27.4-SNAPSHOT
- deepTools v3.5.1
- DESeq2 v1.40.2
- R v4.3.1

### Example and Usage
In the current version, bpCon only supports experimental data with two biological replicates, such as two inputs and their corresponding IP groups.

#### 1. Pre-processing the raw sequencing data
Before inputting data into bpCon, removing sequencing adapters and low-quality bases is necessary. Here, we use cutadapt v4.1 for an example. Example code:

>`cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -g ACACTCTTTCCCTACACGACGCTCTTCCGATCT -G GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT --cores 2 -n 2 -q 20,20 --max-n 0.25 --trim-n -m 36 -o XXX.rmAP.R1.fastq -p XXX.rmAP.R2.fastq XXX.R1.fastq XXX.R2.fastq 1> XXX.log 2> XXX.log`

#### 2. Remove contamination
Befor peak calling procedues, bpCon detect and remove cell contamination such as Mycoplasma, Rickettsia and Acholeplasma and rRNA from the sequecing data. The sequencing of contaminations or rRNAs should be prepared by yourself accoring to your studied species and potential contamination types. 

2.1 With these sequences, we first generated an index for them with hisat2-build:

>`hisat2-build cell_contamination.fas cell_contamination `

>`hisat2-build rRNA.fas rRNA `

2.2 Then, pipeline Remove_Contamination.pl of bpCon will remove the reads derived form the potential contaminations automatically and give reports for the proportion of reads that removed:

>`perl Remove_Contamination.pl -cores 6 -m 8 -s RF -i_C ./Cell_Contamination/cell_contamination -i_r ./Cell_Contamination/rRNA -f_1 XXX1.rmAP.R1.fastq,XXX2.rmAP.R1.fastq,XXX3.rmAP.R1.fastq,XXX4.rmAP.R1.fastq,XXX5.rmAP.R1.fastq,XXX6.rmAP.R1.fastq,XXX7.rmAP.R1.fastq,XXX8.rmAP.R1.fastq -f_2 XXX1.rmAP.R2.fastq,XXX2.rmAP.R2.fastq,XXX3.rmAP.R2.fastq,XXX4.rmAP.R2.fastq,XXX5.rmAP.R2.fastq,XXX6.rmAP.R2.fastq,XXX7.rmAP.R2.fastq,XXX8.rmAP.R2.fastq -label XXX1,XXX2,XXX3,XXX4,XXX5,XXX6,XXX7,XXX8 -o XXX`

#### 3. Mapped to genome
After the cleaning steps, we then mapped the clean data to the genome (Using hg38 as example).

3.1 Build genome index using hisat2-build. Using the python scripts included in the HISAT2 package, extract splice-site and exon information from the gene annotation file

>`hisat2_extract_splice_sites.py gencode.v39.annotation.gtf > human.ss`

>`hisat2_extract_exons.py gencode.v39.annotation.gtf > human.exon`

>`hisat2-build -p 12 --ss human.ss --exon human.exon hg38.fa human_tran`

3.2 Alignment the clean reads to the genome.

>`perl Mapped_to_Genome.pl -cores 6 -m 8 -s RF -i_G ./human/human_tran -dedup yes -picard ./software/picard.jar -genome ./human/hg38.fa -f_1 XXX1.non_r.R1.fastq,XXX2.non_r.R1.fastq,XXX3.non_r.R1.fastq,XXX4.non_r.R1.fastq,XXX5.non_r.R1.fastq,XXX6.non_r.R1.fastq,XXX7.non_r.R1.fastq,XXX8.non_r.R1.fastq -f_2 XXX1.non_r.R2.fastq,XXX2.non_r.R2.fastq,XXX3.non_r.R2.fastq,XXX4.non_r.R2.fastq,XXX5.non_r.R2.fastq,XXX6.non_r.R2.fastq,XXX7.non_r.R2.fastq,XXX8.non_r.R2.fastq -label XXX1,XXX2,XXX3,XXX4,XXX5,XXX6,XXX7,XXX8 -o XXX`

#### 4. Call peaks
4.1 bpCon can works on any regions of the genome given by the users. Here, we also provied a simple script to get exons of the transcripts of the human or mouse genome.

>`perl ./Scripts/single_base_peaks_calling_P1.1.pl ./human/gencode.v39.annotation.gtf hg38_exon.bed hg38_gene.bed`

4.2 Call peaks among the given regions.
>`perl Call_Peaks.pl -cores 6 -m 24 -bed hg38_exon.bed -r normal -adj_p 0.05 -fc 0.58 -CPM 0 -o UC001 -f_IN XXX1.forward.sort.bam,XXX2.forward.sort.bam -f_IP XXX5.forward.sort.bam,XXX6.forward.sort.bam -r_IN XXX1.reverse.sort.bam,XXX2.reverse.sort.bam -r_IP XXX5.reverse.sort.bam,XXX6.reverse.sort.bam -c 15719959,15201522,24482557,26558549`

#### 5. Call peak difference
bpCon highly depended on the enrichment of IP relative to Input samples, so the  data for inputs and IPs need to correspond one-to-one.

>`perl Peak_Difference.pl -cores 6 -bed hg38_exon.bed -r normal -adj_p 0.001 -fc 0.58 -CPM 1 -o UiPSM_vs_UC001 -peak_C UC001.peak -peak_T UiPSM.peak -f1_IN XXX1.forward.sort.bam,XXX2.forward.sort.bam -f1_IP XXX5.forward.sort.bam,XXX6.forward.sort.bam -f2_IN XXX3.forward.sort.bam,XXX4.forward.sort.bam -f2_IP XXX7.forward.sort.bam,XXX8.forward.sort.bam -r1_IN XXX1.reverse.sort.bam,XXX2.reverse.sort.bam -r1_IP XXX5.reverse.sort.bam,XXX6.reverse.sort.bam -r2_IN XXX3.reverse.sort.bam,XXX4.reverse.sort.bam -r2_IP XXX7.reverse.sort.bam,XXX8.reverse.sort.bam -c 15719959,15201522,24482557,26558549,17379154,11392213,28406168,29231147`

### Maintainers and Contributing
boCon is developed and maintained by Yushuai Wang (wang_yushuai@gzlab.ac.cn).

### License
This pipeline is released under MIT license.
