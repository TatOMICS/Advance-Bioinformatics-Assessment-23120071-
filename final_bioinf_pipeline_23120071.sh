#!/bin/bash

# 1. Install required tools for the project

## 1.1. Install Anaconda
cd ~/
wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh
chmod +x ./Anaconda3-2022.10-Linux-x86_64.sh
bash ./Anaconda3-2022.10-Linux-x86_64.sh
source ~/.bashrc

## 1.2 Install required packages with Anaconda 

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install samtools
conda install bwa
conda install freebayes
conda install picard
conda install bedtools
conda install trimmomatic
conda install fastqc
conda install vcflib
conda install -c bioconda trimmomatic # I know it it's a duplicate. if I don't do it won't work. 

## 1.3. Install remaining needed software
sudo apt install bwa

## 1.4. Initial directory structure (I keep the other folders inside the assesment folder so it's easy to delete)
cd ~
mkdir ngs_assesment_FINAL
mkdir ngs_assesment_FINAL/dnaseq
cd ngs_assesment_FINAL/dnaseq
mkdir data meta results logs


## 1.5. Create a simple README description of the project
cd ~/ngs_assesment_FINAL/dnaseq/data
touch README.md
echo “KCL finally assessment for the module Advance Bioinformatics, Year 2024, Wish me luck.” > README.md

## Use 'cat README.md' to check what’s inside.


# 2. Download the files.

## 2.1 Create the untrimmed and trimmed files
cd ~/ngs_assesment_FINAL/dnaseq/data
mkdir untrimmed_fastq
mkdir trimmed_fastq

## 2.2 Download the files FastQ R1 and R2 and bed file
## Files are inside the 'data' fodler.
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

## 2.3 Move the raw fastq files to  the'untrimmed_fastq' directory.
## Leave the bed file inside the 'data' directory
mv NGS0001.R1.fastq.qz NGS0001.R2.fastq.qz ~/ngs_assesment_FINAL/dnaseq/data/untrimmed_fastq

## 2.4 Convert .qz FASTQ files to .gz files
cd ~/ngs_assesment_FINAL/dnaseq/data/untrimmed_fastq
mv NGS0001.R1.fastq.qz NGS0001.R1.fastq.gz
mv NGS0001.R2.fastq.qz NGS0001.R2.fastq.gz

## 2.5 Download reference genome FASTA file into a reference file (data > reference)
## Create the 'Reference' folder first.
mkdir -p ~/ngs_assesment_FINAL/dnaseq/data/Reference
cd ~/ngs_assesment_FINAL/dnaseq/data/Reference
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

## 2.6  Generate index of reference genome fasta file for alignment.
## This step takes ages 3959.920 seC.
bwa index hg19.fa.gz

## 2.7 Download & setup Annovar database for future variant anotation.
## https://www.openbioinformatics.org/annovar/annovar_download_form.php
## Once logged a link is provided.
## http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
cd ~/ngs_assesment_FINAL
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
## Unpack and set up the file.
tar -zxvf annovar.latest.tar.gz

## 2.8 Download Annovar db
cd annovar
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/

## If you get a "permission denied" error message when you try to download the DBs, this might be because of annotate_variation.pl not being set as an executable file.
##  Please fix this but running "chmod +x annotate_variation.pl" in your temrinal and try again

# 3 FastQC quality assessment on raw data
cd ~/ngs_assesment_FINAL/dnaseq/data/untrimmed_fastq
## Note the extra parameter we specified for 4 threads
fastqc -t 4 *.fastq.gz

## 3.1 Move FastQC results to a new directory
mkdir ~/ngs_assesment_FINAL/dnaseq/results/fastqc_untrimmed_reads
mv *fastqc* ~/ngs_assesment_FINAL/dnaseq/results/fastqc_untrimmed_reads

## 3.2 Download R 1 and R2 html using FileZilla !!!

# 4 Perform Trimmomatic (triming) on raw sequencing data
## Trimmomatic trims adapters and drops reads below 50 bp in length, removes Illumina adapters and trims end of reads bases if below a threshold quality (25).
cd ~/ngs_assesment_FINAL/dnaseq/data/untrimmed_fastq

trimmomatic PE \
-threads 4 \
-phred33 \
~/ngs_assesment_FINAL/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq.gz ~/ngs_assesment_FINAL/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq.gz \
-baseout ~/ngs_assesment_FINAL/dnaseq/data/trimmed_fastq/NGS0001_chr22m_trimmed_R \
ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
TRAILING:25 MINLEN:50

## 4.1. Untrimmed data is not longer needed. We will need space later.
cd ~/ngs_assesment_FINAL/dnaseq/data/untrimmed_fastq
rm NGS0001.R1.fastq.gz
rm NGS0001.R2.fastq.gz

# 5. Trimmed sequencing data quality assesment using FastQC.
cd ~/ngs_assesment_FINAL/dnaseq/data/trimmed_fastq
fastqc -t 4 NGS0001_chr22m_trimmed_R_1P NGS0001_chr22m_trimmed_R_2P

## 5.1  Move FastQC to the trimmed results directory
mkdir ~/ngs_assesment_FINAL/dnaseq/results/fastqc_trimmed_reads
mv *fastqc* ~/ngs_assesment_FINAL/dnaseq/results/fastqc_trimmed_reads

# 6 Alignment of reads to reference genome
mkdir ~/ngs_assesment_FINAL/dnaseq/data/aligned_data

# 6.1 alignment of trimmed reads to  hg19 reference genome, using BWA-MEM algorithm.
# 6.2 read group info obtained from the raw read FASTQ files.

bwa mem -t 4 -v 1 \
-R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.WES01\tSM:WES01\tPL:ILLUMINA\tLB:nextera-wes01-blood\tDT:2017-02-23\tPU:HWI-D00119' \
-I 250,50 \
~/ngs_assesment_FINAL/dnaseq/data/Reference/hg19.fa.gz \
~/ngs_assesment_FINAL/dnaseq/data/trimmed_fastq/NGS0001_chr22m_trimmed_R_1P \
~/ngs_assesment_FINAL/dnaseq/data/trimmed_fastq/NGS0001_chr22m_trimmed_R_2P \
> ~/ngs_assesment_FINAL/dnaseq/data/aligned_data/NGS0001.sam

# It’s taking ages !!!!

# The file is now created inside /home/ubuntu/ngs_assesment_FINAL/dnaseq/data/aligned_data/NGS0001.sam

# 7. Convert, process and index SAM file
cd ~/ngs_assesment_FINAL/dnaseq/data/aligned_data
# Sam to BAM
samtools view -h -b NGS0001.sam > NGS0001.bam
# sort the bam file
samtools sort NGS0001.bam > NGS0001_sorted.bam
# generate index
samtools index NGS0001_sorted.bam

# remove unneeded data to save disk space
rm NGS0001.sam

# 7.1 Mark duplicate reads
# this step locates and marks duplicated reads in the BAM file
cd ~/ngs_assesment_FINAL/dnaseq/data/aligned_data
picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt

# 7.2 Filter BAM based on mapping quality and bitwise flags using samtools  (Should I called marked_filtered)??????
samtools view -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam
samtools index NGS0001_sorted_filtered.bam

#.7.3 Generate alignment statistic
cd ~/ngs_assesment_FINAL/dnaseq/data/aligned_data

# 7.4 make a new directory for these stats
mkdir aligment_stats

# 7.5  Generate flagstats
samtools flagstat NGS0001_sorted_filtered.bam > flagstats_output.txt
mv flagstats_output.txt ~/ngs_assesment_FINAL/dnaseq/data/aligned_data/aligment_stats/

# 7.6 view the BAM file in the command line
cd ~/ngs_assesment_FINAL/dnaseq/data/aligned_data
samtools view -h NGS0001_sorted_filtered.bam | less # Perhaps i need to get rid of this for the automation

# 7.7  Generate alignment statistics per chromosome with idxstats
samtools idxstats NGS0001_sorted_filtered.bam > idxstats_output.txt
mv idxstats_output.txt ~/ngs_assesment_FINAL/dnaseq/data/aligned_data/aligment_stats/

# 7.8 Determine the distribution of insert sizes between read pairs with Picard tools
java -jar /home/ubuntu/anaconda3/pkgs/picard-2.18.29-0/share/picard-2.18.29-0/picard.jar CollectInsertSizeMetrics \
  I=/home/ubuntu/ngs_assesment_FINAL/dnaseq/data/aligned_data/NGS0001_sorted_filtered.bam \
  O=insert_size_metrics.txt \
  H=insert_size_histogram.pdf

## Download from FileZilla (inside aligened_data)

# 7.9 Calculate Depth of coverage.
# First get the all BAM regions that overlap with the input, then use this output to calculate the coverage.
bedtools intersect -bed -a NGS0001_sorted_filtered.bam -b $HOME/ngs_assesment_FINAL/dnaseq/data/annotation.bed \
| bedtools coverage -d -a $HOME/ngs_assesment_FINAL/dnaseq/data/annotation.bed -b - \
> coverageBed_output.txt

# 8. Variant calling with Freebayes
cd ~/ngs_assesment_FINAL/dnaseq/data
zcat ~/ngs_assesment_FINAL/dnaseq/data/Reference/hg19.fa.gz > ~/ngs_assesment_FINAL/dnaseq/data/Reference/hg19.fa
samtools faidx ~/ngs_assesment_FINAL/dnaseq/data/Reference/hg19.fa

# 8.1 FreeBayes to report putative variants in the sequencing data compared to the reference allelereebayes \
freebayes \
--bam ~/ngs_assesment_FINAL/dnaseq/data/aligned_data/NGS0001_sorted_filtered.bam \
--fasta-reference ~/ngs_assesment_FINAL/dnaseq/data/Reference/hg19.fa \
--vcf ~/ngs_assesment_FINAL/dnaseq/results/NGS0001_chr22m.vcf
# compress VCF
bgzip ~/ngs_assesment_FINAL/dnaseq/results/NGS0001_chr22m.vcf
# Index the vcf file
tabix -p vcf ~/ngs_assesment_FINAL/dnaseq/results/NGS0001_chr22m.vcf.gz

# remove unneeded data to save disk space
rm ~/ngs_assesment_FINAL/dnaseq/data/Reference/hg19.fa

# 9.Filtering the VCF
# remove "bad" variant calls from the VCF.
# "QUAL > 1" removes low quality data, "QUAL / AO > 10" additional contribution of each obs should be
# 	10 log units (~ Q10 per read), "SAF > 0 & SAR > 0" reads on both strands, "RPR > 1 & RPL > 1" at
# 	least two reads “balanced” to each side of the site.
vcffilter -f 'QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1' \
~/ngs_assesment_FINAL/dnaseq/results/NGS0001_chr22m.vcf.gz > \
~/ngs_assesment_FINAL/dnaseq/results/NGS0001_chr22m_filtered.vcf.gz

# 9.1 use bedtools to filter VCF for regions present in the annotation bed file
bedtools intersect -header -wa -a ~/ngs_assesment_FINAL/dnaseq/results/NGS0001_chr22m_filtered.vcf.gz \
  -b ~/ngs_assesment_FINAL/dnaseq/data/annotation.bed \
  > ~/ngs_assesment_FINAL/dnaseq/results/NGS0001_filtered_in_bedfile.vcf

# 9.2 compress filtered VCF (Not sure if required)
bgzip -f ~/ngs_assesment_FINAL/dnaseq/results/NGS0001_filtered_in_bedfile.vcf

# 9.3 Indexing of the vcf.
tabix -p vcf ~/ngs_assesment_FINAL/dnaseq/results/NGS0001_filtered_in_bedfile.vcf.gz

# 10. Variant Annotation with snpEff
## Install snpEff
mkdir ~/snpEff
cd ~/snpEff
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd snpEff
## Config snpEff
java -Xmx4g -jar snpEff.jar download hg19

# Annotate the variants using snpEff for functional effects prediction
java -Xmx4g -jar snpEff.jar hg19 \
    ~/ngs_assesment_FINAL/dnaseq/results/NGS0001_filtered_in_bedfile.vcf.gz \
    > ~/ngs_assesment_FINAL/dnaseq/results/NGS0001_filtered_in_bedfile_snpEff.vcf

# 10.1 Compress and index the snpEff annotated VCF
bgzip -f ~/ngs_assesment_FINAL/dnaseq/results/NGS0001_filtered_in_bedfile_snpEff.vcf
tabix -p vcf ~/ngs_assesment_FINAL/dnaseq/results/NGS0001_filtered_in_bedfile_snpEff.vcf.gz

# 10.2 Convert VCF to ANNOVAR input (if still needed after snpEff annotation)
~/ngs_assesment_FINAL/annovar/convert2annovar.pl \
-format vcf4 \
~/ngs_assesment_FINAL/dnaseq/results/NGS0001_filtered_in_bedfile_snpEff.vcf.gz \
> ~/ngs_assesment_FINAL/dnaseq/results/NGS0001_chr22m_filtered_chr22.avinput


# 10.3 Annovar for annotation of the variants with database frequencies/functional consequences inclusifn snp138
~/ngs_assesment_FINAL/annovar/table_annovar.pl ~/ngs_assesment_FINAL/dnaseq/results/NGS0001_chr22m_filtered_chr22.avinput \
/home/ubuntu/ngs_assesment_FINAL/annovar/humandb/ -buildver hg19 \
-out ~/ngs_assesment_FINAL/dnaseq/results/NGS0001_filtered_bedfile -remove \
-protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro,snp138 -operation g,g,f,f,f,f \
-otherinfo -nastring . -csvout

# 10.4 Use Filezilla to download the file.








