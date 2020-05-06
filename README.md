# Population-and-Phylogenetics-analysis-of-Brassica-oleracea
Analyses run for Mabry et al. (in prep) The Evolutionary History of Wild and Domesticated Brassica oleracea (Brassicaceae)

## 1. MAPPING TO THE REFERENCE - STAR 2-PASS ALIGNMENT
All of the following GATK scripts from https://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq

Reference : [ftp://ftp.ensemblgenomes.org/pub/plants/release-41/fasta/brassica_oleracea/dna/Brassica_oleracea.BOL.dna.toplevel.fa.gz]

```bash
module load java/openjdk/java-1.8.0-openjdk
module load gatk/gatk-4.0.1.1
module load star/star-2.5.2b
```

```bash
mkdir B_oleracea_ref
sbatch -N 1 -n 12 -p BioCompute --mem=80000 -t 2-00:00:00 --wrap="STAR --runMode genomeGenerate --genomeDir /group/pireslab/mmabry_hpc/Brassica_oleracea/B_oleracea_ref --genomeFastaFiles Brassica_oleracea.BOL.dna.toplevel.fa  --runThreadN 12"
```

```bash
mkdir runDIR
cd runDIR
```

#### A. Shell script for first pass : RunSTAR.sh
```bash
#! /bin/bash

#SBATCH -p BioCompute,hpc5  # partition
#SBATCH -J STAR_reads  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)
#SBATCH --mem=80000 #memory

module load java/openjdk/java-1.8.0-openjdk
module load star/star-2.5.2b

PREFIX=$(echo $1 | awk -F_ '{print $1"_"$2"_"}')


STAR --genomeDir /group/pireslab/mmabry_hpc/Brassica_oleracea/B_oleracea_ref --readFilesIn ${PREFIX}R1_001.fastq ${PREFIX}R2_001.fastq --runThreadN 12 --outFileNamePrefix /group/pireslab/mmabry_hpc/Brassica_oleracea/runDIR/${PREFIX}
```
##### A.1 To run them all for first pass use
```bash
for file in *_R1_001.fastq; do sbatch RunSTAR.sh $file; done
```

#### B. Set up second pass and make shell script : RunSTAR_Pass2.sh
* concatenate all the SJ.out.tab files into one file, perform the filtering, and use only this filtered file for --sjdbFileChrStartEnd

```bash
#genomeDir=/path/to/B_oleracea_ref_Pass2

mkdir B_oleracea_ref_Pass2

sbatch -N 1 -n 12 -p BioCompute --mem=80000 -t 2-00:00:00 --wrap="STAR --runMode genomeGenerate --genomeDir /group/pireslab/mmabry_hpc/Brassica_oleracea/runDIR/B_oleracea_ref_Pass2 --genomeFastaFiles Brassica_oleracea.BOL.dna.toplevel.fa --sjdbFileChrStartEnd /group/pireslab/mmabry_hpc/Brassica_oleracea/runDIR/all_SJ.out.tab --sjdbOverhang 75 --limitSjdbInsertNsj 1023698 --runThreadN 12"
```
```bash
mkdir runDIR_pass2
```
```bash
#! /bin/bash

#SBATCH -p BioCompute,hpc5  # partition
#SBATCH -J STAR_reads2  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)
#SBATCH --mem=80000 #memory

module load java/openjdk/java-1.8.0-openjdk
module load star/star-2.5.2b

PREFIX=$(echo $1 | awk -F_ '{print $1"_"$2"_"}')


STAR --genomeDir /group/pireslab/mmabry_hpc/Brassica_oleracea/B_oleracea_ref_Pass2 --readFilesIn ${PREFIX}R1_001.fastq ${PREFIX}R2_001.fastq --runThreadN 12 --outFileNamePrefix /group/pireslab/mmabry_hpc/Brassica_oleracea/runDIR_pass2/${PREFIX}
```
##### B.1 To run them all for second pass use
```bash
for file in *_R1_001.fastq; do sbatch RunSTAR_Pass2.sh $file; done
```

## 2. ADD READ GROUPS, SORT, MARK DUPLICATES, AND CREATE INDEX
```bash
module load java/openjdk/java-1.8.0-openjdk
module load picard-tools/picard-2.7.1
#lewis version of Picard did not work, could not find jar file, so I had to install it myself
```
#### A. Shell script to create read groups: ReadGroups1.sh
```bash
#! /bin/bash

#SBATCH -p BioCompute,hpc5,Lewis  # partition
#SBATCH -J ReadGroups  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)
#SBATCH --mem=80000 #memory 

module load java/openjdk/java-1.8.0-openjdk

PREFIX=$(echo $1 | awk -F_ '{print $1"_"$2"_"}')

java -jar /home/mmabry/software/picard/picard.jar AddOrReplaceReadGroups I=${PREFIX}Aligned.out.sam O=${PREFIX}rg_added_sorted.bam SO=coordinate RGID=${PREFIX} RGLB=library1 RGPL=illumina RGPU=AHJKFJBGX5 RGSM=${PREFIX}

#RGID = ${PREFIX} # Read Group ID Default value: 1. This option can be set to 'null' to clear the default value.
#RGLB = ???????? # Read Group library Required. Lib1
#RGPL = illumina # Read Group platform (e.g. illumina, solid) Required.
#RGPU = ???????? # Read Group platform unit (eg. run barcode) Required: lane1:AHJKFJBGX5
#RGSM = ${PREFIX} # Read Group sample name Required.
```

##### A.1 To run them all 
```bash
for file in *_Aligned.out.sam; do sbatch ReadGroups1.sh $file; done
```
#### B. Shell Script to mark duplicates : MarkDups.sh
```bash
#! /bin/bash

#SBATCH -p BioCompute,hpc5,Lewis  # partition
#SBATCH -J MarkDups  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)
#SBATCH --mem=80000 #memory

module load java/openjdk/java-1.8.0-openjdk

PREFIX=$(echo $1 | awk -F_ '{print $1"_"$2"_"}')

java -jar /home/mmabry/software/picard/picard.jar MarkDuplicates I=${PREFIX}rg_added_sorted.bam O=${PREFIX}dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${PREFIX}output.metrics
```
##### B.1 To run them all 
```bash
for file in *_rg_added_sorted.bam; do sbatch MarkDups.sh $file; done
```

## 3. SPLIT'N'TRIM AND REASSIGN MAPPING QUALITIES
```bash
module load java/openjdk/java-1.8.0-openjdk
module load gatk/gatk-3.8
#CANNOT USE GATK4 for RNAseq yet
```
#### A. Make another index for GATK to work with
```bash
sbatch -N 1 -n 12 -p BioCompute --mem=80000 -t 2-00:00:00 --wrap="java -jar /home/mmabry/software/picard/picard.jar CreateSequenceDictionary R= /group/pireslab/mmabry_hpc/Brassica_oleracea/Brassica_oleracea.BOL.dna.toplevel.fa O= /group/pireslab/mmabry_hpc/Brassica_oleracea/Brassica_oleracea.BOL.dna.toplevel.fa.dict"

module load samtools/samtools-1.7
sbatch -N 1 -n 12 -p BioCompute --mem=80000 -t 2-00:00:00 --wrap="samtools faidx /group/pireslab/mmabry_hpc/Brassica_oleracea/Brassica_oleracea.BOL.dna.toplevel.fa"
```
#### B. Shell script to run split and trim : SPLITNTRIM.sh
```bash
#! /bin/bash

#SBATCH -p BioCompute,hpc5,Lewis  # partition
#SBATCH -J SplitNTrim  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)
#SBATCH --mem=80000 #memory

module load java/openjdk/java-1.8.0-openjdk
module load gatk/gatk-3.8

PREFIX=$(echo $1 | awk -F_ '{print $1"_"$2"_"}')


java -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar -T SplitNCigarReads -R /group/pireslab/mmabry_hpc/Brassica_oleracea/Brassica_oleracea.BOL.dna.toplevel.fa -I ${PREFIX}dedupped.bam -o ${PREFIX}split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
```
##### B.1 to run them all use
```bash
for file in *_dedupped.bam; do sbatch SplitNTrim.sh $file; done
```

## 4. VARIANT CALLING USING GVCF

#### A. Shell script : VarCall_GVCF.sh
```bash
#! /bin/bash

#SBATCH -p BioCompute,hpc5,Lewis  # partition
#SBATCH -J VarCall  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)
#SBATCH --mem=80000 #memory

module load java/openjdk/java-1.8.0-openjdk
module load gatk/gatk-3.8

PREFIX=$(echo $1 | awk -F_ '{print $1"_"$2"_"}')

java -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar -T HaplotypeCaller -R /group/pireslab/mmabry_hpc/Brassica_oleracea/Brassica_oleracea.BOL.dna.toplevel.fa -I ${PREFIX}split.bam --emitRefConfidence GVCF -o ${PREFIX}output.raw.snps.indels.g.vcf.gz
```
##### A.1 To run them all
```bash
for file in *_split.bam; do sbatch VarCall_GVCF.sh $file; done
```

#### B. Combine samples to make a GVCF file using Shell script : GVCF_all.sh
```bash
#! /bin/bash

#SBATCH -p BioCompute,hpc5,Lewis  # partition
#SBATCH -J GVCF_all  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)
#SBATCH --mem=80000 #memory

module load java/openjdk/java-1.8.0-openjdk
module load gatk/gatk-3.8

 java -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar \
    -T GenotypeGVCFs \
    -R /group/pireslab/mmabry_hpc/Brassica_oleracea/Brassica_oleracea.BOL.dna.toplevel.fa \
    --variant lane1/10_S10_output.raw.snps.indels.g.vcf \
    --variant lane1/11_S11_output.raw.snps.indels.g.vcf \
    --variant lane1/12_S12_output.raw.snps.indels.g.vcf \
    --variant lane1/13_S13_output.raw.snps.indels.g.vcf \
    --variant lane1/14_S14_output.raw.snps.indels.g.vcf \
    --variant lane1/15_S15_output.raw.snps.indels.g.vcf \
    --variant lane1/17_S16_output.raw.snps.indels.g.vcf \
    --variant lane1/18_S17_output.raw.snps.indels.g.vcf \
    --variant lane1/19_S18_output.raw.snps.indels.g.vcf \
    --variant lane1/1_S1_output.raw.snps.indels.g.vcf \
    --variant lane1/20_S19_output.raw.snps.indels.g.vcf \
    --variant lane1/21_S20_output.raw.snps.indels.g.vcf \
    --variant lane1/22_S21_output.raw.snps.indels.g.vcf \
    --variant lane1/23_S22_output.raw.snps.indels.g.vcf \
    --variant lane1/24_S23_output.raw.snps.indels.g.vcf \
    --variant lane1/25_S24_output.raw.snps.indels.g.vcf \
    --variant lane1/2_S2_output.raw.snps.indels.g.vcf \
    --variant lane1/3_S3_output.raw.snps.indels.g.vcf \
    --variant lane1/4_S4_output.raw.snps.indels.g.vcf \
    --variant lane1/5_S5_output.raw.snps.indels.g.vcf \
    --variant lane1/6_S6_output.raw.snps.indels.g.vcf \
    --variant lane1/7_S7_output.raw.snps.indels.g.vcf \
    --variant lane1/8_S8_output.raw.snps.indels.g.vcf \
    --variant lane1/9_S9_output.raw.snps.indels.g.vcf \
    --variant lane2/26_S1_output.raw.snps.indels.g.vcf \
    --variant lane2/27_S2_output.raw.snps.indels.g.vcf \
    --variant lane2/28_S3_output.raw.snps.indels.g.vcf \
    --variant lane2/31_S4_output.raw.snps.indels.g.vcf \
    --variant lane2/34_S5_output.raw.snps.indels.g.vcf \
    --variant lane2/35_S6_output.raw.snps.indels.g.vcf \
    --variant lane2/37_S7_output.raw.snps.indels.g.vcf \
    --variant lane2/38_S8_output.raw.snps.indels.g.vcf \
    --variant lane2/40_S9_output.raw.snps.indels.g.vcf \
    --variant lane2/41_S10_output.raw.snps.indels.g.vcf \
    --variant lane2/42_S11_output.raw.snps.indels.g.vcf \
    --variant lane2/43_S12_output.raw.snps.indels.g.vcf \
    --variant lane2/44_S13_output.raw.snps.indels.g.vcf \
    --variant lane2/45_S14_output.raw.snps.indels.g.vcf \
    --variant lane2/46_S15_output.raw.snps.indels.g.vcf \
    --variant lane2/47_S16_output.raw.snps.indels.g.vcf \
    --variant lane2/48_S17_output.raw.snps.indels.g.vcf \
    --variant lane2/49_S18_output.raw.snps.indels.g.vcf \
    --variant lane2/50_S19_output.raw.snps.indels.g.vcf \
    --variant lane2/51_S20_output.raw.snps.indels.g.vcf \
    --variant lane2/52_S21_output.raw.snps.indels.g.vcf \
    --variant lane2/53_S22_output.raw.snps.indels.g.vcf \
    --variant lane2/54_S23_output.raw.snps.indels.g.vcf \
    --variant lane2/55_S24_output.raw.snps.indels.g.vcf \
    --variant lane3/56_S1_output.raw.snps.indels.g.vcf \
    --variant lane3/57_S2_output.raw.snps.indels.g.vcf \
    --variant lane3/58_S3_output.raw.snps.indels.g.vcf \
    --variant lane3/59_S4_output.raw.snps.indels.g.vcf \
    --variant lane3/61_S5_output.raw.snps.indels.g.vcf \
    --variant lane3/62_S6_output.raw.snps.indels.g.vcf \
    --variant lane3/63_S7_output.raw.snps.indels.g.vcf \
    --variant lane3/64_S8_output.raw.snps.indels.g.vcf \
    --variant lane3/65_S9_output.raw.snps.indels.g.vcf \
    --variant lane3/66_S10_output.raw.snps.indels.g.vcf \
    --variant lane3/67_S11_output.raw.snps.indels.g.vcf \
    --variant lane3/68_S12_output.raw.snps.indels.g.vcf \
    --variant lane3/69_S13_output.raw.snps.indels.g.vcf \
    --variant lane3/70_S14_output.raw.snps.indels.g.vcf \
    --variant lane3/71_S15_output.raw.snps.indels.g.vcf \
    --variant lane3/72_S16_output.raw.snps.indels.g.vcf \
    --variant lane3/73_S17_output.raw.snps.indels.g.vcf \
    --variant lane3/74_S18_output.raw.snps.indels.g.vcf \
    --variant lane3/75_S19_output.raw.snps.indels.g.vcf \
    --variant lane3/76_S20_output.raw.snps.indels.g.vcf \
    --variant lane3/77_S21_output.raw.snps.indels.g.vcf \
    --variant lane3/78_S22_output.raw.snps.indels.g.vcf \
    --variant lane3/79_S23_output.raw.snps.indels.g.vcf \
    --variant lane3/80_S24_output.raw.snps.indels.g.vcf \
    --variant lane4/100_S19_output.raw.snps.indels.g.vcf \
    --variant lane4/101_S20_output.raw.snps.indels.g.vcf \
    --variant lane4/102_S21_output.raw.snps.indels.g.vcf \
    --variant lane4/103_S22_output.raw.snps.indels.g.vcf \
    --variant lane4/104_S23_output.raw.snps.indels.g.vcf \
    --variant lane4/106_S24_output.raw.snps.indels.g.vcf \
    --variant lane4/81_S1_output.raw.snps.indels.g.vcf \
    --variant lane4/82_S2_output.raw.snps.indels.g.vcf \
    --variant lane4/84_S3_output.raw.snps.indels.g.vcf \
    --variant lane4/85_S4_output.raw.snps.indels.g.vcf \
    --variant lane4/86_S5_output.raw.snps.indels.g.vcf \
    --variant lane4/87_S6_output.raw.snps.indels.g.vcf \
    --variant lane4/88_S7_output.raw.snps.indels.g.vcf \
    --variant lane4/89_S8_output.raw.snps.indels.g.vcf \
    --variant lane4/90_S9_output.raw.snps.indels.g.vcf \
    --variant lane4/91_S10_output.raw.snps.indels.g.vcf \
    --variant lane4/92_S11_output.raw.snps.indels.g.vcf \
    --variant lane4/93_S12_output.raw.snps.indels.g.vcf \
    --variant lane4/94_S13_output.raw.snps.indels.g.vcf \
    --variant lane4/95_S14_output.raw.snps.indels.g.vcf \
    --variant lane4/96_S15_output.raw.snps.indels.g.vcf \
    --variant lane4/97_S16_output.raw.snps.indels.g.vcf \
    --variant lane4/98_S17_output.raw.snps.indels.g.vcf \
    --variant lane4/99_S18_output.raw.snps.indels.g.vcf \
    --variant lane5/107_S1_output.raw.snps.indels.g.vcf \
    --variant lane5/108_S2_output.raw.snps.indels.g.vcf \
    --variant lane5/109_S3_output.raw.snps.indels.g.vcf \
    --variant lane5/110_S4_output.raw.snps.indels.g.vcf \
    --variant lane5/111_S5_output.raw.snps.indels.g.vcf \
    --variant lane5/113_S6_output.raw.snps.indels.g.vcf \
    --variant lane5/114_S7_output.raw.snps.indels.g.vcf \
    --variant lane5/115_S8_output.raw.snps.indels.g.vcf \
    --variant lane5/116_S9_output.raw.snps.indels.g.vcf \
    --variant lane5/117_S10_output.raw.snps.indels.g.vcf \
    --variant lane5/118_S11_output.raw.snps.indels.g.vcf \
    --variant lane5/120_S12_output.raw.snps.indels.g.vcf \
    --variant lane5/121_S13_output.raw.snps.indels.g.vcf \
    --variant lane5/122_S14_output.raw.snps.indels.g.vcf \
    --variant lane5/123_S15_output.raw.snps.indels.g.vcf \
    --variant lane5/124_S16_output.raw.snps.indels.g.vcf \
    --variant lane5/125_S17_output.raw.snps.indels.g.vcf \
    --variant lane5/126_S18_output.raw.snps.indels.g.vcf \
    --variant lane5/128_S19_output.raw.snps.indels.g.vcf \
    --variant lane5/129_S20_output.raw.snps.indels.g.vcf \
    --variant lane5/130_S21_output.raw.snps.indels.g.vcf \
    --variant lane5/131_S22_output.raw.snps.indels.g.vcf \
    --variant lane5/132_S23_output.raw.snps.indels.g.vcf \
    --variant lane5/133_S24_output.raw.snps.indels.g.vcf \
    --variant lane6/134_S1_output.raw.snps.indels.g.vcf \
    --variant lane6/135_S2_output.raw.snps.indels.g.vcf \
    --variant lane6/136_S3_output.raw.snps.indels.g.vcf \
    --variant lane6/137_S4_output.raw.snps.indels.g.vcf \
    --variant lane6/139_S5_output.raw.snps.indels.g.vcf \
    --variant lane6/140_S6_output.raw.snps.indels.g.vcf \
    --variant lane6/141_S7_output.raw.snps.indels.g.vcf \
    --variant lane6/142_S8_output.raw.snps.indels.g.vcf \
    --variant lane6/143_S9_output.raw.snps.indels.g.vcf \
    --variant lane6/144_S10_output.raw.snps.indels.g.vcf \
    --variant lane6/145_S11_output.raw.snps.indels.g.vcf \
    --variant lane6/147_S12_output.raw.snps.indels.g.vcf \
    --variant lane6/148_S13_output.raw.snps.indels.g.vcf \
    --variant lane6/150_S14_output.raw.snps.indels.g.vcf \
    --variant lane6/151_S15_output.raw.snps.indels.g.vcf \
    --variant lane6/152_S16_output.raw.snps.indels.g.vcf \
    --variant lane6/155_S17_output.raw.snps.indels.g.vcf \
    --variant lane6/157_S18_output.raw.snps.indels.g.vcf \
    --variant lane6/159_S19_output.raw.snps.indels.g.vcf \
    --variant lane6/160_S20_output.raw.snps.indels.g.vcf \
    --variant lane6/161_S21_output.raw.snps.indels.g.vcf \
    --variant lane6/162_S22_output.raw.snps.indels.g.vcf \
    --variant lane6/163_S23_output.raw.snps.indels.g.vcf \
    --variant lane6/167_S24_output.raw.snps.indels.g.vcf \
    --variant lane7/169_S1_output.raw.snps.indels.g.vcf \
    --variant lane7/170_S2_output.raw.snps.indels.g.vcf \
    --variant lane7/171_S3_output.raw.snps.indels.g.vcf \
    --variant lane7/172_S4_output.raw.snps.indels.g.vcf \
    --variant lane7/173_S5_output.raw.snps.indels.g.vcf \
    --variant lane7/174_S6_output.raw.snps.indels.g.vcf \
    --variant lane7/175_S7_output.raw.snps.indels.g.vcf \
    --variant lane7/176_S8_output.raw.snps.indels.g.vcf \
    --variant lane7/177_S9_output.raw.snps.indels.g.vcf \
    --variant lane7/178_S10_output.raw.snps.indels.g.vcf \
    --variant lane7/179_S11_output.raw.snps.indels.g.vcf \
    --variant lane7/180_S12_output.raw.snps.indels.g.vcf \
    --variant lane7/181_S13_output.raw.snps.indels.g.vcf \
    --variant lane7/182_S14_output.raw.snps.indels.g.vcf \
    --variant lane7/183_S15_output.raw.snps.indels.g.vcf \
    --variant lane7/184_S16_output.raw.snps.indels.g.vcf \
    --variant lane7/185_S17_output.raw.snps.indels.g.vcf \
    --variant lane7/186_S18_output.raw.snps.indels.g.vcf \
    --variant lane7/187_S19_output.raw.snps.indels.g.vcf \
    --variant lane7/188_S20_output.raw.snps.indels.g.vcf \
    --variant lane7/189_S21_output.raw.snps.indels.g.vcf \
    --variant lane7/190_S22_output.raw.snps.indels.g.vcf \
    --variant lane7/191_S23_output.raw.snps.indels.g.vcf \
    --variant lane7/192_S24_output.raw.snps.indels.g.vcf \
    --variant lane8/119_S7_output.raw.snps.indels.g.vcf \
    --variant lane8/127_S8_output.raw.snps.indels.g.vcf \
    --variant lane8/138_S9_output.raw.snps.indels.g.vcf \
    --variant lane8/153_S10_output.raw.snps.indels.g.vcf \
    --variant lane8/154_S11_output.raw.snps.indels.g.vcf \
    --variant lane8/168_S12_output.raw.snps.indels.g.vcf \
    --variant lane8/193_S13_output.raw.snps.indels.g.vcf \
    --variant lane8/194_S14_output.raw.snps.indels.g.vcf \
    --variant lane8/199_S15_output.raw.snps.indels.g.vcf \
    --variant lane8/203_S16_output.raw.snps.indels.g.vcf \
    --variant lane8/204_S17_output.raw.snps.indels.g.vcf \
    --variant lane8/207_S18_output.raw.snps.indels.g.vcf \
    --variant lane8/208_S19_output.raw.snps.indels.g.vcf \
    --variant lane8/212_S20_output.raw.snps.indels.g.vcf \
    --variant lane8/213_S21_output.raw.snps.indels.g.vcf \
    --variant lane8/214_S22_output.raw.snps.indels.g.vcf \
    --variant lane8/220_S23_output.raw.snps.indels.g.vcf \
    --variant lane8/222_S24_output.raw.snps.indels.g.vcf \
    --variant lane8/30_S1_output.raw.snps.indels.g.vcf \
    --variant lane8/32_S2_output.raw.snps.indels.g.vcf \
    --variant lane8/36_S3_output.raw.snps.indels.g.vcf \
    --variant lane8/39_S4_output.raw.snps.indels.g.vcf \
    --variant lane8/60_S5_output.raw.snps.indels.g.vcf \
    --variant lane8/83_S6_output.raw.snps.indels.g.vcf \
    --variant lane9/105_S3_output.raw.snps.indels.g.vcf \
    --variant lane9/146_S4_output.raw.snps.indels.g.vcf \
    --variant lane9/149_S5_output.raw.snps.indels.g.vcf \
    --variant lane9/156_S6_output.raw.snps.indels.g.vcf \
    --variant lane9/164_S7_output.raw.snps.indels.g.vcf \
    --variant lane9/16_S1_output.raw.snps.indels.g.vcf \
    --variant lane9/216_S8_output.raw.snps.indels.g.vcf \
    --variant lane9/221_S9_output.raw.snps.indels.g.vcf \
    --variant lane9/226_S10_output.raw.snps.indels.g.vcf \
    --variant lane9/227_S11_output.raw.snps.indels.g.vcf \
    --variant lane9/230_S12_output.raw.snps.indels.g.vcf \
    --variant lane9/231_S13_output.raw.snps.indels.g.vcf \
    --variant lane9/232_S14_output.raw.snps.indels.g.vcf \
    --variant lane9/233_S15_output.raw.snps.indels.g.vcf \
    --variant lane9/234_S16_output.raw.snps.indels.g.vcf \
    --variant lane9/33_S2_output.raw.snps.indels.g.vcf \
    --variant lane10/112_S2_output.raw.snps.indels.g.vcf \
    --variant lane10/165_S3_output.raw.snps.indels.g.vcf \
    --variant lane10/195_S4_output.raw.snps.indels.g.vcf \
    --variant lane10/196_S5_output.raw.snps.indels.g.vcf \
    --variant lane10/198_S6_output.raw.snps.indels.g.vcf \
    --variant lane10/205_S7_output.raw.snps.indels.g.vcf \
    --variant lane10/209_S8_output.raw.snps.indels.g.vcf \
    --variant lane10/215_S9_output.raw.snps.indels.g.vcf \
    --variant lane10/218_S10_output.raw.snps.indels.g.vcf \
    --variant lane10/219_S11_output.raw.snps.indels.g.vcf \
    --variant lane10/224_S12_output.raw.snps.indels.g.vcf \
    --variant lane10/228_S13_output.raw.snps.indels.g.vcf \
    --variant lane10/229_S14_output.raw.snps.indels.g.vcf \
    --variant lane10/236_S15_output.raw.snps.indels.g.vcf \
    --variant lane10/237_S16_output.raw.snps.indels.g.vcf \
    --variant lane10/238_S17_output.raw.snps.indels.g.vcf \
    --variant lane10/29_S1_output.raw.snps.indels.g.vcf \
    -o GVCF_all_output.vcf
```
#### C. Make a data table using Shell script: VarientsTable.sh
```bash
#! /bin/bash

#SBATCH -J VarientsTB
#SBATCH -o VarientsTB.o%J
#SBATCH -e VarientsTB.e%J
#SBATCH --partition CLUSTER
#SBATCH --ntasks 1
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 20G 
 

java -jar /share/apps/GATK/current/GenomeAnalysisTK.jar \
     -R /home/mem2c2/data/Brassica_oleracea.v2.1.fa \
     -T VariantsToTable \
     -V GVCF_all_output.vcf \
     -F CHROM -F POS -F ID -F QUAL -F AC -F QD -F FS \
     -o Bo_results2.table
```
* to plot these use R 
* to check SNP number
```bash
grep -v "^#" GVCF_all_output.vcf | wc -l
```

## 5. FILTERING VARIANTS (VCF)
```bash
#! /bin/bash

#SBATCH -J FilterVCF
#SBATCH -o FilterVCF.o%J
#SBATCH -e FilterVCF.e%J
#SBATCH --partition CLUSTER
#SBATCH --ntasks 1
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 20G

module load java

java -jar /share/apps/GATK/current/GenomeAnalysisTK.jar -T VariantFiltration -R /home/mem2c2/data/Brassica_oleracea.v2.1.fa -V /home/mem2c2/data/GVCF_all_output.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o /home/mem2c2/data/FilteredGVCF_allBo.vcf
```
```bash
! /bin/bash

#SBATCH -J SelectVarients
#SBATCH -o SelectVarients.o%J
#SBATCH -e SelectVarients.e%J
#SBATCH --partition CLUSTER
#SBATCH --ntasks 1
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 20G

module load java

 java -jar /share/apps/GATK/current/GenomeAnalysisTK.jar \
   -R /home/mem2c2/data/Brassica_oleracea.v2.1.fa \
   -T SelectVariants \
   -V FilteredGVCF_allBo.vcf \
   -o GVCF_selectVar.vcf \
   -env \
   -ef
```

## 6. REMOVE ALL SCAFFOLDS, ADD SNP ID, AND FILTER FOR MISSING DATA
```bash
module load samtools/samtools-1.8
bgzip GVCF_selectVar.vcf
tabix GVCF_selectVar.vcf.gz
```
```bash
##then create regions file (Bo_noScaff.txt) 
grep -v '^#' GVCF_selectVar.vcf | awk '{print $1"\t"$2}'
#then delete scaffolds from the bottom
```
#### A. Remove Scaffolds
```bash
#! /bin/bash

#SBATCH -J RemoveScaffs
#SBATCH -o RemoveScaffs.o%J
#SBATCH -e RemoveScaffs.e%J
#SBATCH --partition CLUSTER
#SBATCH --ntasks 1
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 20G

module load bcftools-1.8

bcftools view GVCF_selectVar.vcf.gz -R Bo_noScaff.txt -o GVCF_selectVar_noScaff.vcf
```
#### B. Add ID
```bash
#! /bin/bash

#SBATCH -J AddID
#SBATCH -o AddID.o%J
#SBATCH -e AddID.e%J
#SBATCH --partition CLUSTER
#SBATCH --ntasks 1
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 20G

module load bcftools-1.7 

bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' GVCF_selectVar_noScaff.vcf -o Bo_noScaff_ID_FiltMiss.vcf
```

#### C. Filter Missing, on depth, and indels
```bash
#! /bin/bash

#SBATCH -p BioCompute,hpc5  # partition
#SBATCH -J FilterVCF  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)
#SBATCH --mem=80G #memory

module load vcftools/vcftools-v0.1.14

vcftools --recode --recode-INFO-all --vcf GVCF_selectVar_noScaff.vcf --max-missing 0.4 --min-meanDP 5 --remove-indels --remove-indv 188_S20_ --out GVCF_selectVar_noScaff_FiltMiss_depth_snps_no188.vcf
```

#### D. Check number of SNPs
```bash
grep -v '^#' Bo_noScaff_ID_FiltMiss.vcf.recode.vcf | wc -l
```

## 7. PLINK CONVERSION and LD filtering 
```bash
#! /bin/bash

#SBATCH -p BioCompute,hpc5  # partition
#SBATCH -J VCF_PLINK_LD  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)
#SBATCH --mem=80G #memory


module load plink/plink-1.90b 

plink --vcf GVCF_selectVar_noScaff_FiltMiss_depth_snps_no188.vcf --make-bed --double-id --allow-extra-chr

plink --bfile plink --indep 40kb 5 2 --allow-extra-chr

plink --bfile plink --extract plink.prune.in --make-bed --out pruned50 --allow-extra-chr

plink --bfile pruned50 --recode vcf-fid --out GVCF_selectVar_noScaff_FiltMiss_depth_snps_no188_LD40kb.vcf --allow-extra-chr
```

## 8. PCAngsd [http://www.popgen.dk/software/index.php/PCAngsd#Quick_start]
#### A. sort and index using samtools before running, use bam files from mapping to genome in step 1. Shell Script : SortIndex.sh
```bash
#! /bin/bash

#SBATCH -p BioCompute,hpc5,Lewis  # partition
#SBATCH -J Sort  # give the job a custom name
#SBATCH -o Sort-%j.out  # give the job output a custom name
#SBATCH -t 01:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)
#SBATCH --mem=80G #memory

module load samtools/samtools-1.9 

PREFIX=$(echo $1 | awk -F. '{print $1}')

samtools sort ${PREFIX}.bam -o ${PREFIX}.sorted.bam

samtools index ${PREFIX}.sorted.bam
```
##### A.1 to run them all in a loop
```bash
for file in *_split.bam; do sbatch SortIndex.sh $file; done
```
```bash
git clone https://github.com/Rosemeis/pcangsd.git
cd pcangsd/
pip install --user -r python_packages.txt
```
#### B. Prepare input files using BAM and ANGSD
```bash
ls *sorted.bam > bam.filelist
```
#### C. Run Angsd
```bash
#! /bin/bash

#SBATCH -p BioCompute  # partition
#SBATCH -J Angsd  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 54  # number of cores (AKA tasks)
#SBATCH --mem=250G #memory

module load htslib/htslib-1.9


/home/mmabry/software/angsd-0.925/angsd -GL 1 -out B_oleracea.angsd -nThreads 14 -doGlf 2 -doMajorMinor 1 -doMaf 2 -minMapQ 30 -SNP_pval 1e-6 -bam bam.filelist

#/home/mmabry/software/angsd-0.925/angsd -GL 1 -out B_oleracea_1_9.angsd -nThreads 14 -doGlf 2 -doMajorMinor 1 -doMaf 2 -minMapQ 30 -rf regions.txt -SNP_pval 1e-6 -bam bam.filelist # this will run the regions listed in the regions.txt file, for example only chromosomes 1-9

#/home/mmabry/software/angsd-0.925/angsd -GL 1 -out B_oleracea.angsd -nThreads 14 -doGlf 2 -doMajorMinor 1 -doMaf 2 -SNP_pval 1e-6 -r C1: -bam bam.filelist   #just Chromosome 1
```
* #to run just chromosomes 1-9 make separate region list file and use -rf in command above
```bash
C1:
C2:
C3:
C4:
C5:
C6:
C7:
C8:
C9:
```
#### D. Run PCangsd
```bash
#! /bin/bash

#SBATCH -p BioCompute,hpc5  # partition
#SBATCH -J PCAngsd  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)
#SBATCH --mem=80G #memory


export PYTHONPATH=/home/mmabry/.local/lib/python2.7/site-packages/:$PYTHONPATH
module load py-numpy/py-numpy-1.14.2-openblas-python-2.7.14

#python /home/mmabry/software/pcangsd-0.97/pcangsd.py -beagle B_oleracea_allC_noCol_title_angsd.beagle.gz -admix -inbreed 2 -o B_oleracea_PCangsd -threads 14

python /home/mmabry/software/pcangsd-0.97/pcangsd.py -beagle B_oleracea_allC_noCol_title_angsd.beagle.gz -e 23 -admix -selection 2 -inbreed 2 -o B_oleracea_PCangsd_Select_Admix -threads 14 -sites_save
```


## 9. SNPhylo [http://chibba.pgml.uga.edu/snphylo/]
#### A. Remove C's in front of Chromosome number
```bash
sed -r 's/C([1-9])/\1/g' Bo_noScaff_ID_FiltMiss.vcf.recode.vcf > Bo_noScaff_ID_FiltMiss.noC.vcf
```
#### B. Shell Script for SNPhylo 
```bash
#! /bin/bash

#SBATCH -p BioCompute,hpc5,Lewis  # partition
#SBATCH -J SNPhylo  # give the job a custom name
#SBATCH -o SNPhylo-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)
#SBATCH --mem=80G #memory

module load snphylo/snphylo-2016-02-04-python-2.7.14-tk


snphylo.sh -v Bo_noScaff_ID_FiltMiss_LD50prunded.snps.no188.noC.rename.vcf -a 9 -l 0.1 -m 0.01 -M 0.4 -b -B 1000 -o 238_S17
```

## 10. FastSTRUCTURE [https://rajanil.github.io/fastStructure/]
#### A. Shell script : FastStucture.sh
```bash
#! /bin/bash

#SBATCH -p BioCompute,hpc5,Lewis  # partition
#SBATCH --account=biosci
#SBATCH -J FastStucture  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)
#SBATCH --mem=80000 #memory

module load miniconda3/miniconda3-4.3.30 
source activate fastStructure-condaenv

python /storage/hpc/data/mmabry/src/fastStructure/structure.py -K 4 --input=plink --output=wild_domest --full --seed=100
```
#### B. Choose K Shell script : ChooseK.sh
```bash
#! /bin/bash

#SBATCH -p BioCompute,hpc5,Lewis  # partition
#SBATCH --account=biosci
#SBATCH -J chooseK  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)
#SBATCH --mem=80000 #memory

module load miniconda3/miniconda3-4.3.30 
source activate fastStructure-condaenv

python /storage/hpc/data/mmabry/src/fastStructure/chooseK.py --input=wild_domest
```

## 11. IQ-tree POMO [http://www.iqtree.org/doc/Polymorphism-Aware-Models]
#### A. Use only the pure samples (as determined using FastStructure) in a vcf
```bash
module load vcftools/vcftools-0.1.17

vcftools --recode --vcf Bo_noScaff_ID_FiltMiss_LD50prunded.snps.vcf --keep KeepPureSample_all.txt --out Bo_noScaff_ID_FiltMiss_LD50prunded.snps.vcf
```
#### B. Get file format correct
```bash
#! /bin/bash

#SBATCH -J PDGSpider
#SBATCH -o PDGSpider.o%J
#SBATCH -e PDGSpider.e%J
#SBATCH --partition CLUSTER
#SBATCH --ntasks 1
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 20G

java -Xmx2000m -Xms512m -jar /home/mem2c2/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile Bo_noScaff_ID_FiltMiss_LD50prunded.snps.vcf -inputformat VCF -outputfile Bo_noScaff_ID_FiltMiss.vcf.recode.fasta -outputformat FASTA 
```
###### example template
```bash
# spid-file generated: Wed May 29 21:15:03 CDT 2019

# VCF Parser questions
PARSER_FORMAT=VCF

# Only output SNPs with a phred-scaled quality of at least:
VCF_PARSER_QUAL_QUESTION=
# Select population definition file:
VCF_PARSER_POP_FILE_QUESTION=
# What is the ploidy of the data?
VCF_PARSER_PLOIDY_QUESTION=
# Do you want to include a file with population definitions?
VCF_PARSER_POP_QUESTION=
# Output genotypes as missing if the phred-scale genotype quality is below:
VCF_PARSER_GTQUAL_QUESTION=
# Do you want to include non-polymorphic SNPs?
VCF_PARSER_MONOMORPHIC_QUESTION=TRUE
# Only output following individuals (ind1, ind2, ind4, ...):
VCF_PARSER_IND_QUESTION=
# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
VCF_PARSER_REGION_QUESTION=
# Output genotypes as missing if the read depth of a position for the sample is below:
VCF_PARSER_READ_QUESTION=
# Take most likely genotype if "PL" or "GL" is given in the genotype field?
VCF_PARSER_PL_QUESTION=
# Do you want to exclude loci with only missing data?
VCF_PARSER_EXC_MISSING_LOCI_QUESTION=

# FASTA Writer questions
WRITER_FORMAT=FASTA

# Numeric SNP data: enter the integer that codes for nucleotide A:
FASTA_WRITER_CODE_A_QUESTION=
# Specify the locus you want to write to the FASTA file:
FASTA_WRITER_LOCUS_COMBINATION_QUESTION=
# Numeric SNP data: enter the integer that codes for nucleotide G:
FASTA_WRITER_CODE_G_QUESTION=
# Numeric SNP data: enter the integer that codes for nucleotide C:
FASTA_WRITER_CODE_C_QUESTION=
# Do you want to save sequences on a single line?
FASTA_SINGLE_LINE_QUESTION=TRUE
# Numeric SNP data: enter the integer that codes for nucleotide T:
FASTA_WRITER_CODE_T_QUESTION=
# Specify the DNA locus you want to write to the FASTA file, write "CONCAT" for concatenation or "SEPARATE" to separate the loci:
FASTA_WRITER_CONCATENATE_QUESTION=
# Specify which data type should be included in the FASTA file  (FASTA can only analyze one data type per file):
FASTA_WRITER_DATA_TYPE_QUESTION=SNP
# Do you want to save haploid sequences (consensus sequence with ambiguity codes is taken if ploidy is higher)?
FASTA_WRITER_HAPLOID_QUESTION=TRUE
```
#### C. Rename population header to match the format below
```bash
>costata-33_S2_ | population:costata
```
#### D. Remove and remaining ? and replace with N
```bash
sed 's/?/N/g' Bo_noScaff_ID_FiltMiss_LD50prunded.snps.fasta > Bo_noScaff_ID_FiltMiss_LD50prunded.snps2.fasta 
```
#### E. Run Fasta to Counts script
```bash
#! /bin/bash

#SBATCH -J Fasta2counts
#SBATCH -o Fasta2counts.o%J
#SBATCH -e Fasta2counts.e%J
#SBATCH --partition CLUSTER
#SBATCH --ntasks 1
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 20G

module load Python/Python-3.6.4 

/home/mem2c2/cflib/scripts/FastaToCounts.py Bo_noScaff_ID_FiltMiss.recode.fasta.gz Bo_noScaff_ID_FiltMiss.recode.cf.gz
```
#### F. Run POMO
```bash
#! /bin/bash

#SBATCH -J POMO_BS
#SBATCH -o POMO_BS.o%J
#SBATCH -e POMO_BS.e%J
#SBATCH --partition CLUSTER
#SBATCH --ntasks 1
#SBATCH --nodes 1
#SBATCH --cpus-per-task 14
#SBATCH --mem 40G
#SBATCH --time 2-00:00:00

module load iqtree/iqtree-1.6.10

iqtree -s Bo_noScaff_ID_FiltMiss_LD50prunded.snps2.cf -m GTR+P -bb 1000 -o rupestris -pre BOOT
```

## 12. TreeMix [https://speciationgenomics.github.io/Treemix/]
#### A. First make .clust file by hand using a text editing program such as BBedit  (from .fam file), see example below
```bash
58_S3_ 58_S3_ sabellica
59_S4_ 59_S4_ sabellica
5_S5_ 5_S5_ palmifolia
```
#### B. Then use the plink2treemix.py script from the authors of treemix to prepare plink data for treemix
```bash
#! /bin/bash

#SBATCH -p BioCompute  # partition
#SBATCH -J plink2treemix  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)
#SBATCH --mem=80000 #memory 

##make sure plink2treemix.py is in the same directory as when you are running this

python plink2treemix.py plink.frq.strat.gz treemix.frq.gz
```
#### C. Run treemix
```bash
#! /bin/bash

#SBATCH -p BioCompute  # partition
#SBATCH -J treemix  # give the job a custom name
#SBATCH -o results-%j.out  # give the job output a custom name
#SBATCH -t 2-00:00:00  # two day time limit

#SBATCH -N 1  # number of nodes
#SBATCH -n 14  # number of cores (AKA tasks)
#SBATCH --mem=80000 #memory 


for i in {2..9}

do
    /home/mem2c2/treemix/src/treemix -i treemix.frq.gz -root montana -k 300 -m $i -o Bo_300_${i}_Montana

done
```
#### D. Plot results in R
```R
source("/Users/mem2c2/OneDrive - University of Missouri/Computer/Projects/BoleraceaPopGen/TreeMix/plotting_funcs.R")
library(ggplot2)

#read in file
setwd("/Users/mem2c2/OneDrive - University of Missouri/Computer/Projects/BoleraceaPopGen/TreeMix/70pure/")

pdf("70PureSamples_IndWild.pdf", width = 25, height = 40)
par(mfrow=c(3,3))

m.2 = plot_tree("Bo_treemix_70Pure_2_indivWild", cex = 1)
m.3 = plot_tree("Bo_treemix_70Pure_3_indivWild", cex = 1)
m.4 = plot_tree("Bo_treemix_70Pure_4_indivWild", cex = 1)
m.5 = plot_tree("Bo_treemix_70Pure_5_indivWild", cex = 1)
m.6 = plot_tree("Bo_treemix_70Pure_6_indivWild", cex = 1)
m.7 = plot_tree("Bo_treemix_70Pure_7_indivWild", cex = 1)
m.8 = plot_tree("Bo_treemix_70Pure_8_indivWild", cex = 1)
m.9 = plot_tree("Bo_treemix_70Pure_9_indivWild", cex = 1)

dev.off()

pdf("PlotResid.pdf", width = 8, height = 8)
par(mar=c(5,5,1,1)+0.1)
r.2 = plot_resid("Bo_treemix_70Pure_2_indivWild", pop_order = "popList.txt")
r.3 = plot_resid("Bo_treemix_70Pure_3_indivWild", pop_order = "popList.txt")
r.4 = plot_resid("Bo_treemix_70Pure_4_indivWild", pop_order = "popList.txt")
r.5 = plot_resid("Bo_treemix_70Pure_5_indivWild", pop_order = "popList.txt", cex = 0.8)
r.6 = plot_resid("Bo_treemix_70Pure_6_indivWild", pop_order = "popList.txt")
r.7 = plot_resid("Bo_treemix_70Pure_7_indivWild", pop_order = "popList.txt")
r.8 = plot_resid("Bo_treemix_70Pure_8_indivWild", pop_order = "popList.txt")
r.9 = plot_resid("Bo_treemix_70Pure_9_indivWild", pop_order = "popList.txt")

dev.off()

##check how much varience is explaned 
var_2 <- get_f("Bo_treemix_70Pure_2_indivWild")
var_3 <- get_f("Bo_treemix_70Pure_3_indivWild")
var_4 <- get_f("Bo_treemix_70Pure_4_indivWild")
var_5 <- get_f("Bo_treemix_70Pure_5_indivWild")
var_6 <- get_f("Bo_treemix_70Pure_6_indivWild")
var_7 <- get_f("Bo_treemix_70Pure_7_indivWild")
var_8 <- get_f("Bo_treemix_70Pure_8_indivWild")
var_9 <- get_f("Bo_treemix_70Pure_9_indivWild")


VarienceExplained <- c(var_2, var_3, var_4, var_5, var_6, var_7, var_8, var_9)
MigrationEdges <- c("2", "3", "4", "5", "6", "7", "8", "9")

varTable <- data.frame(MigrationEdges, VarienceExplained)

ggplot(varTable, aes(x= MigrationEdges, y = VarienceExplained)) +
  theme(text = element_text(size = 30)) +
  geom_text(aes(label = migrationEdges, size = 30)) +
  theme(axis.title.x=element_text(size = 15), axis.text.x=element_blank()) +
  theme(axis.title.y=element_text(size = 15), axis.text.y=element_text(size = 15)) +
  ylim(0.9,1) +
  theme(legend.position = "none")
```

## 13. F-statistics

## 14. Salmon [https://combine-lab.github.io/salmon/getting_started/]

## 15. tximport and DEseq2

## 16. WGCNA [https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/]



