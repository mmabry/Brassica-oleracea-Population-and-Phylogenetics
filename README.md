# Population-and-Phylogenetics-analysis-of-Brassica-oleracea
Analyses run for Mabry et al. (in prep) The Evolutionary History of Wild and Domesticated Brassica oleracea (Brassicaceae)

## 1. MAPPING TO THE REFERENCE

Reference : [ftp://ftp.ensemblgenomes.org/pub/plants/release-41/fasta/brassica_oleracea/dna/Brassica_oleracea.BOL.dna.toplevel.fa.gz]

#### STAR 2-pass Alignment

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

###### Shell script for first pass : RunSTAR.sh
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
###### To run them all for first pass use
```bash
for file in *_R1_001.fastq; do sbatch RunSTAR.sh $file; done
```

###### Set up second pass
concatenate all the SJ.out.tab files into one file, perform the filtering, and use only this filtered file for --sjdbFileChrStartEnd

```bash
#genomeDir=/path/to/B_oleracea_ref_Pass2

mkdir B_oleracea_ref_Pass2

sbatch -N 1 -n 12 -p BioCompute --mem=80000 -t 2-00:00:00 --wrap="STAR --runMode genomeGenerate --genomeDir /group/pireslab/mmabry_hpc/Brassica_oleracea/runDIR/B_oleracea_ref_Pass2 --genomeFastaFiles Brassica_oleracea.BOL.dna.toplevel.fa --sjdbFileChrStartEnd /group/pireslab/mmabry_hpc/Brassica_oleracea/runDIR/all_SJ.out.tab --sjdbOverhang 75 --limitSjdbInsertNsj 1023698 --runThreadN 12"
```
then make run directory
```bash
mkdir runDIR_pass2
```
###### Shell script for second pass : RunSTAR_Pass2.sh
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
###### To run them all for second pass use
```bash
for file in *_R1_001.fastq; do sbatch RunSTAR_Pass2.sh $file; done
```

## 2. ADD READ GROUPS, SORT, MARK DUPLICATES, AND CREATE INDEX
```bash
module load java/openjdk/java-1.8.0-openjdk
module load picard-tools/picard-2.7.1
#lewis version of Picard did not work, could not find jar file, so I had to install it myself
```
###### Shell script to create read groups: ReadGroups1.sh
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

###### To run them all 
```bash
for file in *_Aligned.out.sam; do sbatch ReadGroups1.sh $file; done
```
###### Shell Script to mark duplicates : MarkDups.sh
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
###### To run them all 
```bash
for file in *_rg_added_sorted.bam; do sbatch MarkDups.sh $file; done
```

