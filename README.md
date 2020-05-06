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

###### make shell script for first pass
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

