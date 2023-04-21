#!/bin/bash

#Tools
TRIMMOMATIC=/genomics/software/Trimmomatic-0.39/trimmomatic-0.39.jar
CUTADAPT=/home/xiao/.local/bin/cutadapt ##v4.1
FASTX_COLLAPSER=/genomics/software/fastx_toolkit-0.0.13/fastx_collapser
BOWTIE2=/genomics/software/bowtie2-2.4.2-linux-x86_64/bowtie2
SAMTOOLS=/genomics/software/samtools-1.15.1/samtools

#Reference
ADAPTERS=/SOFTWARE/bbmap-38.96/resources/adapters.fa
NCRNA=/genomics/reference/small-rna/Rfam/rfam_rRNA-tRNA-snRNA-snoRNA ## Rfam: The RNA families database (Release 14.8)
REPEAT=/genomics/reference/small-rna/msRepDB/Hordeum_vulgare_repeats ## msRepDB (https://msrepdb.cbrc.kaust.edu.sa/pages/msRepDB/index.html), accessed on 20/07/2022
MIRNA=/genomics/reference/small-rna/miRBase/mature-miRNA ## The miRBase Sequence Database -- Release 22.1

INPUT="fastqs"
OUTPUT="aligned-miRNA"
mkdir -p $OUTPUT

THREADS=60
ulimit -n 4096

# Sample command to align all Fastq in a dir: (ECHO ONLY)
for f in $INPUT/*R1_*.fastq.gz; do
  R1=`basename $f`
  sample=${R1%.fastq.gz}

  echo "$INPUT/$R1"
  #echo "$OUTPUT/$sample.clean.fastq.gz"
  
  #removal of low-quality reads and adapter sequences
  ##Note: "adapters.fa" is from bbmap: /SOFTWARE/bbmap-38.96/resources/adapters.fa 
  echo "## step1: adapter and low-quality removal is running.."
  java -jar $TRIMMOMATIC SE -threads $THREADS $INPUT/$R1 $OUTPUT/$sample.clean.fastq.gz ILLUMINACLIP:$ADAPTERS:2:30:10 SLIDINGWINDOW:10:20

  #poly-A tails and 4 bases trimmed, too short/long removed
  echo "## step2: poly-A tails removal is runnning.."
  $CUTADAPT --cores=$THREADS -u 4 -a A{8} -m 18 -M 34 -o $OUTPUT/$sample.clean.polyAtrimmed.fastq.gz $OUTPUT/$sample.clean.fastq.gz

  #fastq converted into fasta format and collapsed into unique tags
  echo "## step3: collapsing into unique tags is runnning.."
  gzip -dc $OUTPUT/$sample.clean.polyAtrimmed.fastq.gz | paste - - - - | sed 's/^@/>/g' | cut -f 1,2 | tr '\t' '\n' > $OUTPUT/$sample.clean.polyAtrimmed.fasta
  #$FASTX_COLLAPSER -i Sample-1-25ng_S5_L002_R1_001.clean.polyAtrimmed.fasta -o unique_tags.fa
  #sed -i 's/-/_x/' $OUTPUT/$sample.clean.polyAtrimmed.fasta

  #rRNA, tRNA, snRNA, and snoRNA removal
  echo "## step4: ncRNA removal is running.."
  $BOWTIE2 --threads $THREADS -x $NCRNA -f $OUTPUT/$sample.clean.polyAtrimmed.fasta -S $OUTPUT/$sample.ncRNA.sam --un $OUTPUT/$sample.ncRNA.filtered.fa

  #repeats removal
  echo "## step5: repeat removal is running.."
  $BOWTIE2 --threads $THREADS -x $REPEAT -f $OUTPUT/$sample.ncRNA.filtered.fa -S $OUTPUT/$sample.ncRNA.repeats.sam --un $OUTPUT/$sample.ncRNA.repeats.filtered.fa

  #aligned against the known plant miRNAs from uhe miRBase database
  echo "## step6: aligning to the mirbase miRNA is running.."
  $BOWTIE2 --threads $THREADS -x $MIRNA -N 1 -f $OUTPUT/$sample.ncRNA.repeats.filtered.fa -S $OUTPUT/$sample.miRNA.sam --un $OUTPUT/$sample.known-miRNA.unaligned.fa
  $SAMTOOLS view -F 4 $OUTPUT/$sample.miRNA.sam > $OUTPUT/$sample.miRNA.alinged.sam
  grep -v "^@" $OUTPUT/$sample.miRNA.alinged.sam | awk -F"\t" '{if($3!="*") print}' |cut -f 1,3 | sort -u > $OUTPUT/$sample.known_miRNAs.txt

  echo "## all steps done.."
done
