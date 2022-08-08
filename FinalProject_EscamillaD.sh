#!/bin/bash
##To run the script you need to provide the folder with the raw data and reference genome folder
dir0=/home/i2gds/resources/DataFinalProject/MUS
##Setting out directories
dir1=$HOME/FinalProject_EscamillaD;
dir2=$dir1/00_rawdata;
dir3=$dir1/01_processed;
dir4=$dir1/02_mapping;
dir5=$dir1/03_quantification;
dir6=$dir5/01_cufflinks;
dir7=$dir5/02_cuffmerge;
dir8=$dir5/03_cuffquant;
dir9=$dir5/04_cuffdiff;
dir10=$dir7/Assembly;
dir11=$dir1/04_results

##Making directory FInalProject
cd $HOME
   mkdir FinalProject_EscamillaD

##Making subdirectories within LinuxProject for downstream steps
cd $dir1
   mkdir -p $dir1/{00_rawdata,01_processed,02_mapping,03_quantification/{01_cufflinks,02_cuffmerge/Assembly,03_cuffquant,04_cuffdiff},04_results}

##Copying fastq files from home/i2gds/resources/DataFinalProject/MUS/subset directory using a bash script
cd $dir0
for f in *.fastq; do
   cp -r $f $dir2/${f%.fastq}.fastq
done

##Copying Illumina adapater to the raw data file. Don't forget to replace the folder path for the adapter fasta file
cp /home/i2gds/resources/DataFinalProject/IlluminaAdapters_V2.fasta $dir2/IlluminaAdapters_V2.fasta

##Copying reference files to  my FinalProject/00_rawdata directory using a bash script
cd $dir0
   cp *.{fa.gz,gtf.gz} $dir2

##Unziping fasta.gz files
cd $dir2 
for f in *.fa.gz; do
   gunzip $dir2/${f%.fa.gz}.fa.gz
done

##unziping gtf.gz file
cd $dir2
for f in *.gtf.gz; do
   gunzip $dir2/${f%.gtf.gz}.gtf.gz
done

##processing  the reads with fastq-mcf 
cd $dir2
for f1 in *.fastq; do
f2=${f1%.fastq}
   echo "$f1"
   echo "$f2"
   fastq-mcf -q 30 -l 50 -o $dir3/${f2}.fastq $dir2/IlluminaAdapters_V2.fasta $dir2/$f1
done
   
##making subdirectories for each named file under the 02_mapping directory
cd $dir2
for f in *.fastq; do
   mkdir $dir4/"${f%.fastq}"
done

##Making subdirectories for each named file under the 03_cuffquant directory
cd $dir2
for f in *.fastq; do
   mkdir $dir8/"${f%.fastq}"
done

##Making text file assembly_list
cd $dir10
touch assembly_list.txt

##Fastq-Stats for rawdata and putting the resulst in the  corresponding stats.txt files in 00_rawdata
cd $dir2
for f1 in *.fastq; do
   f2=${f1%.fastq}
   echo "$f1"
   echo "$f2"
   fastq-stats $f1>$dir2/${f2}.stats.txt;
done

##Fastq-Stats for processed data and putting the resulst in the  corresponding stats.txt files in 01_processed
cd $dir3
for f1 in *.fastq; do
   f2=${f1%.fastq}
   echo "$f1"
   echo "$f2"
   fastq-stats $f1>$dir3/${f2}.stats.txt;
done

##Generation of genome index with STAR BE SURE TO SET THE .fa AND .GTF files
STAR --runThreadN 1 \
     --runMode genomeGenerate \
     --genomeDir $dir4 \
     --genomeFastaFiles $dir2/Mus_musculus.GRCm38.dna_rm.primary_assembly.fa \
     --sjdbGTFfile $dir2/Mus_musculus.GRCm38.86.gtf;

##Second mapping reads via STAR
cd $dir3
for f1 in *.fastq; do
f2=${f1%.fastq}
   echo "$f1"
   echo "$f2"
STAR --runThreadN 1 \
     --genomeDir $dir4 \
     --readFilesIn $dir3/$f1 \
     --outSAMstrandField intronMotif \
     --outFileNamePrefix $dir4/$f2/$f2 \
     --outSAMtype BAM SortedByCoordinate;
done

##Using samtools to get the number of mapped reads
cd $dir4
for f1 in $(find -name '*Aligned.sortedByCoord.out.bam'); do
   echo "$f1";
f3=${f1%Aligned.sortedByCoord.out.bam}.stats.txt
   touch $f3
   echo "$f3";

samtools view -F 4 $f1 | wc -l > $f3;
done

##Using Cufflinks for transcriptome assembly BE SURE TO SPECIFY GTF FILE
cd $dir2
bamsuffix=Aligned.sortedByCoord.out.bam;
for f1 in *.fastq; do
f2=${f1%.fastq}
echo $f2;
cufflinks --GTF-guide $dir2/Mus_musculus.GRCm38.86.gtf \
      	   --library-type fr-unstranded \
          --output-dir $dir6/$f2 \
          --num-threads 6 \
          $dir4/$f2/$f2$bamsuffix;
done

##Copying paths of transcript.gtf files into assembly_list.txt
cd $dir6
echo $dir6 | printf '%s\n' "$PWD"/*/transcripts.gtf > $dir10/assembly_list.txt

##Merging the transcript asseemblies INPUT GTF AND REFERENCE FILES
cuffmerge -o $dir7 -g $dir2/Mus_musculus.GRCm38.86.gtf -s $dir2/Mus_musculus.GRCm38.dna_rm.primary_assembly.fa $dir10/assembly_list.txt

##Quantify and normalize the number of reads using cuffquant INPUT REFERENCE GENOME
cd $dir2
bamsuffix=Aligned.sortedByCoord.out.bam
for f1 in *.fastq; do
	f2=${f1%.fastq}
	echo $f2;
	cuffquant -o $dir8/$f2 -b $dir2/Mus_musculus.GRCm38.dna_rm.primary_assembly.fa -p 8 -u -library-type fr-firststrand $dir7/merged.gtf $dir4/$f2/$f2$bamsuffix;
done

#Identification of differential gene expression with cuffdiff
cuffdiff -o $dir9 -L Brain,Heart,Liver -b $dir2/Mus_musculus.GRCm38.dna_rm.primary_assembly.fa -p 12 $dir7/merged.gtf $dir8/Brain_1/abundances.cxb,$dir8/Brain_2/abundances.cxb $dir8/Heart_1/abundances.cxb,$dir8/Heart_2/abundances.cxb $dir8/Liver_1/abundances.cxb,$dir8/Liver_2/abundances.cxb;

#Run the r script to analysi the data
Rscript /home/i2gds/dianae91/RFinalProject_Escamilla.R 

exit;
