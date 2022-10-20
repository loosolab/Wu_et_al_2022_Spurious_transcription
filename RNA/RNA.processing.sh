REMOVE THIS LATER:
\yua9_Fan_RNAseq_July2018\rnaseq_2_replicates

## Samples:
prefix	fastq
ctrl_1	ctrl_1_R1.fastq
ctrl_2	ctrl_2_R1.fastq
ko-tet3_1	ko-tet3_1_R1.fastq
ko-tet3_2	ko-tet3_2_R1.fastq

## Trimming
java -Djava.io.tmpdir=./trimmomatic_temp -jar /mnt/software/x86_64/packages/trimmomatic/0.39/trimmomatic-0.39.jar SE -threads 8 ./raw/${sample_prefix}_R1.fastq trimmomatic/${sample_prefix}_R1_trim.fastq HEADCROP:0 LEADING:0 TRAILING:0 SLIDINGWINDOW:5:15 CROP:500 AVGQUAL:0 MINLEN:15

## RNA preprocessing

## RNA mammping and counting

## RNA differential analysis
