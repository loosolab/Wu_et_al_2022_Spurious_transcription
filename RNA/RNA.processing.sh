REMOVE THIS LATER:
\yua9_Fan_RNAseq_July2018\rnaseq_2_replicates

## Samples:
prefix	fastq
ctrl_1	ctrl_1_R1.fastq
ctrl_2	ctrl_2_R1.fastq
ko-tet3_1	ko-tet3_1_R1.fastq
ko-tet3_2	ko-tet3_2_R1.fastq

## Trimming
java -Djava.io.tmpdir=./trimmomatic_temp -jar /mnt/software/x86_64/packages/trimmomatic/0.39/trimmomatic-0.39.jar SE -threads 8 ./raw/${sample_prefix}_R1.fastq ./trim/${sample_prefix}_R1_trim.fastq HEADCROP:0 LEADING:0 TRAILING:0 SLIDINGWINDOW:5:15 CROP:500 AVGQUAL:0 MINLEN:15

## Mapping
STAR --genomeDir /mnt/flatfiles/organisms/new_organism/mus_musculus/101/index_star --runThreadN 16 --readFilesIn ./trim/${sample_prefix}_R1_trim.fastq --outReadsUnmapped Fastx --outFileNamePrefix ./star/${sample_prefix} --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --genomeLoad LoadAndKeep --outFilterMismatchNoverLmax 0.1 --outFilterScoreMinOverLread 0.9 --outFilterMatchNminOverLread 0.9 --outFilterMatchNmin 20 --alignIntronMax 200000 --alignMatesGapMax 2000 --alignEndsProtrude 10 ConcordantPair --outMultimapperOrder Random --limitOutSAMoneReadBytes 10000000 --sjdbOverhang 100 --outFilterMultimapNmax 1

## Filter mapped reads
filter: use_rrna = no, use_chrm = no...
## mark rrna read ids
# todo

## remove rrna ids
java -Xmx24g -jar /mnt/software/x86_64/packages/picard/2.21.7/picard.jar FilterSamReads I=./star/${sample_prefix}.bam O=./star/${sample_prefix}.filter.bam READ_LIST_FILE=./star/${sample_prefix}_remove_ids.txt FILTER=excludeReadList


## RNA differential analysis
