## Reference GTF/FASTA: 
Ensembl mus musculus release 101 (mm10)

## Samples:
prefix	fastq
ctrl_1_pol2s5p	ctrl_1_pol2s5p_R1.fastq
ctrl_1_input	ctrl_1_input_R1.fastq
ctrl_2_pol2s5p	ctrl_2_pol2s5p_R1.fastq
ctrl_2_input	ctrl_2_input_R1.fastq
mut_1_pol2s5p	mut_1_pol2s5p_R1.fastq
mut_1_input	mut_1_input_R1.fastq
mut_2_pol2s5p	mut_2_pol2s5p_R1.fastq
mut_2_input	mut_2_input_R1.fastq

## Trimming: Trimmomatic
java -jar /mnt/software/x86_64/packages/trimmomatic/0.39/trimmomatic-0.39.jar SE -threads 8 ./raw/${sample_prefix}_R1.fastq ./trim/${sample_prefix}_R1.fastq HEADCROP:0 LEADING:0 TRAILING:0 SLIDINGWINDOW:5:15 CROP:500 AVGQUAL:0 MINLEN:15

## Mapping: STAR
/mnt/software/x86_64/packages/star/2.7.9a/STAR --genomeDir /mnt/flatfiles/organisms/new_organism/mus_musculus/101/index_star --runThreadN 16 --readFilesIn ./trim/${sample_prefix}_R1.fastq --outReadsUnmapped Fastx --outFileNamePrefix ./star/${sample_prefix} --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --genomeLoad LoadAndKeep --outFilterMismatchNoverLmax 0.2 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outFilterMatchNmin 20 --alignEndsProtrude 10 ConcordantPair --alignMatesGapMax 2000 --limitOutSAMoneReadBytes 10000000 --outMultimapperOrder Random --sjdbOverhang 100 --alignEndsType EndToEnd --alignIntronMax 1 --alignSJDBoverhangMin 999 --outFilterMultimapNmax 1

## Deduplication: PICARD
java -Xmx4g -jar /mnt/software/x86_64/packages/picard/2.21.7/picard.jar MarkDuplicates I=./star/${sample_prefix}.bam O=./star/${sample_prefix}_nodup.bam REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT

## Generate BigWig normalized to mapped reads per sample: python package deeptools 3.5.1 / bamCoverage
/mnt/software/x86_64/packages/python/3.8.0-miniconda-4.9.2-buster/bin/bamCoverage -b ./star/${sample_prefix}_nodup_filter.bam -o ./star/${sample_prefix}_nodup_filter.bw -p 16 --binSize 25 --smoothLength 75 --normalizeUsing RPKM --outFileFormat bigwig

