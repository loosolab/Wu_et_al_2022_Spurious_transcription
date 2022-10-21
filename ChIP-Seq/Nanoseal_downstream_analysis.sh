

#downstream anaylsis on Nanoseal-seq

#merge the bw files accroding to genotype
bigwigCompare --numberOfProcessors 11 \
              --outFileName TETKO_in_libnorm.bw \
              --outFileFormat bigwig \
              --bigwig1 KO_1_spikein_fragnorm_frag.bw \
              --bigwig2 KO_3_spikein_fragnorm_frag.bw \
              --operation mean
              
              
bigwigCompare --numberOfProcessors 11 \
              --outFileName Ctrl_in_libnorm.bw \
              --outFileFormat bigwig \
              --bigwig1 WT_2_spikein_fragnorm_frag.bw \
              --bigwig2 WT_3_spikein_fragnorm_frag.bw \
              --operation mean


#quantile seperation by gene expression (RPKM in RNA-seq)
#RNA_seq_RPKM_Q1.csv
#RNA_seq_RPKM_Q2.csv
#RNA_seq_RPKM_Q3.csv
#RNA_seq_RPKM_Q4.csv

#annotation the gene region in mm10
join -1 1 -2 4 -o 2.1,2.2,2.3,2.4,2.5,2.6 \
<(sort -k 1 RNA_seq_RPKM_Q1.csv) \
<(sort -k 4 Mus_musculus.GRCm38.100_genes.bed) \
| sed "s/ /\t/g" > RNA_seq_RPKM_Q1.bed

join -1 1 -2 4 -o 2.1,2.2,2.3,2.4,2.5,2.6 \
<(sort -k 1 RNA_seq_RPKM_Q2.csv) \
<(sort -k 4 Mus_musculus.GRCm38.100_genes.bed) \
| sed "s/ /\t/g" > RNA_seq_RPKM_Q2.bed

join -1 1 -2 4 -o 2.1,2.2,2.3,2.4,2.5,2.6 \
<(sort -k 1 RNA_seq_RPKM_Q3.csv) \
<(sort -k 4 Mus_musculus.GRCm38.100_genes.bed) \
| sed "s/ /\t/g" > RNA_seq_RPKM_Q3.bed

join -1 1 -2 4 -o 2.1,2.2,2.3,2.4,2.5,2.6 \
<(sort -k 1 RNA_seq_RPKM_bottom5to25.csv) \
<(sort -k 4 Mus_musculus.GRCm38.100_genes.bed.csv) \
| sed "s/ /\t/g" > RNA_seq_RPKM_bottom5to25.bed

join -1 1 -2 4 -o 2.1,2.2,2.3,2.4,2.5,2.6 \
<(sort -k 1 RNA_seq_RPKM_bottom5.csv) \
<(sort -k 4 Mus_musculus.GRCm38.100_genes.bed.csv) \
| sed "s/ /\t/g" > RNA_seq_RPKM_bottom5.bed



#remove black list region
bedtools intersect -v -a RNA_seq_RPKM_Q1.bed \
                             -b mm10.blacklist.bed \
                             > RNA_seq_RPKM_Q1_rmbl.bed
bedtools intersect -v -a RNA_seq_RPKM_Q2.bed \
                             -b mm10.blacklist.bed \
                             > RNA_seq_RPKM_Q2_rmbl.bed
bedtools intersect -v -a RNA_seq_RPKM_Q3.bed \
                             -b mm10.blacklist.bed \
                             > RNA_seq_RPKM_Q3_rmbl.bed                            
bedtools intersect -v -a RNA_seq_RPKM_bottom5to25.bed \
                             -b mm10.blacklist.bed \
                             > RNA_seq_RPKM_bottom5to25_rmbl.bed
bedtools intersect -v -a RNA_seq_RPKM_bottom5.bed \
                             -b mm10.blacklist.bed \
                             > RNA_seq_RPKM_bottom5_rmbl.bed


computeMatrix scale-regions \
-S Ctrl_in_libnorm.bw \
TETKO_in_libnorm.bw \
-R RNA_seq_RPKM_Q1_rmbl.bed \
RNA_seq_RPKM_Q2_rmbl.bed \
RNA_seq_RPKM_Q3_rmbl.bed \
RNA_seq_RPKM_bottom5to25_rmbl.bed \
--outFileName Nanoseal_byRNAquantile_10kb_RmBl.mat.gz \
--regionBodyLength 10000 \
--beforeRegionStartLength 10000 \
--afterRegionStartLength 10000 \
--missingDataAsZero \
-p 11 \
-bl mm10.blacklist.bed


plotProfile -m Nanoseal_byRNAquantile_10kb_RmBl.mat.gz \
              --plotHeight 15 \
              --plotWidth 15 \
              --yMin 0 \
              --yMax 20 \
              -out Nanoseal_byRNAquantile_10kb_RmBl.pdf \
              --outFileNameData Nanoseal_byRNAquantile_10kb_RmBl \
              --regionsLabel "Top25%" "25%-50%" "50%-75%" "75%-95%" \
              --samplesLabel "Ctrl." "TET3KO" \
              --perGroup
              

computeMatrix scale-regions \
-S Ctrl_in_libnorm.bw \
TETKO_in_libnorm.bw \
-R RNA_seq_RPKM_bottom5_rmbl.bed \
--outFileName Nanoseal_bottom5percent_10kb_RmBl.mat.gz \
--regionBodyLength 10000 \
--beforeRegionStartLength 10000 \
--afterRegionStartLength 10000 \
--missingDataAsZero \
-p 11 \
-bl mm10.blacklist.bed


plotProfile -m Nanoseal_bottom5percent_10kb_RmBl.mat.gz \
              --plotHeight 15 \
              --plotWidth 15 \
              --yMin 0 \
              --yMax 20 \
              -out Nanoseal_byRNAquantile_10kb_RmBl.pdf \
              --outFileNameData Nanoseal_bottom5percent_10kb_RmBl \
              --regionsLabel "Bottom5%" \
              --samplesLabel "ctrl." "TET3KO" \
              --perGroup


#Count 5hmC signal at the genebody, TSS, and TES region


#annotate TSS region +- 200bp

awk 'OFS="\t" { print $1, $2-200, $2+200, $4,$5,$6,$7,$8}' RNA_seq_RPKM_Q1.bed >  RNA_seq_RPKM_Q1.TSS200.bed
awk 'OFS="\t" { print $1, $2-200, $2+200, $4,$5,$6,$7,$8}' RNA_seq_RPKM_Q2.bed >  RNA_seq_RPKM_Q2.TSS200.bed
awk 'OFS="\t" { print $1, $2-200, $2+200, $4,$5,$6,$7,$8}' RNA_seq_RPKM_Q3.bed >  RNA_seq_RPKM_Q3.TSS200.bed
awk 'OFS="\t" { print $1, $2-200, $2+200, $4,$5,$6,$7,$8}' RNA_seq_RPKM_Q4.bed >  RNA_seq_RPKM_Q4.TSS200.bed

#annotate TES region +-200bp

awk 'OFS="\t" { print $1, $3-200, $3+200, $4,$5,$6,$7,$8}' RNA_seq_RPKM_Q1.bed >  RNA_seq_RPKM_Q1.TES200.bed
awk 'OFS="\t" { print $1, $3-200, $3+200, $4,$5,$6,$7,$8}' RNA_seq_RPKM_Q2.bed >  RNA_seq_RPKM_Q2.TES200.bed
awk 'OFS="\t" { print $1, $3-200, $3+200, $4,$5,$6,$7,$8}' RNA_seq_RPKM_Q3.bed >  RNA_seq_RPKM_Q3.TES200.bed
awk 'OFS="\t" { print $1, $3-200, $3+200, $4,$5,$6,$7,$8}' RNA_seq_RPKM_Q4.bed >  RNA_seq_RPKM_Q4.TES200.bed


sed -i 's/\t*$//' RNA_seq_RPKM_Q1.TSS200.bed
sed -i 's/\t*$//' RNA_seq_RPKM_Q2.TSS200.bed
sed -i 's/\t*$//' RNA_seq_RPKM_Q3.TSS200.bed
sed -i 's/\t*$//' RNA_seq_RPKM_bottom5to25.TSS200.bed
sed -i 's/\t*$//' RNA_seq_RPKM_bottom5.TSS200.bed

sed -i 's/\t*$//' RNA_seq_RPKM_Q1.TES200.bed
sed -i 's/\t*$//' RNA_seq_RPKM_Q2.TES200.bed
sed -i 's/\t*$//' RNA_seq_RPKM_Q3.TES200.bed
sed -i 's/\t*$//' RNA_seq_RPKM_bottom5to25.TES200.bed
sed -i 's/\t*$//' RNA_seq_RPKM_bottom5.TES200.bed


#remove blacklist region
bedtools intersect -v -a RNA_seq_RPKM_Q1.TSS200.bed \
                             -b mm10.blacklist.bed \
                             > RNA-seq_Q1_mm10_filtered.TSS200_rmbl.bed
bedtools intersect -v -a RNA_seq_RPKM_Q2.TSS200.bed \
                             -b mm10.blacklist.bed \
                             > RNA-seq_Q2_mm10_filtered.TSS200_rmbl.bed
bedtools intersect -v -a RNA_seq_RPKM_Q3.TSS200.bed \
                             -b mm10.blacklist.bed \
                             > RNA-seq_Q3_mm10_filtered.TSS200_rmbl.bed
bedtools intersect -v -a RNA_seq_RPKM_bottom5to25.TSS200.bed \
                             -b mm10.blacklist.bed \
                             > RNAbottom5to25percent.TSS200_rmbl.bed
bedtools intersect -v -a RNA_seq_RPKM_bottom5.TSS200.bed \
                             -b mm10.blacklist.bed \
                             > RNAbottom5percent.TSS200_rmbl.bed              
        
bedtools intersect -v -a RNA_seq_RPKM_Q1.TES200.bed \
                             -b mm10.blacklist.bed \
                             > RNA-seq_Q1_mm10_filtered.TES200_rmbl.bed
bedtools intersect -v -a RNA_seq_RPKM_Q2.TES200.bed \
                             -b mm10.blacklist.bed \
                             > RNA-seq_Q2_mm10_filtered.TES200_rmbl.bed
bedtools intersect -v -a RNA_seq_RPKM_Q3.TES200.bed \
                             -b mm10.blacklist.bed \
                             > RNA-seq_Q3_mm10_filtered.TES200_rmbl.bed
bedtools intersect -v -a RNA_seq_RPKM_bottom5to25.TES200.bed \
                             -b mm10.blacklist.bed \
                             > RNAbottom5to25percent.TES200_rmbl.bed
bedtools intersect -v -a RNA_seq_RPKM_bottom5.TES200.bed \
                             -b mm10.blacklist.bed \
                             > RNAbottom5percent.TES200_rmbl.bed

#convert bed files to saf files
find . -name '*_rmbl.bed' | parallel awk \'OFS=\"\\t\"\ {print \$1\".\"\$2\".\"\$3\,\$1\,\$2\,\$3\,\".\"}\'\ \{}\ \>\{.}.saf


#count the reads in genomic region based on bam reads

#-----------------------genebody----------------------------
featureCounts -T 11 \
-F SAF \
-p \
-a RNA_seq_RPKM_Q1_rmbl.saf \
-o Nanoseal_featureCounts_gene_body_Q1.txt \
Ctrl_1.bam Ctrl_2.bam TETKO_1.bam TETKO_2.bam

featureCounts -T 11 \
-F SAF \
-p \
-a RNA_seq_RPKM_Q2_rmbl.saf \
-o Nanoseal_featureCounts_gene_body_Q2.txt \
Ctrl_1.bam Ctrl_2.bam TETKO_1.bam TETKO_2.bam

featureCounts -T 11 \
-F SAF \
-p \
-a RNA_seq_RPKM_Q3_rmbl.saf \
-o Nanoseal_featureCounts_gene_body_Q3.txt \
Ctrl_1.bam Ctrl_2.bam TETKO_1.bam TETKO_2.bam

featureCounts -T 11 \
-F SAF \
-p \
-a RNA_seq_RPKM_bottom5to25_rmbl.saf \
-o Nanoseal_featureCounts_gene_body_bottom5to25percent.txt \
Ctrl_1.bam Ctrl_2.bam TETKO_1.bam TETKO_2.bam

featureCounts -T 11 \
-F SAF \
-p \
-a RNA_seq_RPKM_bottom5_rmbl.saf \
-o Nanoseal_featureCounts_gene_body_bottom5percent.txt \
Ctrl_1.bam Ctrl_2.bam TETKO_1.bam TETKO_2.bam



#-------------------TSS +-200bp-------------------

featureCounts -T 11 \
-F SAF \
-p \
-a RNA-seq_Q1_mm10_filtered.TSS200_rmbl.saf \
-o Nanoseal_featureCounts_TSS_Q1.txt \
Ctrl_1.bam Ctrl_2.bam TETKO_1.bam TETKO_2.bam

featureCounts -T 11 \
-F SAF \
-p \
-a RNA-seq_Q2_mm10_filtered.TSS200_rmbl.saf \
-o Nanoseal_featureCounts_TSS_Q2.txt \
Ctrl_1.bam Ctrl_2.bam TETKO_1.bam TETKO_2.bam

featureCounts -T 11 \
-F SAF \
-p \
-a RNA-seq_Q3_mm10_filtered.TSS200_rmbl.saf \
-o Nanoseal_featureCounts_TSS_Q3.txt \
Ctrl_1.bam Ctrl_2.bam TETKO_1.bam TETKO_2.bam

featureCounts -T 11 \
-F SAF \
-p \
-a RNAbottom5to25percent.TSS200_rmbl.saf \
-o Nanoseal_featureCounts_TSS_bottom5to25percent.txt \
Ctrl_1.bam Ctrl_2.bam TETKO_1.bam TETKO_2.bam

featureCounts -T 11 \
-F SAF \
-p \
-a RNAbottom5percent.TSS200_rmbl.saf \
-o Nanoseal_featureCounts_TSS_RNAbottom5percent.txt \
Ctrl_1.bam Ctrl_2.bam TETKO_1.bam TETKO_2.bam

#---------------------TES +-200bp--------------------

featureCounts -T 11 \
-F SAF \
-p \
-a RNA-seq_Q1_mm10_filtered.TES200_rmbl.saf \
-o Nanoseal_featureCounts_TES_Q1.txt \
Ctrl_1.bam Ctrl_2.bam TETKO_1.bam TETKO_2.bam

featureCounts -T 11 \
-F SAF \
-p \
-a RNA-seq_Q2_mm10_filtered.TES200_rmbl.saf \
-o Nanoseal_featureCounts_TES_Q2.txt \
Ctrl_1.bam Ctrl_2.bam TETKO_1.bam TETKO_2.bam

featureCounts -T 11 \
-F SAF \
-p \
-a RNA-seq_Q3_mm10_filtered.TES200_rmbl.saf \
-o Nanoseal_featureCounts_TES_Q3.txt \
Ctrl_1.bam Ctrl_2.bam TETKO_1.bam TETKO_2.bam

featureCounts -T 11 \
-F SAF \
-p \
-a RNAbottom5to25percent.TES200_rmbl.saf \
-o Nanoseal_featureCounts_TES_bottom5to25percent.txt \
Ctrl_1.bam Ctrl_2.bam TETKO_1.bam TETKO_2.bam

featureCounts -T 11 \
-F SAF \
-p \
-a RNAbottom5percent.TES200_rmbl.saf \
-o Nanoseal_featureCounts_TES_RNAbottom5percent.txt \
Ctrl_1.bam Ctrl_2.bam TETKO_1.bam TETKO_2.bam


  #whole genome profile by genebody

computeMatrix scale-regions \
-S Ctrl_in_libnorm.bw \
TETKO_in_libnorm.bw  \
-R Mus_musculus.GRCm38.100_genes.bed \
--outFileName Nanoseal_global_10kb.mat.gz \
--regionBodyLength 10000 \
--beforeRegionStartLength 10000 \
--afterRegionStartLength 10000 \
--missingDataAsZero \
-p 11 \
-bl /mnt/bigfiles/CnR/2019_Oct/reference_data/mm10.blacklist.bed

plotProfile -m Nanoseal_global_10kb.mat.gz \
              --plotHeight 15 \
              --plotWidth 15 \
              -out Nanoseal_global_10kb.pdf \
              --outFileNameData Nanoseal_global_10kb \
              --samplesLabel "ctrl." "TET3KO" \
              --perGroup

#ploting according to the gene anotation, grouped by PolII ChIP-seq groups

computeMatrix scale-regions \
-S Ctrl_in_libnorm.bw \
TETKO_in_libnorm.bw \
-R ChIP-a.bed \
ChIP-b.bed \
ChIP-c.bed \
ChIP-d.bed \
--outFileName Nanoseal_byPolIIChIP_3kb_RmBl.mat.gz \
--regionBodyLength 5000 \
--beforeRegionStartLength 3000 \
--afterRegionStartLength 3000 \
--missingDataAsZero \
-p 11 \
-bl /mnt/bigfiles/CnR/2019_Oct/reference_data/mm10.blacklist.bed

plotProfile -m Nanoseal_byPolIIChIP_3kb_RmBl.mat.gz \
              --plotHeight 15 \
              --plotWidth 15 \
              -out Nanoseal_byPolIIChIP_3kb_RmBl.pdf \
              --outFileNameData Nanoseal_byPolIIChIP_3kb_RmBl \
              --regionsLabel "a" "b" "c" "d" \
              --samplesLabel "ctrl." "TET3KO" \
              --perGroup####use UCSC annotation of intron and exon from RefSeq mm10#####

cut -f 1,2,3,4,5,6 UCSC_mm10_refSeq_exons.rtf >UCSC_mm10_refSeq_exons.bed
cut -f 1,2,3,4,5,6 UCSC_mm10_refSeq_introns.rtf >UCSC_mm10_refSeq_introns.bed

#remove the exons and introns smaller than 200 bp

awk  '$3-$2 > 200 {print}' UCSC_mm10_refSeq_exons.bed | sort  -k1 -k2 > UCSC_mm10_refSeq_exons_bt200np.bed
wc -l UCSC_mm10_refSeq_exons.bed
#470749
wc -l UCSC_mm10_refSeq_exons_bt200np.bed
#115477


awk  '$3-$2 > 200 {print}' UCSC_mm10_refSeq_introns.bed | sort  -k1 -k2 > UCSC_mm10_refSeq_introns_bt200np.bed
wc -l UCSC_mm10_refSeq_introns.bed
#421571
wc -l UCSC_mm10_refSeq_introns_bt200np.bed
#363279


computeMatrix scale-regions \
-S Ctrl_in_libnorm.bw \
TETKO_in_libnorm.bw  \
-R UCSC_mm10_refSeq_exons_bt200np.bed \
--outFileName Nanoseal_global_ucsc_exon_bt200np.mat.gz \
--regionBodyLength 5000 \
--beforeRegionStartLength 200 \
--afterRegionStartLength 200 \
--missingDataAsZero \
-bl /mnt/bigfiles/CnR/2019_Oct/reference_data/mm10.blacklist.bed \
-p 11

plotProfile -m Nanoseal_global_ucsc_exon_bt200np.mat.gz \
              --plotHeight 15 \
              --plotWidth 15 \
              -out Nanoseal_global_ucsc_exon_bt200np.pdf \
              --outFileNameData Nanoseal_global_ucsc_exon_bt200np \
              --samplesLabel "ctrl." "TET3KO" \
              --perGroup

              
              
computeMatrix scale-regions \
-S Ctrl_in_libnorm.bw \
TETKO_in_libnorm.bw  \
-R UCSC_mm10_refSeq_introns_bt200np.bed \
--outFileName Nanoseal_global_ucsc_intron_bt200np.mat.gz \
--regionBodyLength 5000 \
--beforeRegionStartLength 200 \
--afterRegionStartLength 200 \
--missingDataAsZero \
-bl /mnt/bigfiles/CnR/2019_Oct/reference_data/mm10.blacklist.bed \
-p 11

plotProfile -m Nanoseal_global_ucsc_intron_bt200np.mat.gz \
              --plotHeight 15 \
              --plotWidth 15 \
              -out Nanoseal_global_ucsc_intron_bt200np.pdf \
              --outFileNameData Nanoseal_global_ucsc_intron_bt200np \
              --samplesLabel "ctrl." "TET3KO" \
              --perGroup#spurious gene and radomn non-spirous genes profile grouped by genotype

#generate random non-spurious expressed gene list with the same size of spurious gene group
cat non_spurious_gene.bed.csv |shuf | head -515 > non_spurious_gene_random1.bed
cat non_spurious_gene.bed.csv |shuf | head -515 > non_spurious_gene_random2.bed

computeMatrix scale-regions \
-S Ctrl_in_libnorm.bw \
TETKO_in_libnorm.bw \
-R 515_spurious_gene.bed \
non_spurious_gene_random1.bed \
non_spurious_gene_random2.bed \
--outFileName Nanoseal_bySpurious_gene_10kb_mm10.mat.gz \
--regionBodyLength 10000 \
--beforeRegionStartLength 10000 \
--afterRegionStartLength 10000 \
--missingDataAsZero \
-p 11 \
-bl /mnt/bigfiles/CnR/2019_Oct/reference_data/mm10.blacklist.bed

plotProfile -m Nanoseal_bySpurious_gene_10kb_mm10.mat.gz \
              --plotHeight 15 \
              --plotWidth 15 \
              -out Nanoseal_bySpurious_gene_10kb_mm10.pdf \
              --outFileNameData Nanoseal_bySpurious_gene_100kb_mm10 \
              --regionsLabel "spurious genes" "random geneset 1" "random geneset 2" \
              --samplesLabel "ctrl." "TET3KO" \
              --perGroup
              
#5hmC at spurious gene promotor, intron, and exon

##promotors
computeMatrix reference-point \
-S Ctrl_in_libnorm.bw \
TETKO_in_libnorm.bw  \
-R 515_spurious_gene.bed \
--outFileName Nanoseal_spurious_promoter.mat.gz \
--beforeRegionStartLength 2000 \
--afterRegionStartLength 2000 \
--missingDataAsZero \
-bl /mnt/bigfiles/CnR/2019_Oct/reference_data/mm10.blacklist.bed \
-p 11 


plotProfile -m Nanoseal_spurious_promoter.mat.gz \
              --plotHeight 15 \
              --plotWidth 15 \
              -out Nanoseal_spurious_promoter.pdf \
              --outFileNameData Nanoseal_spurious_promoter \
              --samplesLabel "ctrl." "TET3KO" \
              --perGroup
              

###exons and introns

###find the transcript annotation of spurious genes
##annotation file: EnsemblGeneToTranscripts.txt
sed -E 's/("([^"]*)")?,/\2\t/g' EnsemblGeneToTranscripts.txt > EnsemblGeneToTranscripts.tsv
sed -i '1d' EnsemblGeneToTranscripts.tsv

join -1 4 -2 1 <(sort -k 4 515_spurious_gene.bed.csv) <(sort EnsemblGeneToTranscripts.tsv) > 515_spurious_transcripts.csv

sed -E 's/("([^"]*)")?_/\2\t/g' UCSC_mm10_Ensembl_exons.txt > UCSC_mm10_Ensembl_exons_annotation_split.bed
sed -E 's/("([^"]*)")?_/\2\t/g' UCSC_mm10_Ensembl_introns.txt > UCSC_mm10_Ensembl_introns_annotation_split.bed

wc -l 515_spurious_transcripts.csv
#3220

join -1 9 -2 4 -o 2.1,2.2,2.3,2.4,2.5,2.6 <(sort -k 9 515_spurious_transcripts.csv) <(sort -k 4 UCSC_mm10_Ensembl_exons_annotation_split.bed) | sed "s/ /\t/g" > 515_spurious_exons.bed

join -1 9 -2 4 -o 2.1,2.2,2.3,2.4,2.5,2.6 <(sort -k 9 515_spurious_transcripts.csv) <(sort -k 4 UCSC_mm10_Ensembl_introns_annotation_split.bed) | sed "s/ /\t/g" > 515_spurious_introns.bed



#filter out regions smaller than 200bp

awk  '$3-$2 > 200 {print}' 515_spurious_exons.bed | sort  -k1 -k2 > 515_spurious_exons_bt200np.bed
awk  '$3-$2 > 200 {print}' 515_spurious_introns.bed | sort  -k1 -k2 > 515_spurious_introns_bt200np.bed

computeMatrix scale-regions \
-S Ctrl_in_libnorm.bw \
TETKO_in_libnorm.bw  \
-R 515_spurious_exons_bt200np.bed \
--outFileName Nanoseal_spurious_ucsc_exon_bt200np.mat.gz \
--regionBodyLength 5000 \
--beforeRegionStartLength 200 \
--afterRegionStartLength 200 \
--missingDataAsZero \
-bl /mnt/bigfiles/CnR/2019_Oct/reference_data/mm10.blacklist.bed \
-p 11

plotProfile -m Nanoseal_spurious_ucsc_exon_bt200np.mat.gz \
              --plotHeight 15 \
              --plotWidth 15 \
              -out Nanoseal_spurious_ucsc_exon_bt200np.pdf \
              --outFileNameData Nanoseal_spurious_ucsc_exon_bt200np \
              --samplesLabel "ctrl." "TET3KO" \
              --perGroup


computeMatrix scale-regions \
-S Ctrl_in_libnorm.bw \
TETKO_in_libnorm.bw  \
-R 515_spurious_introns_bt200np.bed \
--outFileName Nanoseal_spurious_ucsc_intron_bt200np.mat.gz \
--regionBodyLength 5000 \
--beforeRegionStartLength 200 \
--afterRegionStartLength 200 \
--missingDataAsZero \
-bl /mnt/bigfiles/CnR/2019_Oct/reference_data/mm10.blacklist.bed \
-p 11

plotProfile -m Nanoseal_spurious_ucsc_intron_bt200np.mat.gz \
              --plotHeight 15 \
              --plotWidth 15 \
              -out Nanoseal_spurious_ucsc_intron_bt200np.pdf \
              --outFileNameData Nanoseal_spurious_ucsc_intron_bt200np \
              --samplesLabel "ctrl." "TET3KO" \
              --perGroup



