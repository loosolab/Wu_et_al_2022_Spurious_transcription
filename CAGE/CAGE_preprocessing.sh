REMOVE THIS LATER:
yua35_Fan_Cage_seq_japan



## Basic CAGE analysis was peformed by the RIKEN Institute Japan using the published RECLU pipeline (RECLU: a pipeline to discover reproducible transcriptional start sites and their alternative regulation using capped analysis of gene expression (CAGE), Ohmiya et. al., BMC Genomics. 2014 Apr 25;15:269. doi: 10.1186/1471-2164-15-269)
## code repository: https://osdn.net/projects/reclu/releases/


##CAGE processing
#acquired mapped bam files from RILKEN 
#remove the black list region
find . -name '*.bam*' | parallel 'bedtools intersect -v -b mm10.blacklist.bed -abam {} > {.}.filtered.bam'

#sort bam files
sambamba sort -t 11 -o TET3KO-1_filtered_sorted.bam TET3KO-1_filtered.bam
sambamba sort -t 11 -o TET3KO-2_filtered_sorted.bam TET3KO-2_filtered.bam
sambamba sort -t 11 -o Ctrl-1_filtered_sorted.bam Ctrl-1_filtered.bam
sambamba sort -t 11 -o Ctrl-2_filtered_sorted.bam Ctrl-2_filtered.bam

#conovert bam files to bigwig
bamCoverage -b TET3KO-1_filtered_sorted.bam -o TET3KO-1_filtered.bw -of bigwig -bs 10 -p 11 --normalizeUsing CPM
bamCoverage -b TET3KO-2_filtered_sorted.bam -o TET3KO-2_filtered.bw -of bigwig -bs 10 -p 11 --normalizeUsing CPM
bamCoverage -b Ctrl-1_filtered_sorted.bam -o Ctrl-1_filtered.bw -of bigwig -bs 10 -p 11 --normalizeUsing CPM
bamCoverage -b Ctrl-2_filtered_sorted.bam -o Ctrl-2_filtered.bw -of bigwig -bs 10 -p 11 --normalizeUsing CPM

#merge two replicates into one for visualisation
wiggletools mean Ctrl-1_filtered.bw Ctrl-2_filtered.bw | wigToBigWig stdin mm10.chrom.sizes Ctrl_mean.bw
wiggletools mean TET3KO-1_filtered.bw TET3KO-2_filtered.bw | wigToBigWig stdin mm10.chrom.sizes TET3KO_mean.bw


#visualisation of CAGE signal at spurious and non-spurious genes
computeMatrix reference-point \
-p 11 \
-S Ctrl_mean.bw \
-R 515_spurious_gene.bed.csv non_spurious_gene.bed.csv \
--outFileName matrix_CAGE_TSS_spurious_nonspurious.gz \
-a 1000 \
-b 1000

plotProfile -m matrix_CAGE_TSS_spurious_nonspurious.gz \
      --plotFileFormat svg \
      -out CAGE_TSS_spurious_nonspurious.svg 
      
      
#make_ctss.sh acquired from Takahashi, H., Lassmann, T., Murata, M. et al. 5′ end–centered expression profiling using cap-analysis gene expression and next-generation sequencing. Nat Protoc 7, 542–561 (2012). https://doi.org/10.1038/nprot.2012.005
#clip the tags to 1 nucleotide length at 5'end, outpout bed files
sed -i 's/\r//' make_ctss.sh
bash make_ctss.sh -q 10 TET3KO-1.bam TET3KO-2.bam Ctrl-1.bam Ctrl-2.bam


# motif enrichment by the method of Neri, F., Rapelli, S., Krepelova, A. et al. Intragenic DNA methylation prevents spurious transcription initiation. Nature 543, 72–77 (2017). https://doi.org/10.1038/nature21373


awk '{ total += $4; count++ } END { print total/count }' Ctrl-1.bam.ctss
#8.90398
awk '{ total += $4; count++ } END { print total/count }' Ctrl-2.bam.ctss
#8.07687
awk '{ total += $4; count++ } END { print total/count }' TET3KO-1.bam.ctss
#8.66748
awk '{ total += $4; count++ } END { print total/count }' TET3KO-2.bam.ctss
#8.47664

# the average CAGE tag number at a CTSS is 8, use 8 as the cut-off. CTSS with higher than 8 CAGE tags are considered as highly epxressed CTSS
awk '{ if ($4 >= 8) { print } }' Ctrl-1.bam.ctss > Ctrl-1.bam.ctss_cutoff_8.bed
awk '{ if ($4 >= 8) { print } }' Ctrl-2.bam.ctss > Ctrl-2.bam.ctss_cutoff_8.bed
awk '{ if ($4 >= 8) { print } }' TET3KO-1.bam.ctss > TET3KO-1.bam.ctss_cutoff_8.bed
awk '{ if ($4 >= 8) { print } }' TET3KO-2.bam.ctss > TET3KO-2.bam.ctss_cutoff_8.bed

#merge the CTSS in each replicates
bedops --merge Ctrl-1.bam.ctss_cutoff_8.bed Ctrl-2.bam.ctss_cutoff_8.bed > Ctrl.bam.ctss_cutoff_8.bed
bedops --merge TET3KO-1.bam.ctss_cutoff_8.bed TET3KO-2.bam.ctss_cutoff_8.bed > TET3KO.bam.ctss_cutoff_8.bed

#find the TET3 knockout specific CTSS with range of +-25bp
bedops --range 25 --everything Ctrl.bam.ctss_cutoff_8.bed > Ctrl.bam.ctss_cutoff_8_25bp.bed
bedops --not-element-of 1 Ctrl.bam.ctss_cutoff_8_25bp.bed TET3KO.bam.ctss_cutoff_8.bed > TET3KO_specific_CTSS_filter8.bed

sort -k1,1 -k2,2n TET3KO_specific_CTSS_filter8.bed > TET3KO_specific_CTSS_filter8_sorted.bed

mkdir homer_motif_TETE3KO_filterby_8_sn_25
LC_ALL=C findMotifsGenome.pl TET3KO_specific_CTSS_filter8.bed mm10 homer_motif_TETE3KO_filterby_8_sn_25/ -size 25 -p 11


# annotate putative motif location
cd /homer_motif_TETE3_KO_filterby_8_sn_25

scanMotifGenomeWide.pl known17.motif mm10 -bed > known17.Sp2.motif.mm10.bed
scanMotifGenomeWide.pl known24.motif mm10 -bed > known24.pu.1.motif.mm10.bed


#visualise the CTSS in genome track

#merge the replilcates in Control group
awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' Ctrl-1.bam.ctss.bed > Ctrl-1_CTSS.bedgraph
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s Ctrl-1_CTSS.bedgraph > Ctrl-1_CTSS.bedgraph.tmp
bedtools merge -i Ctrl-1_CTSS.bedgraph.tmp -c 4 -d 0 -o max > Ctrl-1_CTSS_out.bedgraph
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s Ctrl-1_CTSS_out.bedgraph > Ctrl-1_CTSS_sorted.bedgraph
bedGraphToBigWig Ctrl-1_CTSS_sorted.bedgraph mm10.chrom.sizes Ctrl-1_CTSS.bw

awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' Ctrl-2.bam.ctss.bed > Ctrl-2_CTSS.bedgraph
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s Ctrl-2_CTSS.bedgraph > Ctrl-2_CTSS.bedgraph.tmp
bedtools merge -i Ctrl-2_CTSS.bedgraph.tmp -c 4 -d 0 -o max > Ctrl-2_CTSS_out.bedgraph
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s Ctrl-2_CTSS_out.bedgraph > Ctrl-2_CTSS_sorted.bedgraph
bedGraphToBigWig Ctrl-2_CTSS_sorted.bedgraph mm10.chrom.sizes Ctrl-2_CTSS.bw

wiggletools mean Ctrl-1_CTSS.bw Ctrl-2_CTSS.bw | wigToBigWig stdin mm10.chrom.sizes Ctrl_CTSS_mean.bw

#merge the replilcates in TET3KO group
awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' TET3KO-1.bam.ctss.bed > TET3KO-1_CTSS.bedgraph
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s TET3KO-1_CTSS.bedgraph > TET3KO-1_CTSS.bedgraph.tmp
bedtools merge -i TET3KO-1_CTSS.bedgraph.tmp -c 4 -d 0 -o max > TET3KO-1_CTSS_out.bedgraph
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s TET3KO-1_CTSS_out.bedgraph > TET3KO-1_CTSS_sorted.bedgraph
bedGraphToBigWig TET3KO-1_CTSS_sorted.bedgraph mm10.chrom.sizes TET3KO-1_CTSS.bw

awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' TET3KO-2.bam.ctss.bed > TET3KO-2_CTSS.bedgraph
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s TET3KO-2_CTSS.bedgraph > TET3KO-2_CTSS.bedgraph.tmp
bedtools merge -i TET3KO-2_CTSS.bedgraph.tmp -c 4 -d 0 -o max > TET3KO-2_CTSS_out.bedgraph
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s TET3KO-2_CTSS_out.bedgraph > TET3KO-2_CTSS_sorted.bedgraph
bedGraphToBigWig TET3KO-2_CTSS_sorted.bedgraph mm10.chrom.sizes TET3KO-2_CTSS.bw

wiggletools mean TET3KO-1_CTSS.bw TET3KO-2_CTSS.bw | wigToBigWig stdin mm10.chrom.sizes TET3KO_CTSS_mean.bw

#merge the replilcates in Control group with 8 tags as cut-off
awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' Ctrl-1.bam.ctss_cutoff_8.bed > Ctrl-1_CTSS_cutoff_8.bedgraph
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s Ctrl-1_CTSS_cutoff_8.bedgraph > Ctrl-1_CTSS_cutoff_8.bedgraph.tmp
bedtools merge -i Ctrl-1_CTSS_cutoff_8.bedgraph.tmp -c 4 -d 0 -o max > Ctrl-1_CTSS_cutoff_8_out.bedgraph
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s Ctrl-1_CTSS_cutoff_8_out.bedgraph > Ctrl-1_CTSS_cutoff_8_sorted.bedgraph
bedGraphToBigWig Ctrl-1_cutoff_8_CTSS_sorted.bedgraph mm10.chrom.sizes Ctrl-1_cutoff_8_CTSS.bw

awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' Ctrl-2.bam.ctss_cutoff_8.bed > Ctrl-2_CTSS_cutoff_8.bedgraph
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s Ctrl-2_CTSS_cutoff_8.bedgraph > Ctrl-2_CTSS_cutoff_8.bedgraph.tmp
bedtools merge -i Ctrl-2_CTSS_cutoff_8.bedgraph.tmp -c 4 -d 0 -o max > Ctrl-2_CTSS_cutoff_8_out.bedgraph
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s Ctrl-2_CTSS_cutoff_8_out.bedgraph > Ctrl-2_CTSS_cutoff_8_sorted.bedgraph
bedGraphToBigWig Ctrl-2_cutoff_8_CTSS_sorted.bedgraph mm10.chrom.sizes Ctrl-2_cutoff_8_CTSS.bw

wiggletools mean Ctrl-1_cutoff_8_CTSS.bw Ctrl-2_cutoff_8_CTSS.bw | wigToBigWig stdin mm10.chrom.sizes Ctrl_CTSS_cutoff_8_mean.bw

#merge the replilcates in TET3KO group with 8 tags as cut-off
awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' TET3KO-1.bam.ctss_cutoff_8.bed > TET3KO-1_CTSS_cutoff_8.bedgraph
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s TET3KO-1_CTSS_cutoff_8.bedgraph > TET3KO-1_CTSS_cutoff_8.bedgraph.tmp
bedtools merge -i TET3KO-1_CTSS_cutoff_8.bedgraph.tmp -c 4 -d 0 -o max > TET3KO-1_CTSS_cutoff_8_out.bedgraph
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s TET3KO-1_CTSS_cutoff_8_out.bedgraph > TET3KO-1_CTSS_cutoff_8_sorted.bedgraph
bedGraphToBigWig TET3KO-1_cutoff_8_CTSS_sorted.bedgraph mm10.chrom.sizes TET3KO-1_cutoff_8_CTSS.bw

awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$5}' TET3KO-2.bam.ctss_cutoff_8.bed > TET3KO-2_CTSS_cutoff_8.bedgraph
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s TET3KO-2_CTSS_cutoff_8.bedgraph > TET3KO-2_CTSS_cutoff_8.bedgraph.tmp
bedtools merge -i TET3KO-2_CTSS_cutoff_8.bedgraph.tmp -c 4 -d 0 -o max > TET3KO-2_CTSS_cutoff_8_out.bedgraph
LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -s TET3KO-2_CTSS_cutoff_8_out.bedgraph > TET3KO-2_CTSS_cutoff_8_sorted.bedgraph
bedGraphToBigWig TET3KO-2_cutoff_8_CTSS_sorted.bedgraph mm10.chrom.sizes TET3KO-2_cutoff_8_CTSS.bw

wiggletools mean TET3KO-1_cutoff_8_CTSS.bw TET3KO-2_cutoff_8_CTSS.bw | wigToBigWig stdin mm10.chrom.sizes TET3KO_CTSS_cutoff_8_mean.bw


# calculate the CTSS at genbody without exon1
# Mus_musculus.GRCm38.100_genewithoutexon1.bed
# Count overlap of CTSS on genebody without exon1. Only CTSS at the same strand with gene annotation (sense CTSS) are kept.
find . -name '*ctss.bed' | parallel 'bedtools coverage -counts -s -a Mus_musculus.GRCm38.100_genewithoutexon1.bed.csv -b {} > {.}_quantif_withoutexon1.sense.bed'

#visualisation proceeded with R


#CAGE signal according to the PolII pSer5 ChIP-seq group: a, b, c, d defined in Fig2b 
computeMatrix reference-point \
-p 11 \
-S Ctrl_mean.bw \
-R ChIP-a.bed ChIP-b.bed ChIP-c.bed ChIP-d.bed \
--outFileName matrix_CAGE_TSS_toChIPgroup.gz \
-a 1000 \
-b 1000 

plotProfile -m matrix_CAGETSS_toChIPgroup.gz \
      --plotFileFormat svg \
      -out matrix_CAGETSS_toChIPgroup.svg
      

