
#!/bin/sh
#make_ctss.sh acquired from Takahashi, H., Lassmann, T., Murata, M. et al. 5′ end–centered expression profiling using cap-analysis gene expression and next-generation sequencing. Nat Protoc 7, 542–561 (2012). https://doi.org/10.1038/nprot.2012.005
if [ $# -eq 0 ]                  
then
 cat <<EOF
Usage is : $0 -q <mapping quality cutoff> <map1.bam> <map2.bam> ... 
EOF
 exit 1;
fi

QCUT= 

while getopts q: opt
do
case ${opt} in
q) QCUT=${OPTARG};;
*) usage;;
esac
done

if [ "${QCUT}" = "" ]; then QCUT=10; fi

for var in "$@"
do
if [[ $var =~ bam$ ]]; then
foo=$var
file=${foo##*/}
base=${file%%.*}
echo working on: $base

TMPFILE="/tmp/$(basename $0).$RANDOM.txt"
   samtools view  -F 4 -u -q $QCUT -b $var |  bamToBed -i stdin > $TMPFILE
   cat  ${TMPFILE} \
| awk 'BEGIN{OFS="\t"}{if($6=="+"){print $1,$2,$5}}' \
| sort -k1,1 -k2,2n \
| groupBy -i stdin -g 1,2 -c 3 -o count \
| awk -v x="$base" 'BEGIN{OFS="\t"}{print $1,$2,$2+1,  x  ,$3,"+"}' >> $var.ctss.bed

cat  ${TMPFILE} \
| awk 'BEGIN{OFS="\t"}{if($6=="-"){print $1,$3,$5}}' \
| sort -k1,1 -k2,2n \
| groupBy -i stdin -g 1,2 -c 3 -o count \
| awk -v x="$base" 'BEGIN{OFS="\t"}{print $1,$2-1,$2, x  ,$3,"-"}' >> $var.ctss.bed

rm $TMPFILE
fi
done


