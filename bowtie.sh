REF=$1
READS=$2

echo "Reference is" $REF 
echo "Reads are" $READS

#bowtie2-2.2.5/bowtie2-build $REF $REF-bowtie-index
bowtie2-2.2.5/bowtie2 -x $REF-bowtie-index -U $READS -S $REF.$READS.sam
samtools view -Sb $REF.$READS.sam > $REF.$READS.bam
#https://broadinstitute.github.io/picard/explain-flags.htmls
samtools view -c -f 4 $REF.$READS.bam > reads-unmapped.count.txt
samtools view -c -F 4 $REF.$READS.bam > reads-mapped.count.txt
samtools sort $REF.$READS.bam $REF.$READS.bam.sorted
#SRR172903.NC005008.bam.sorted.bam
samtools index $REF.$READS.bam.sorted.bam
samtools idxstats $REF.$READS.bam.sorted.bam > reads.by.contigs.txt