# TODO make these values command line arguments
#REF='chr22.fa'
#VARIANTS='chr22.fa.10.vcf'
#FASTQ1='/home/mkirsche/chr22.fa.10.1'
#FASTQ2='/home/mkirsche/chr22.fa.10.2'
#BAMFILE='inserts.bam'

REF=$1
VARIANTS=$2
FASTQ1=$3
FASTQ2=$4
BAMFILE=$5

javac *.java
java FastaInsert $REF $VARIANTS $REF.inserts

bwa index $REF
bwa mem -o aln.sam $REF $FASTQ1 $FASTQ2

samtools view -h -b aln.sam > aln_unsorted.bam
samtools sort aln_unsorted.bam > aln.bam
samtools index aln.bam
