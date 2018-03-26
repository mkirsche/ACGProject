REF=$1
VARIANTS=$2
FASTQ1=$3
FASTQ2=$4
BAMFILE='inserts.bam'

javac *.java
java FastaInsert $REF $VARIANTS $REF.inserts

bwa index $REF.inserts
bwa mem -o inserts.sam $REF.inserts $FASTQ1 $FASTQ2

samtools view -h -b inserts.sam > ref_unsorted.bam
samtools sort ref_unsorted.bam > $BAMFILE
samtools index $BAMFILE
