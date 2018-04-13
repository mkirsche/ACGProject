REF='chr22.fa.inserts'
VARIANTS='chr22.fa.10.vcf'
BAMFILE='aln.bam'

javac *.java
rm -r $OUTDIR
mkdir $OUTDIR
flank=50
while read p || [[ -n $p ]]; do
  if [[ ${p:0:1} == "#" ]] ; then continue; fi
  tokens=( $p )
  chr=${tokens[0]}
  pos=${tokens[1]}
  echo $chr:$pos
  info=${tokens[7]}
  lenSubstring=`echo $info | grep -o "SVLEN=.*"`
  len=${lenSubstring:6}
  samtools view -h $BAMFILE $chr:$(($pos - $flank))-$(($pos + $flank)) > $chr.$pos.sam
  # TODO Run java program TestEnds to get the soft clipped reads and see if they align to the ends of the sequence
done < $VARIANTS

