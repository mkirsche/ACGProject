REF='chr22.fa.inserts'
VARIANTS='chr22.fa.10.vcf'
BAMFILE='inserts.bam'
OUTDIR='output'
NEWVCF=$VARIANTS'.filtered'

javac *.java
rm -r $OUTDIR
mkdir $OUTDIR
insPos=0
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
  end=$(($insPos + $len))
  samtools view -h $BAMFILE $chr:$(($pos - $flank))-$(($pos + $flank)) > $chr.$pos.sam
  samtools view -b $chr.$pos.sam > $chr.$pos.bam
  samtools view -b $BAMFILE insertion:$insPos-$end > $chr.$pos.ins.bam
  samtools merge $chr.$pos.all.bam $chr.$pos.bam $chr.$pos.ins.bam
  rm $chr.$pos.bam
  rm $chr.$pos.ins.bam
  samtools bam2fq $chr.$pos.all.bam > $OUTDIR/$chr.$pos.fq
  rm $chr.$pos.all.bam
  # Now I have relevant reads in fq file
  #echo 'Running spades'
  spades.py -s output/$chr.$pos.fq -o $OUTDIR/$chr.$pos -t 8 &> /dev/null
  cat $OUTDIR/$chr.$pos'/contigs.fasta' | bioawk -c fastx '{ print length($seq), $name }' | sort -k1,1rn
  longest=( `cat $OUTDIR/$chr.$pos'/contigs.fasta' | bioawk -c fastx '{ print length($seq), $seq }' | sort -k1,1rn | head -1` )
  assemblyInsert=${longest[1]}
  before=( `samtools faidx $REF $chr:$(($pos - $flank))-$pos` )
  after=( `samtools faidx $REF $chr:$(($pos + 1))-$(($pos + $flank))` )
  samtools faidx $REF 'insertion':"$insPos"-"$end" > $OUTDIR/$chr.$pos'/ref.fa'
  insert=`cat $OUTDIR/$chr.$pos'/ref.fa' | bioawk -c fastx '{print $seq}'`
  refInsert=${before[1]}${insert}${after[1]}
  refInsert=`echo $refInsert | awk '{print toupper($0)}'`
  echo '>test' > $OUTDIR/$chr.$pos'/ref.fa'
  echo $refInsert>> $OUTDIR/$chr.$pos'/ref.fa'
  #echo 'Reference: '
  #echo $refInsert
  java Containment $refInsert $OUTDIR/$chr.$pos'/contigs.fasta'
done < $VARIANTS

