for F in $(cat data/bam1.list)
do
    TOTDEPTH=$(samtools view -H $F | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}')
    DEPTH=$(samtools depth -a $F |  awk '{sum+=$3} END { print "Average = ",sum/$TOTDEPTH}')
    printf $DEPTH >> try_depth2.txt
done