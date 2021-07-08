#script to run to extract SRA code from bam list of individuals in the reference panel

for F in $(cat ../data/bam1.list)
do
    basename $F _dedup_overlapclipped.bam >> ../data/identifiers.txt
done