#usage: ./16S_extraction.sh gff fna
#gfffile is the file generated from barrnap with 16S positions
#fastafile is the fasta file with the genome sequence to extract 16S rRNA sequence from
#_dfiles ="gff/*.gff"
INPUT_DIR=$gff
for f in $INPUT_DIR/*.gff
do
	grep '16S' "$f" > "${$f}"_16S.gff   # update signature
done
