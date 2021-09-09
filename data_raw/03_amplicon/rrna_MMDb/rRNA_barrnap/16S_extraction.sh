#usage: ./16S_extraction.sh gff fna
#gfffile is the file generated from barrnap with 16S positions
#fastafile is the fasta file with the genome sequence to extract 16S rRNA sequence from


gfffile=$1
fastafile=$2

        grep '16S' $gfffile > 16S-gff.gff;
        bedtools getfasta -fi $fastafile -bed 16S-gff.gff -fo 16S-fasta.fna;
        grep -m 1 ">" 16S-fasta.fna|sed 's/>//g' > 16S-id.txt;
        xargs samtools faidx 16S-fasta.fna < 16S-id.txt > $fastafile-16S.fna

rm 16S-gff.gff
rm 16S-fasta.fna
rm 16S-fasta.fna.fai
rm 16S-id.txt
        exit 0;