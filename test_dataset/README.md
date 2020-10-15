# Gencode
The current version is v35

## GTF file
```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz
```

## Genome file
Nucleotide sequence of the GRCh38.p13 genome assembly version on all regions, including reference chromosomes, scaffolds, assembly patches and haplotypes
```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/GRCh38.p13.genome.fa.gz
```
## mRNAs file
Nucleotide sequences of all transcripts on the reference chromosomes

```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.transcripts.fa.gz
```


## Generation of test data

Genomic region  chr3:52287089-52809009 (around BAP1 gene)

```
#pick the fasta region
gzip -dc GRCh38.p13.genome.fa.gz > GRCh38.p13.genome.fa
samtools faidx GRCh38.p13.genome.fa chr3:52287089-52809009  -o GRCh38.chr3_52287089_52809009.fa
#pick the gene models from gtf
gzip -dc gencode.v35.annotation.gtf.gz | awk '{if($1 == "chr3" && $4 >= 52287089 && $5 <= 52809009){print $0}}' > gencode.v35.chr3_52287089_52809009.gtf
grep "	gene	"  gencode.v35.chr3_52287089_52809009.gtf | grep "protein_coding" | awk '{print $10}' > protein_coding_genes.txt
awk 'BEGIN{while(getline <"protein_coding_genes.txt"){a[$1]=1}}{if($10 in a){print $0}}' gencode.v35.chr3_52287089_52809009.gtf > gencode.v35.chr3_52287089_52809009.protein_coding.gtf
grep "	transcript	" gencode.v35.chr3_52287089_52809009.protein_coding.gtf | awk '{print $10" "$12}' | sort -u -k1,1 | sed 's/"//g;s/;//g' | awk '{print $2}'| sort -u  > protein_coding_transcripts.txt
#format transcript file
gzip -dc gencode.v35.transcripts.fa.gz | awk '{if($0 ~/>/){split($0,a,"|"); print a[1]" "a[2]}else{print $0}}' > gencode.v35.transcripts.fa
#index
samtools faidx gencode.v35.transcripts.fa -r protein_coding_transcripts.txt > gencode.v35.chr3_52287089_52809009.transcripts.fa
#we have to donwload wgsim
git clone https://github.com/lh3/wgsim.git
cd wgsim
gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm
cd ..
#we generate 10 samples with upto 10 fusion transcripts,it do not take into account the transcript orientation
perl generate_chimeric.pl -a gencode.v35.chr3_52287089_52809009.transcripts.fa -b T000
#we have to rescale the coordinates
awk -F"\t" 'BEGIN{OFS="\t";}{$4-=52287088; $5-=52287088; print $0}' gencode.v35.chr3_52287089_52809009.protein_coding.gtf > gencode.v35.chr3_52287089_52809009.protein_coding_rescaled.gtf
#we link the test genome and gtf annotation
ln -s gencode.v35.chr3_52287089_52809009.protein_coding_rescaled.gtf genome.gtf
ln -s GRCh38.chr3_52287089_52809009.fa genome.fa
```

