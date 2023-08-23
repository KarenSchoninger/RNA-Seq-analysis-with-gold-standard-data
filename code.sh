# In the downloaded data folder there are the following files: features.gff, genome.fa and transcripts.fa.

# The reference genome:
IDX=refs/genome.fa

# Build the genome index:
hisat2-build $IDX $IDX

# Index the reference genome with samtools:
samtools faidx $IDX

# Generate the alignments for BORED:
hisat2 -x refs/genome.fa -1 reads/BORED_1_R1.fq -2 reads/BORED_1_R2.fq | samtools sort > reads/BORED_1.bam
hisat2 -x refs/genome.fa -1 reads/BORED_2_R1.fq -2 reads/BORED_2_R2.fq | samtools sort > reads/BORED_2.bam
hisat2 -x refs/genome.fa -1 reads/BORED_3_R1.fq -2 reads/BORED_3_R2.fq | samtools sort > reads/BORED_3.bam

# Generate the alignments for EXCITED:
hisat2 -x refs/genome.fa -1 reads/EXCITED_1_R1.fq -2 reads/EXCITED_1_R2.fq | samtools sort > reads/EXCITED_1.bam
hisat2 -x refs/genome.fa -1 reads/EXCITED_2_R1.fq -2 reads/EXCITED_2_R2.fq | samtools sort > reads/EXCITED_2.bam
hisat2 -x refs/genome.fa -1 reads/EXCITED_3_R1.fq -2 reads/EXCITED_3_R2.fq | samtools sort > reads/EXCITED_3.bam

# To visualize the alignments in IGV:
cat reads/ids.txt | parallel "bedtools genomecov -ibam bam/{}.bam -split -bg  > bam/{}.bg"
cat reads/ids.txt | parallel "bedGraphToBigWig bam/{}.bg  ${IDX}.fai bam/{}.bw"




