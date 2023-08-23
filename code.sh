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

# Count features: 
featureCounts -p -a refs/features.gff -o counts.txt bam/BORED_1.bam

# Show the first five lines of the counts.txt file:
cat counts.txt | head -5

# Standardizing the count matrix:
featureCounts -p -a refs/features.gff -o counts.txt \
    reads/BORED_1.bam reads/BORED_2.bam reads/BORED_3.bam \
    reads/EXCITED_1.bam reads/EXCITED_2.bam reads/EXCITED_3.bam

# Classification based RNA-Seq with kallisto and salmon

#Create and activate a new environment:
mamba create -y -n salmon
conda activate salmon
# Install the software:
mamba install salmon parallel

#Prepare the transcriptome for classification:

REF=refs/transcripts.fa

# The  salmon index:
IDX=idx/salmon.idx

# Build the index with salmon.
salmon index -t ${REF} -i ${IDX}

#Running the classification:
mkdir -p salmon

# Run a salmon quantification.
salmon quant -i ${IDX} -l A --validateMappings -1 reads/BORED_1_R1.fq -2 reads/BORED_1_R2.fq  -o salmon/BORED_1
cat salmon/BORED_1/quant.sf | head  | column -t

# automate the process:
cat reads/ids.txt | parallel -j 4 "salmon quant -i ${IDX} -l A --validateMappings -1 reads/{}_R1.fq -2 reads/{}_R2.fq  -o salmon/{}"

#combine de counts:
Rscript code/combine_transcripts.r                        #the script combine_transcripts.r is available in this repository.

# Differential analysis with DESeq2:
RScript code/deseq2.R                                     #the script combine_transcripts.r is available in this repository.





















