#### EXERCISE 1 #####

# user: testing / passwd: testing

## create directory unix_training
mkdir unix_training

## move into directory and create two directories
cd unix_training/
mkdir data
mkdir results

## copy paste the fastq files to your directory
cd data
cp /home/itg.be/pmonsieurs/data/fastq/* ./

## check the file size of the fastq files
du HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset_1.fastq.gz
du -sh HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset_1.fastq.gz
du -sh *.fastq.gz

## inspect fastq file
less HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset_1.fastq.gz

## zip and unzip a fastq file, and do less or wc
gunzip HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset_1.fastq.gz
less HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset_1.fastq
wc HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset_1.fastq
gzip HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset_1.fastq


#### EXERCISE 2 ####

## print the help file for BWA
bwa 
bwa index
bwa mem

## downloading reference genome via wget and index for BWA
cd ~/unix_training/data/
wget https://plasmodb.org/common/downloads/Current_Release/PvivaxP01/fasta/data/PlasmoDB-55_PvivaxP01_Genome.fasta
bwa index PlasmoDB-55_PvivaxP01_Genome.fasta

## run default BWA (no seed and threads added)
bwa mem -R "@RG\tID:001\tSM:001\tPL:ILLUMINA" \
        -o ../results/HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.sam \
        PlasmoDB-55_PvivaxP01_Genome.fasta \
        HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset_R1.fastq.gz \
        HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset_R2.fastq.gz
samtools sort -@2 HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.sam -o HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.bam

## run default BWA with piping


## run BWA with a seed of 30 & threads 2 and write results to the 
## "results" directory
bwa mem -R "@RG\tID:001\tSM:001\tPL:ILLUMINA" \
        -k 30 \
        -t 2 \
        -o ../results/HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.sam \
        PlasmoDB-55_PvivaxP01_Genome.fasta \
        HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset_R1.fastq.gz \
        HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset_R2.fastq.gz
samtools sort -@2 HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.sam -o HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.bam

## do the same but now pipe the output, so that you do not 
## have to write a big same file but rather write directly
## to a .bam file
bwa mem -R "@RG\tID:001\tSM:001\tPL:ILLUMINA" \
        -k 30 \
        -t 2 \
        PlasmoDB-55_PvivaxP01_Genome.fasta \
        HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset_R1.fastq.gz \
        HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset_R2.fastq.gz \
        | samtools sort -@2 -o ../results/HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.bam -


## defined the R1 fastq-file as a variable to 
## partially automatize the process
fastq_file_R1=HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset_R1.fastq.gz
fastq_file_R2=${fastq_file_R1/_R1.fastq.gz/_R2.fastq.gz}
bam_file=${fastq_file_R1/_R1.fastq.gz/.bam}

echo $fastq_file_R1
echo $fastq_file_R2
echo $bam_file

## now use those variable to run the BWA command
bwa mem -R "@RG\tID:001\tSM:001\tPL:ILLUMINA" \
        -k 30 \
        -t 2 \
        PlasmoDB-55_PvivaxP01_Genome.fasta \
        ${fastq_file_R1} \
        ${fastq_file_R2} \
        | samtools sort -@2 -o ../results/${bam_file} -



## you can run this for all 5 pairs of fastq files manual
## or you can use a for loop
for fastq_file_R1 in *_R1.fastq.gz
do
    file_prefix=${fastq_file_R1%_R1.fastq.gz}
    fastq_file_R2=${fastq_file_R1/_R1.fastq.gz/_R2.fastq.gz}
    bam_file=${fastq_file_R1/_R1.fastq.gz/.bam}
    bwa mem -R "@RG\tID:001\tSM:001\tPL:ILLUMINA" \
        -k 30 \
        -t 2 \
        PlasmoDB-55_PvivaxP01_Genome.fasta \
        $fastq_file_R1 \
        $fastq_file_R2 \
        | samtools sort -@2 -o ../results/${bam_file} -
done




## remove duplicates from the .bam file
cd ../results/
java -jar /home/itg.be/pmonsieurs/software/picard/picard.jar \
     MarkDuplicates REMOVE_DUPLICATES=true \
     I=HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset.bam \
     O=HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset.markdups.bam \
     M=HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset.markdups_metrics.txt

## or again, do this in a for loop
for bam_file in *.subset.bam
do
    file_prefix=${bam_file%.bam}
    java -jar /home/itg.be/pmonsieurs/software/picard/picard.jar \
        MarkDuplicates REMOVE_DUPLICATES=true \
        I=${bam_file} \
        O=${file_prefix}.markdups.bam \
        M=${file_prefix}.markdups_metrics.txt
done

## optional: select proper paired reads 
samtools view -@2 -bf 0x2 \
    HJMTYDSX2_104610-001-040_GTCAGTTG-CGATTGGA_L004.subset.markdups.bam \
    > HJMTYDSX2_104610-001-040_GTCAGTTG-CGATTGGA_L004.subset.markdups.proper_paired.bam

for bam_file in *.markdups.bam
do
    file_prefix=${bam_file%.bam}
    samtools view -@2 -bf 0x2 \
        ${bam_file} \
        > ${file_prefix}.proper_paired.bam
done


## indexing of bam files
samtools index HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset.bam
for bam_file in *.bam
do
    samtools index $bam_file
done

## run flagstat to get some aligning statistics
samtools flagstat HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.subset.bam




#### EXERCISE 3 ####

## index the reference genome
samtools faidx ../data/PlasmoDB-55_PvivaxP01_Genome.fasta

## run old-fashioned
gatk HaplotypeCaller \
    -R ../data/PlasmoDB-55_PvivaxP01_Genome.fasta \
    -I HJMTYDSX2_104610-001-040_GTCAGTTG-CGATTGGA_L004.subset.markdups.proper_paired.bam \
    -O HJMTYDSX2_104610-001-040_GTCAGTTG-CGATTGGA_L004.vcf \
    --native-pair-hmm-threads 2

## run new approach for unified genotyping
gatk HaplotypeCaller \
    -R ../data/PlasmoDB-55_PvivaxP01_Genome.fasta \
    -I HJMTYDSX2_104610-001-040_GTCAGTTG-CGATTGGA_L004.subset.markdups.proper_paired.bam \
    -O HJMTYDSX2_104610-001-040_GTCAGTTG-CGATTGGA_L004.vcf \
    -ERC GVCF \
    --native-pair-hmm-threads 2

## example of variables
bam_file=HJMTYDSX2_104610-001-040_GTCAGTTG-CGATTGGA_L004.subset.markdups.proper_paired.bam
echo $bam_file

## example of a substitution
vcf_file=${bam_file/.subset.markdups.proper_paired.bam/.g.vcf}
echo $vcf_file

## example of a for loop
for bam_file in *.proper_paired.bam
do
    echo ${bam_file}
done

## example of a for loop with substitution
for bam_file in *.proper_paired.bam
do
    echo ${bam_file}
    vcf_file=${bam_file/.subset.markdups.proper_paired.bam/.g.vcf}
    echo $vcf_file
done



## run in a for loop
for bam_file in *.proper_paired.bam
do
    vcf_file=${bam_file/.subset.markdups.proper_paired.bam/.g.vcf}
    gatk HaplotypeCaller \
        -R ../data/PlasmoDB-55_PvivaxP01_Genome.fasta \
        -I ${bam_file} \
        -O ${vcf_file} \
        -ERC GVCF \
        --native-pair-hmm-threads 4
done

## merge the different .g.vcf files
gatk CombineGVCFs \
    -R ../data/PlasmoDB-55_PvivaxP01_Genome.fasta \
    --variant HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004.g.vcf \
    --variant HJMTYDSX2_104610-001-021_CGGATTGA-ACGTCGTT_L004.g.vcf \
    --variant HJMTYDSX2_104610-001-037_AGCTCCTA-CAACACAG_L004.g.vcf\
    --variant HJMTYDSX2_104610-001-040_GTCAGTTG-CGATTGGA_L004.g.vcf \
    --variant HJMTYDSX2_104610-001-046_TAGTTGCG-TCTGTCGT_L004.g.vcf \
    -O merged.g.vcf.gz

## do unified genotyping
gatk GenotypeGVCFs \
    -V merged.g.vcf.gz \
    -O merged.vcf.gz \
    -R ../data/PlasmoDB-55_PvivaxP01_Genome.fasta



##### SNPs #####
## select SNPs and do basic filtering on the SNPs based
## on the GATK recommended settings
gatk SelectVariants \
    -V merged.vcf.gz \
    -select-type SNP \
    -O merged.snps.vcf.gz

## add filter information
gatk VariantFiltration \
    -V merged.snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O merged.filter_added.snps.vcf.gz

## select only SNPs that passed the filter
gatk SelectVariants \
    -R ../data/PlasmoDB-55_PvivaxP01_Genome.fasta \
    -V merged.filter_added.snps.vcf.gz \
    --exclude-filtered true \
    -O merged.filtered.snps.vcf.gz




##### INDELS #####
## select indels
gatk SelectVariants \
    -V merged.vcf.gz \
    -select-type INDEL \
    -O merged.indels.vcf.gz


## add filter information
gatk VariantFiltration \
    -V merged.indels.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O merged.filter_added.indels.vcf.gz

## select only indels that passed the filter passed filter
gatk SelectVariants \
    -R ../data/PlasmoDB-55_PvivaxP01_Genome.fasta \
    -V merged.filter_added.indels.vcf.gz \
    --exclude-filtered true \
    -O merged.filtered.indels.vcf.gz


#### merge SNPs and INDELS again
## merge the SNPs and Indels again to produce one output file
java -jar /home/itg.be/pmonsieurs/software/picard/picard.jar SortVcf \
    I=merged.filtered.snps.vcf.gz \
    I=merged.filtered.indels.vcf.gz \
    O=merged.filtered.vcf.gz