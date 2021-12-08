## copy paste 5 full set data files (bam without
## human reads + index file) to the training directory
cd /user/antwerpen/205/vsc20587/scratch/plasmodium_moi/data/wgs_vivax2021_bam
cp HJMTYDSX2_104610-001-040_GTCAGTTG-CGATTGGA_L004_NoHuman_CoordSorted.bam* /user/antwerpen/205/vsc20587/scratch/ext_training_plasmodium/data/full/
cp HJMTYDSX2_104610-001-021_CGGATTGA-ACGTCGTT_L004_NoHuman_CoordSorted.bam* /user/antwerpen/205/vsc20587/scratch/ext_training_plasmodium/data/full/
cp HJMTYDSX2_104610-001-046_TAGTTGCG-TCTGTCGT_L004_NoHuman_CoordSorted.bam* /user/antwerpen/205/vsc20587/scratch/ext_training_plasmodium/data/full/
cp HJMTYDSX2_104610-001-001_AGTCTCAC-CGTCCATT_L004_NoHuman_CoordSorted.bam* /user/antwerpen/205/vsc20587/scratch/ext_training_plasmodium/data/full/
cp HJMTYDSX2_104610-001-037_AGCTCCTA-CAACACAG_L004_NoHuman_CoordSorted.bam* /user/antwerpen/205/vsc20587/scratch/ext_training_plasmodium/data/full/


## only extract the first chromosome from the .bam files
## using samtools
module load SAMtools
bam_out_dir=/user/antwerpen/205/vsc20587/scratch/ext_training_plasmodium/data/subset
mkdir $bam_out_dir
cd /user/antwerpen/205/vsc20587/scratch/ext_training_plasmodium/data/full/
for bam_file in *.bam
do
    bam_out_prefix=${bam_file%_NoHuman_CoordSorted.bam}
    bam_out_file=${bam_out_dir}/${bam_out_prefix}.subset.bam
    echo $bam_out_file
    samtools view $bam_file PvP01_01_v1 -b > $bam_out_file
done

## only select proper paired reads, to make sure that 
## R1 and R2 are containing the same amount of reads
for bam_file in *.subset.bam
do
    bam_file_prefix=${bam_file%.bam}
    # samtools view -bf 0x2 $bam_file > ${bam_file_prefix}.proper_paired.bam
    samtools sort ${bam_file_prefix}.proper_paired.bam -o ${bam_file_prefix}.proper_paired.sorted.bam
done

## convert the subset bam files to fastq files
fastq_dir=/user/antwerpen/205/vsc20587/scratch/ext_training_plasmodium/data/fastq
# mkdir $fastq_dir
cd ${bam_out_dir}
for bam_file in *.proper_paired.sorted.bam
do
    sample=${bam_file%.proper_paired.sorted.bam}
    echo $sample
    samtools fastq -@ 1 -1 ${fastq_dir}/${sample}_R1.fastq.gz -2 ${fastq_dir}/${sample}_R2.fastq.gz ${bam_file} -0 ${fastq_dir}/${sample}_R0.fastq.gz
done


## curation of data: now only align read to the first 
## chromosome to make sure that first and second read
## fastq file contain same amount of reads
cd /user/antwerpen/205/vsc20587/scratch/ext_training_plasmodium/data/
fastq_dir=/user/antwerpen/205/vsc20587/scratch/ext_training_plasmodium/data/fastq/
bwa_dir=/user/antwerpen/205/vsc20587/scratch/ext_training_plasmodium/data/bwa_chrom1/
mkdir $bwa_dir
cd $bwa_dir
wget https://plasmodb.org/common/downloads/Current_Release/PvivaxP01/fasta/data/PlasmoDB-55_PvivaxP01_Genome.fasta

## select only first chromosome from fasta file
module load BWA
sed '/>PvP01_01_v2/,/>PvP01_02_v2/!d;/PvP01_02_v2/q' PlasmoDB-55_PvivaxP01_Genome.fasta > PlasmoDB-55_PvivaxP01_Genome.PvP01_01.fasta
bwa index PlasmoDB-55_PvivaxP01_Genome.PvP01_01.fasta

## map fastq to the reference genome
cd $fastq_dir
for fastq_file_R1 in *_R1.fastq.gz
do
    echo $fastq_file_R1
    file_prefix=${fastq_file_R1%_R1.fastq.gz}
    fastq_file_R2=${file_prefix}_R2.fastq.gz
    # bwa_file=${bwa_dir}/${file_prefix}.bam
    bwa mem ${bwa_dir}/PlasmoDB-55_PvivaxP01_Genome.PvP01_01.fasta ${fastq_file_R1} ${fastq_file_R2} | samtools sort -o ${bwa_dir}/${file_prefix}.bam
done


## create a set of cured fastq files, which should
## result in the same amount of reads for R1 and 
## R2 files
