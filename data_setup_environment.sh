## copy paste the fastq file from the calcua
## to the test server
cd ~/programming/ext_training_plasmodium/data/
scp calcua:/user/antwerpen/205/vsc20587/scratch/ext_training_plasmodium/data/fastq/* ./fastq/
scp -rp fastq/ training:/home/itg.be/pmonsieurs/data/

## to be run on the training server
sudo adduser trainee01 ## with same password as user name


## install the necessary software
mkdir /home/itg.be/pmonsieurs/software
cd /home/itg.be/pmonsieurs/software/

sudo apt-get make
sudo apt-get gcc
sudo apt-get install libz-dev

## install bwa
git clone https://github.com/lh3/bwa
cd bwa; make

## download Picard
cd /home/itg.be/pmonsieurs/software
mkdir picard; cd picard
sudo apt-get install default-jre
wget https://github.com/broadinstitute/picard/releases/download/2.26.6/picard.jar

## install htop
sudo apt-get install htop


## install GATK
wget https://github.com/broadinstitute/gatk/releases/download/4.2.3.0/gatk-4.2.3.0.zip
unzip gatk-4.2.3.0.zip
## create symbolic link to python3 as python is not defined.
## also create symbolic link to gatk in the /usr/local/bin 
## directory
sudo ln -s /usr/bin/python3 /usr/bin/python
sudo ln -s $PWD/gatk /usr/local/bin/gatk



