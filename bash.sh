!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=4
#SBATCH --time=25-60:00:00
#SBATCH --job-name=LDmap

# Load env
. /local/env/envsnakemake-6.0.5.sh
. /local/env/envsingularity-3.8.5.sh
. /local/env/envconda.sh

datadir=$(cat datadir.conf)
scratchdir=$(cat scratch.conf)
dataset=${1}
chrom=${2}
randomid=$(echo $RANDOM | md5sum | head -c 20; echo;)
ncores=4

export OMP_NUM_THREADS=$ncores

# Init pipeline
echo "Create directory ${dataset}_${chrom}_${randomid}"

echo "Build environment"
git clone https://github.com/ThomasBrazier/LDRecombinationMaps-pipeline.git $scratchdir/${dataset}_${chrom}_${randomid}
cd $scratchdir/${dataset}_${chrom}_${randomid}
singularity pull ldhot.sif docker://tombrazier/ldhot:v1.0
singularity pull ldhat.sif docker://tombrazier/ldhat:v1.0
mkdir lk_files
cd lk_files
wget https://github.com/auton1/LDhat/raw/master/lk_files/lk_n100_t0.001.gz
wget https://github.com/auton1/LDhat/raw/master/lk_files/lk_n100_t0.01.gz
gunzip *.gz
cd ..

echo "Copy data"
mkdir $scratchdir/${dataset}_${chrom}_${randomid}/data
mkdir $scratchdir/${dataset}_${chrom}_${randomid}/data/${dataset}
cp $datadir/data/$dataset/* $scratchdir/${dataset}_${chrom}_${randomid}/data/${dataset}/
cp -r $datadir/data/$dataset/structure $scratchdir/${dataset}_${chrom}_${randomid}/data/$dataset/

echo "Run pipeline"
snakemake -s Snakefile -p -j $ncores --configfile data/${dataset}/config.yaml --use-conda --use-singularity --nolock --rerun-incomplete --printshellcmds --config dataset=${dataset} chrom=${chrom} cores=$ncores

echo "Check results"
test -f $scratchdir/${dataset}_${chrom}_${randomid}/data/${dataset}/K*.pop*/ldhot/*.hot_summary.txt.gz && echo "LDhot summary exists"
test -f $scratchdir/${dataset}_${chrom}_${randomid}/data/${dataset}/K*.pop*/ldhot/*.hotspots.txt.gz && echo "LDhot hotspots exists"

echo "Clean temporary files"
bash clean.sh $dataset

echo "Sync results back"
rsync -ah $scratchdir/${dataset}_${chrom}_${randomid}/data/$dataset/ $datadir/data/$dataset/

echo "Clean scratch"
#rm -rf $scratchdir/${dataset}_${chrom}_${randomid}
