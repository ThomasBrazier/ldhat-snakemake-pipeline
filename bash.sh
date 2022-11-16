#!/bin/bash
#SBATCH --mail-user=thomas.brazier@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --mem=60GB
#SBATCH --time=20-60:00:00
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
ncores=16

# Init pipeline
echo "Create directory ${dataset}_${randomid}"

echo "Build environment"
git clone https://github.com/ThomasBrazier/LDRecombinationMaps-pipeline.git $scratchdir/${dataset}_${randomid}
cd $scratchdir/${dataset}_${randomid}
singularity pull ldhot.sif docker://tombrazier/ldhot:v1.0
singularity pull ldhat.sif docker://tombrazier/ldhat:v1.0
wget https://github.com/auton1/LDhat/tree/master/lk_files && gunzip lk_files/*.gz

echo "Copy data"
mkdir $scratchdir/${dataset}_${randomid}/data
cp -v $datadir/data/$dataset $scratchdir/${dataset}_${randomid}/data/
cp -v $datadir/data/$dataset/structure $scratchdir/${dataset}_${randomid}/data/$dataset/

echo "Run pipeline"
snakemake -s Snakefile -p -j $ncores --configfile data/${dataset}/config.yaml --use-conda --use-singularity --nolock --rerun-incomplete --config dataset=${dataset} chrom=${chrom}

echo "Check results"
if [[ -f $scratchdir/${dataset}_${randomid}/data/*/ldhot/*.hot_summary.txt.gz && -f $scratchdir/${dataset}_${randomid}/data/*/ldhot/*.hotspots.txt.gz ]];
then
    echo "Result files exist."
else
    exit 1
fi

echo "Clean temporary files"
bash clean.sh $dataset

echo "Sync results back"
rsync -avh $scratchdir/${dataset}_${randomid}/data/$dataset/ $datadir/data/$dataset/

echo "Clean scratch"
rm -r $scratchdir/${dataset}_${randomid}
