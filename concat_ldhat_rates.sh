#! /bin/bash
# Concatenate all MCMC chains of rates outputs

wdirpop=${1}
dataset=${2}
chromosome=${3}
bpen=${4}

cd $wdirpop/ldhat/${dataset}.${chromosome}/

test -f bpen${bpen}.rates_noheader.txt && rm bpen${bpen}.rates_noheader.txt

n_batch=$(ls | grep "bpen${bpen}.batch_" | grep ".rates.txt" | wc -l)
echo $n_batch

for f in $(seq 1 $(( $n_batch )))
do
  echo "Processing bpen${bpen}.batch_${f}.rates.txt.gz file..."
  zcat bpen${bpen}.batch_${f}.rates.txt | tail -n +2 >> bpen${bpen}.rates_noheader.txt
done
n_snps=$(cat bpen${bpen}.rates_noheader.txt | wc -l)
echo "$n_snps	$n_snps" > bpen${bpen}.rates.txt
cat bpen${bpen}.rates_noheader.txt >> bpen${bpen}.rates.txt
mv bpen${bpen}.rates.txt $wdirpop/ldhat/${dataset}.${chromosome}.bpen${bpen}.rates.txt
