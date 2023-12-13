#! /bin/bash
# Concatenate all MCMC chains of rates outputs

wdirpop=${1}
dataset=${2}
chromosome=${3}
bpen=${4}
bigchunk=${5}
smalloverlap=${6}
chunk={7}

cd $wdirpop/ldhat/${dataset}.${chromosome}/

test -f bpen${bpen}.rates_noheader.txt && rm bpen${bpen}.rates_noheader.txt

n_batch=$(ls | grep "bpen${bpen}.batch_" | grep ".rates.txt" | wc -l)
echo $n_batch

touch bpen${bpen}.rates_noheader.txt

# Remove overlapping SNPs
# SNP positions are in columns
# Cut extra columns
echo "First chunk"
zcat bpen${bpen}.batch_1.rates.txt.gz | tail -n +2 | cut -f $(seq -s , 1 $bigchunk) > bpen${bpen}.rates_noheader.txt

for f in $(seq 2 $(( $nbatch-1 )))
do
  echo "Processing bpen${bpen}.batch_${f}.rates.txt.gz file..."
  zcat bpen${bpen}.batch_${f}.rates.txt.gz | tail -n +2 | cut -f $(seq -s , $(( $smalloverlap+1 )) $bigchunk) > tmp.txt
  paste bpen${bpen}.rates_noheader.txt tmp.txt > tmp2
  cat tmp2 > bpen${bpen}.rates_noheader.txt
done

echo "Last chunk"
zcat bpen${bpen}.batch_${nbatch}.rates.txt.gz | tail -n +2 | cut -f $(seq -s , $(( $smalloverlap+1 )) $chunk) > tmp.txt
paste bpen${bpen}.rates_noheader.txt tmp.txt > tmp2
cat tmp2 > bpen${bpen}.rates_noheader.txt

# Add header
n_snps=$(awk '{print NF}' bpen${bpen}.rates_noheader.txt | sort -nu | tail -n 1)

rm tmp2
rm tmp.txt

echo "$n_snps	$n_snps" > bpen${bpen}.rates.txt
cat bpen${bpen}.rates_noheader.txt >> bpen${bpen}.rates.txt
