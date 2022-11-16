#! /bin/bash
# Clean data repository from temporary files
set=$1

rm data/$set/vcf_newheader.vcf.gz

rm data/$set/K*.pop*/${set}.phased.chromosome.*
rm data/$set/K*.pop*/${set}.pop.vcf*
rm data/$set/K*.pop*/${set}.chromosome.1.vcf*
rm data/$set/K*.pop*/ldhat/${set}.*.new_lk.txt
rm data/$set/K*.pop*/ldhat/${set}.*.locs
rm data/$set/K*.pop*/ldhat/${set}.*.sites
rm data/$set/K*.pop*/ldhot/${set}.*.summary.log
