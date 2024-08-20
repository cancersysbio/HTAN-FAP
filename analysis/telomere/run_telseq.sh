#!/bin/bash

## Author: Aziz

## Run TelSeq of the three cohorts. This is script for FAP multi-region cohort

## WG or WE
method="WE"

oak_path="<set oak path for azizk>"
input_path=${oak_path}/projects/FAP/results/sarek_FuchouTang_wes

bam_path=${oak_path}/projects/FAP/results/sarek_FuchouTang_wes
input_tsv=${oak_path}/projects/FAP/results/sarek_FuchouTang_wes/conf/tumor_normal_pairs.tsv

outdir_main=${input_path}/TelSeq/WES

# Slrum account
slrum_account=default
slrum_p=normal

genome='hg38'

while read tumor normal; do
  
  echo "$tumor"
    
  outdir=${outdir_main}/${tumor}
  mkdir -p ${outdir}
  cd ${outdir}

  tumor_bam=${input_path}/Preprocessing/${tumor}/Recalibrated/${tumor}.recal.bam
  normal_bam=${input_path}/Preprocessing/${normal}/Recalibrated/${normal}.recal.bam 
  
  ## RUN TelSeq Tumor and normal
  sbatch -J $tumor --export=ALL --mem 4G -c 1 --mail-user=azizk@stanford.edu --account=default --partition=interactive --mail-type=FAIL -t 0-8:00 --wrap="${oak_path}/tools/bin/telseq_v0.0.2 ${tumor_bam} -o ${outdir}/${tumor}.telseq.txt -r 150 -m"
  sbatch -J $tumor --export=ALL --mem 4G -c 1 --mail-user=azizk@stanford.edu --account=default --partition=interactive --mail-type=FAIL -t 0-8:00 --wrap="${oak_path}/tools/bin/telseq_v0.0.2 ${normal_bam} -o ${outdir}/${normal}.telseq.txt -r 150 -m"

  done <${input_tsv}

