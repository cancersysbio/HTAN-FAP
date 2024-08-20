#!/bin/bash

## Author: Aziz

## Run TelomereHunter of the three cohorts. This is script for FAP multi-region cohort

## WG or WE
method="WE"

oak_path="<set oak path for azizk>"
input_path=${oak_path}/projects/FAP/results/sarek_FuchouTang_wes

bam_path=${oak_path}/projects/FAP/results/sarek_FuchouTang_wes
input_tsv=${oak_path}/projects/FAP/results/sarek_FuchouTang_wes/conf/tumor_normal_pairs.tsv

 outdir_main=${input_path}/TelomereHunter/WES

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
  
   ## Run TelomereHunter
   sbatch -J $tumor --export=ALL --mem 4G -c 2 --mail-user=azizk@stanford.edu --account=${slrum_account} --partition=${slrum_p} --mail-type=FAIL -t 0-8:00 \
  --wrap="${oak_path}/tools/bin/telomerehunter_v1.1.0 --pid ${tumor} --inputBamTumor ${tumor_bam} --inputBamControl ${normal_bam} --parallel --outPath ${outdir} -pff all --bandingFile /oak/stanford/groups/ccurtis2/isabl/assemblies/GRCh38/annotations/chromosome.band.hg38.bed"

done <${input_tsv}

