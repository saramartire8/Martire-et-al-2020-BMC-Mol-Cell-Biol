#!/bin/bash
#SBATCH --job-name=profile2                              # job name
#SBATCH --partition=super                                 # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                         # number of nodes requested by user
#SBATCH --time=0-24:00:00                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=job_profile.%j.out                         # standard output file name
#SBATCH --error=job_profile.%j.err                         # standard error output file name
#SBATCH --mail-user=sara.martire@utsouthwestern.edu           # specify an email address
#SBATCH --mail-type=ALL                                   # send email when job status change (start, end, abortion and etc.)

module load deeptools

computeMatrix reference-point --referencePoint center -R H3K27ac_going_DOWN.bed -b 3000 -a 3000 --sortRegions no -S /project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/bigwig_normalized/129_p300_V_rep3_norm.nodup.bw /project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/bigwig_normalized/282_p300_V_rep3_norm.nodup.bw --skipZeros -o matrix_p300V-down.gz --outFileSortedRegions matrix_p300V-down.bed --numberOfProcessors max
plotProfile -m matrix_p300V-down.gz -out 129vs282_p300V-down.pdf --perGroup --color blue red --samplesLabel p300V_129 p300V_282 --startLabel center --plotTitle "p300V under down"
