# -*- coding: utf-8 -*-

import sys
import os
import csv
import gzip
import copy
import math
import numpy as np
from scipy.interpolate import CubicSpline
import random
import json

# 形成完整的imputation过程指令文件 HGSVC和cuteSVTrio相同样本之间的比较
def make_imputation_command_file(file_base) :
    chr_ls = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]
    for chr in chr_ls :
        cmd_lines = []
        cmd_lines.append("#!/bin/bash".replace("chr1",chr))
        cmd_lines.append("".replace("chr1",chr))
        cmd_lines.append(". ~/.bashrc".replace("chr1",chr))
        cmd_lines.append("conda activate lx".replace("chr1",chr))
        # 准备工作
        cmd_lines.append("python pick_chr_from_vcf.py chr1".replace("chr1",chr))
        cmd_lines.append("bgzip -f output_real/output_trio/HiFI/hg38/M1_phased/panel/merge.truvari.rectified.chr1.vcf ; tabix -f output_real/output_trio/HiFI/hg38/M1_phased/panel/merge.truvari.rectified.chr1.vcf.gz".replace("chr1",chr))
        cmd_lines.append("bgzip -f answer_vcf/freeze4/merge.truvari.rectified.chr1.vcf ; tabix -f answer_vcf/freeze4/merge.truvari.rectified.chr1.vcf.gz".replace("chr1",chr))
        cmd_lines.append("python pick_samples_from_vcf.py chr1".replace("chr1",chr))
        cmd_lines.append("bgzip -f 1kGP/1kGP_chr1_SNP_INDEL.260.target.vcf ; tabix -f 1kGP/1kGP_chr1_SNP_INDEL.260.target.vcf.gz".replace("chr1",chr))
        cmd_lines.append("python split_vcf_file.py chr1".replace("chr1",chr))
        cmd_lines.append("bgzip -f 1kGP/1kGP_chr1_SV.260.target.vcf ; tabix -f 1kGP/1kGP_chr1_SV.260.target.vcf.gz ; bcftools index 1kGP/1kGP_chr1_SV.260.target.vcf.gz".replace("chr1",chr))
        cmd_lines.append("bgzip -f 1kGP/1kGP_chr1_SV.260.answer.vcf ; tabix -f 1kGP/1kGP_chr1_SV.260.answer.vcf.gz ; bcftools index 1kGP/1kGP_chr1_SV.260.answer.vcf.gz".replace("chr1",chr))
        cmd_lines.append("bcftools concat 1kGP/1kGP_chr1_SNP_INDEL.260.target.vcf.gz 1kGP/1kGP_chr1_SV.260.target.vcf.gz -o 1kGP/1kGP_chr1_all.260.unsorted.vcf".replace("chr1",chr))
        cmd_lines.append("bgzip -f 1kGP/1kGP_chr1_all.260.unsorted.vcf".replace("chr1",chr))
        cmd_lines.append("bcftools sort 1kGP/1kGP_chr1_all.260.unsorted.vcf.gz -o 1kGP/1kGP_chr1_all.260.target.vcf".replace("chr1",chr))
        cmd_lines.append("bgzip -f 1kGP/1kGP_chr1_all.260.target.vcf ; tabix -f 1kGP/1kGP_chr1_all.260.target.vcf.gz".replace("chr1",chr))
        
        # 处理16+15个samples
        cmd_lines.append("conda activate lx".replace("chr1",chr))
        cmd_lines.append("bgzip -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.chr1.overlap.sv.vcf ; tabix -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.chr1.overlap.sv.vcf.gz".replace("chr1",chr))
        cmd_lines.append("bgzip -f answer_vcf/freeze4/samples_31/HGSVC.chr1.overlap.sv.vcf ; tabix -f answer_vcf/freeze4/samples_31/HGSVC.chr1.overlap.sv.vcf.gz".replace("chr1",chr))
        cmd_lines.append("bgzip -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/1kGP.chr1.overlap.SNP_INDEL.vcf ; tabix -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/1kGP.chr1.overlap.SNP_INDEL.vcf.gz".replace("chr1",chr))
        cmd_lines.append("bgzip -f answer_vcf/freeze4/samples_31/1kGP.chr1.overlap.SNP_INDEL.vcf ; tabix -f answer_vcf/freeze4/samples_31/1kGP.chr1.overlap.SNP_INDEL.vcf.gz".replace("chr1",chr))
        cmd_lines.append("bcftools concat output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/1kGP.chr1.overlap.SNP_INDEL.vcf.gz output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.chr1.overlap.sv.vcf.gz -o output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.1kGP.chr1.unsorted.vcf".replace("chr1",chr))
        cmd_lines.append("bgzip -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.1kGP.chr1.unsorted.vcf ; tabix -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.1kGP.chr1.unsorted.vcf.gz".replace("chr1",chr))
        cmd_lines.append("bcftools sort output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.1kGP.chr1.unsorted.vcf.gz -o output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.1kGP.chr1.panel.vcf".replace("chr1",chr))
        cmd_lines.append("bgzip -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.1kGP.chr1.panel.vcf ; tabix -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.1kGP.chr1.panel.vcf.gz".replace("chr1",chr))
        cmd_lines.append("bcftools concat answer_vcf/freeze4/samples_31/1kGP.chr1.overlap.SNP_INDEL.vcf.gz answer_vcf/freeze4/samples_31/HGSVC.chr1.overlap.sv.vcf.gz -o answer_vcf/freeze4/samples_31/HGSVC.1kGP.chr1.unsorted.vcf".replace("chr1",chr))
        cmd_lines.append("bgzip -f answer_vcf/freeze4/samples_31/HGSVC.1kGP.chr1.unsorted.vcf ; tabix -f answer_vcf/freeze4/samples_31/HGSVC.1kGP.chr1.unsorted.vcf.gz".replace("chr1",chr))
        cmd_lines.append("bcftools sort answer_vcf/freeze4/samples_31/HGSVC.1kGP.chr1.unsorted.vcf.gz -o answer_vcf/freeze4/samples_31/HGSVC.1kGP.chr1.panel.vcf".replace("chr1",chr))
        cmd_lines.append("bgzip -f answer_vcf/freeze4/samples_31/HGSVC.1kGP.chr1.panel.vcf ; tabix -f answer_vcf/freeze4/samples_31/HGSVC.1kGP.chr1.panel.vcf.gz".replace("chr1",chr))
        cmd_lines.append("minimac4 --compress-reference output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.1kGP.chr1.panel.vcf.gz > output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.1kGP.chr1.panel.msav -t 32".replace("chr1",chr))
        cmd_lines.append("minimac4 --compress-reference answer_vcf/freeze4/samples_31/HGSVC.1kGP.chr1.panel.vcf.gz > answer_vcf/freeze4/samples_31/HGSVC.1kGP.chr1.panel.msav -t 32".replace("chr1",chr))
        cmd_lines.append("minimac4 -f GT,HDS,GP -t 8 -c 15000000 output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.1kGP.chr1.panel.msav 1kGP/1kGP_chr1_all.260.target.vcf.gz -o output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.1kGP.chr1.minimac4.vcf.gz".replace("chr1",chr))
        cmd_lines.append("minimac4 -f GT,HDS,GP -t 8 -c 15000000 answer_vcf/freeze4/samples_31/HGSVC.1kGP.chr1.panel.msav 1kGP/1kGP_chr1_all.260.target.vcf.gz -o answer_vcf/freeze4/samples_31/HGSVC.1kGP.chr1.minimac4.vcf.gz".replace("chr1",chr))
        cmd_lines.append("python pick_sv_from_imputed_vcf.py output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.1kGP.chr1.minimac4.vcf.gz output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.1kGP.chr1.minimac4.SV.vcf".replace("chr1",chr))
        cmd_lines.append("bgzip -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.1kGP.chr1.minimac4.SV.vcf ; tabix -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.1kGP.chr1.minimac4.SV.vcf.gz".replace("chr1",chr))
        cmd_lines.append("python pick_sv_from_imputed_vcf.py answer_vcf/freeze4/samples_31/HGSVC.1kGP.chr1.minimac4.vcf.gz answer_vcf/freeze4/samples_31/HGSVC.1kGP.chr1.minimac4.SV.vcf".replace("chr1",chr))
        cmd_lines.append("bgzip -f answer_vcf/freeze4/samples_31/HGSVC.1kGP.chr1.minimac4.SV.vcf ; tabix -f answer_vcf/freeze4/samples_31/HGSVC.1kGP.chr1.minimac4.SV.vcf.gz".replace("chr1",chr))
        cmd_lines.append("conda activate truvari".replace("chr1",chr))
        cmd_lines.append("rm -r output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.minimac4.chr1.r2; truvari bench -b 1kGP/1kGP_chr1_SV.260.answer.vcf.gz -c output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.1kGP.chr1.minimac4.SV.vcf.gz -o output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.minimac4.chr1.r2 -p 0 -r 5000 -C 5000 --passonly".replace("chr1",chr))
        cmd_lines.append("rm -r answer_vcf/freeze4/samples_31/HGSVC.minimac4.chr1.r2; truvari bench -b 1kGP/1kGP_chr1_SV.260.answer.vcf.gz -c answer_vcf/freeze4/samples_31/HGSVC.1kGP.chr1.minimac4.SV.vcf.gz -o answer_vcf/freeze4/samples_31/HGSVC.minimac4.chr1.r2 -p 0 -r 5000 -C 5000 --passonly".replace("chr1",chr))
        cmd_lines.append("rm -r output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.minimac4.chr1.easy.r2; truvari bench -b 1kGP/1kGP_chr1_SV.260.answer.vcf.gz -c output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.1kGP.chr1.minimac4.SV.vcf.gz --includebed answer_vcf/genome_stratifications/GRCh38_notinalldifficultregions.chr1.bed -o output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.minimac4.chr1.easy.r2 -p 0 -r 5000 -C 5000 --passonly".replace("chr1",chr))
        cmd_lines.append("rm -r answer_vcf/freeze4/samples_31/HGSVC.minimac4.chr1.easy.r2; truvari bench -b 1kGP/1kGP_chr1_SV.260.answer.vcf.gz -c answer_vcf/freeze4/samples_31/HGSVC.1kGP.chr1.minimac4.SV.vcf.gz --includebed answer_vcf/genome_stratifications/GRCh38_notinalldifficultregions.chr1.bed -o answer_vcf/freeze4/samples_31/HGSVC.minimac4.chr1.easy.r2 -p 0 -r 5000 -C 5000 --passonly".replace("chr1",chr))
        cmd_lines.append("rm -r output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.minimac4.chr1.diff.r2; truvari bench -b 1kGP/1kGP_chr1_SV.260.answer.vcf.gz -c output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.1kGP.chr1.minimac4.SV.vcf.gz --includebed answer_vcf/genome_stratifications/GRCh38_alldifficultregions.chr1.bed -o output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_31/cuteSVTrio.minimac4.chr1.diff.r2 -p 0 -r 5000 -C 5000 --passonly".replace("chr1",chr))
        cmd_lines.append("rm -r answer_vcf/freeze4/samples_31/HGSVC.minimac4.chr1.diff.r2; truvari bench -b 1kGP/1kGP_chr1_SV.260.answer.vcf.gz -c answer_vcf/freeze4/samples_31/HGSVC.1kGP.chr1.minimac4.SV.vcf.gz --includebed answer_vcf/genome_stratifications/GRCh38_alldifficultregions.chr1.bed -o answer_vcf/freeze4/samples_31/HGSVC.minimac4.chr1.diff.r2 -p 0 -r 5000 -C 5000 --passonly".replace("chr1",chr))
        
        # 处理16个samples
        cmd_lines.append("conda activate lx".replace("chr1",chr))
        cmd_lines.append("bgzip -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.chr1.overlap.sv.vcf ; tabix -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.chr1.overlap.sv.vcf.gz".replace("chr1",chr))
        cmd_lines.append("bgzip -f answer_vcf/freeze4/samples_16/HGSVC.chr1.overlap.sv.vcf ; tabix -f answer_vcf/freeze4/samples_16/HGSVC.chr1.overlap.sv.vcf.gz".replace("chr1",chr))
        cmd_lines.append("bgzip -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/1kGP.chr1.overlap.SNP_INDEL.vcf ; tabix -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/1kGP.chr1.overlap.SNP_INDEL.vcf.gz".replace("chr1",chr))
        cmd_lines.append("bgzip -f answer_vcf/freeze4/samples_16/1kGP.chr1.overlap.SNP_INDEL.vcf ; tabix -f answer_vcf/freeze4/samples_16/1kGP.chr1.overlap.SNP_INDEL.vcf.gz".replace("chr1",chr))
        cmd_lines.append("bcftools concat output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/1kGP.chr1.overlap.SNP_INDEL.vcf.gz output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.chr1.overlap.sv.vcf.gz -o output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.1kGP.chr1.unsorted.vcf".replace("chr1",chr))
        cmd_lines.append("bgzip -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.1kGP.chr1.unsorted.vcf ; tabix -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.1kGP.chr1.unsorted.vcf.gz".replace("chr1",chr))
        cmd_lines.append("bcftools sort output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.1kGP.chr1.unsorted.vcf.gz -o output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.1kGP.chr1.panel.vcf".replace("chr1",chr))
        cmd_lines.append("bgzip -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.1kGP.chr1.panel.vcf ; tabix -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.1kGP.chr1.panel.vcf.gz".replace("chr1",chr))
        cmd_lines.append("bcftools concat answer_vcf/freeze4/samples_16/1kGP.chr1.overlap.SNP_INDEL.vcf.gz answer_vcf/freeze4/samples_16/HGSVC.chr1.overlap.sv.vcf.gz -o answer_vcf/freeze4/samples_16/HGSVC.1kGP.chr1.unsorted.vcf".replace("chr1",chr))
        cmd_lines.append("bgzip -f answer_vcf/freeze4/samples_16/HGSVC.1kGP.chr1.unsorted.vcf ; tabix -f answer_vcf/freeze4/samples_16/HGSVC.1kGP.chr1.unsorted.vcf.gz".replace("chr1",chr))
        cmd_lines.append("bcftools sort answer_vcf/freeze4/samples_16/HGSVC.1kGP.chr1.unsorted.vcf.gz -o answer_vcf/freeze4/samples_16/HGSVC.1kGP.chr1.panel.vcf".replace("chr1",chr))
        cmd_lines.append("bgzip -f answer_vcf/freeze4/samples_16/HGSVC.1kGP.chr1.panel.vcf ; tabix -f answer_vcf/freeze4/samples_16/HGSVC.1kGP.chr1.panel.vcf.gz".replace("chr1",chr))
        cmd_lines.append("minimac4 --compress-reference output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.1kGP.chr1.panel.vcf.gz > output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.1kGP.chr1.panel.msav -t 32".replace("chr1",chr))
        cmd_lines.append("minimac4 --compress-reference answer_vcf/freeze4/samples_16/HGSVC.1kGP.chr1.panel.vcf.gz > answer_vcf/freeze4/samples_16/HGSVC.1kGP.chr1.panel.msav -t 32".replace("chr1",chr))
        cmd_lines.append("minimac4 -f GT,HDS,GP -t 8 -c 15000000 output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.1kGP.chr1.panel.msav 1kGP/1kGP_chr1_all.260.target.vcf.gz -o output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.1kGP.chr1.minimac4.vcf.gz".replace("chr1",chr))
        cmd_lines.append("minimac4 -f GT,HDS,GP -t 8 -c 15000000 answer_vcf/freeze4/samples_16/HGSVC.1kGP.chr1.panel.msav 1kGP/1kGP_chr1_all.260.target.vcf.gz -o answer_vcf/freeze4/samples_16/HGSVC.1kGP.chr1.minimac4.vcf.gz".replace("chr1",chr))
        cmd_lines.append("python pick_sv_from_imputed_vcf.py output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.1kGP.chr1.minimac4.vcf.gz output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.1kGP.chr1.minimac4.SV.vcf".replace("chr1",chr))
        cmd_lines.append("bgzip -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.1kGP.chr1.minimac4.SV.vcf ; tabix -f output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.1kGP.chr1.minimac4.SV.vcf.gz".replace("chr1",chr))
        cmd_lines.append("python pick_sv_from_imputed_vcf.py answer_vcf/freeze4/samples_16/HGSVC.1kGP.chr1.minimac4.vcf.gz answer_vcf/freeze4/samples_16/HGSVC.1kGP.chr1.minimac4.SV.vcf".replace("chr1",chr))
        cmd_lines.append("bgzip -f answer_vcf/freeze4/samples_16/HGSVC.1kGP.chr1.minimac4.SV.vcf ; tabix -f answer_vcf/freeze4/samples_16/HGSVC.1kGP.chr1.minimac4.SV.vcf.gz".replace("chr1",chr))
        cmd_lines.append("conda activate truvari".replace("chr1",chr))
        cmd_lines.append("rm -r output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.minimac4.chr1.r2; truvari bench -b 1kGP/1kGP_chr1_SV.260.answer.vcf.gz -c output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.1kGP.chr1.minimac4.SV.vcf.gz -o output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.minimac4.chr1.r2 -p 0 -r 5000 -C 5000 --passonly".replace("chr1",chr))
        cmd_lines.append("rm -r answer_vcf/freeze4/samples_16/HGSVC.minimac4.chr1.r2; truvari bench -b 1kGP/1kGP_chr1_SV.260.answer.vcf.gz -c answer_vcf/freeze4/samples_16/HGSVC.1kGP.chr1.minimac4.SV.vcf.gz -o answer_vcf/freeze4/samples_16/HGSVC.minimac4.chr1.r2 -p 0 -r 5000 -C 5000 --passonly".replace("chr1",chr))
        cmd_lines.append("rm -r output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.minimac4.chr1.easy.r2; truvari bench -b 1kGP/1kGP_chr1_SV.260.answer.vcf.gz -c output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.1kGP.chr1.minimac4.SV.vcf.gz --includebed answer_vcf/genome_stratifications/GRCh38_notinalldifficultregions.chr1.bed -o output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.minimac4.chr1.easy.r2 -p 0 -r 5000 -C 5000 --passonly".replace("chr1",chr))
        cmd_lines.append("rm -r answer_vcf/freeze4/samples_16/HGSVC.minimac4.chr1.easy.r2; truvari bench -b 1kGP/1kGP_chr1_SV.260.answer.vcf.gz -c answer_vcf/freeze4/samples_16/HGSVC.1kGP.chr1.minimac4.SV.vcf.gz --includebed answer_vcf/genome_stratifications/GRCh38_notinalldifficultregions.chr1.bed -o answer_vcf/freeze4/samples_16/HGSVC.minimac4.chr1.easy.r2 -p 0 -r 5000 -C 5000 --passonly".replace("chr1",chr))
        cmd_lines.append("rm -r output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.minimac4.chr1.diff.r2; truvari bench -b 1kGP/1kGP_chr1_SV.260.answer.vcf.gz -c output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.1kGP.chr1.minimac4.SV.vcf.gz --includebed answer_vcf/genome_stratifications/GRCh38_alldifficultregions.chr1.bed -o output_real/output_trio/HiFI/hg38/M1_phased/panel/samples_16/cuteSVTrio.minimac4.chr1.diff.r2 -p 0 -r 5000 -C 5000 --passonly".replace("chr1",chr))
        cmd_lines.append("rm -r answer_vcf/freeze4/samples_16/HGSVC.minimac4.chr1.diff.r2; truvari bench -b 1kGP/1kGP_chr1_SV.260.answer.vcf.gz -c answer_vcf/freeze4/samples_16/HGSVC.1kGP.chr1.minimac4.SV.vcf.gz --includebed answer_vcf/genome_stratifications/GRCh38_alldifficultregions.chr1.bed -o answer_vcf/freeze4/samples_16/HGSVC.minimac4.chr1.diff.r2 -p 0 -r 5000 -C 5000 --passonly".replace("chr1",chr))
        
        with open(file_base+"/imputation_"+chr+".sh",'w') as f:
            for line in cmd_lines :
                f.write(line+"\n")
        f.close()

# 计算cuteSV-Trio和HGSVC的minimac4的imputation输出的R2，分为INS和DEL
def cal_cuteSVTrio_HGSVC_r2_svtype(cuteSVTrio_base,HGSVC_base,write_file) :
    af_range_ls = [0,0.01,0.02,0.05,0.1,0.2,0.35,0.5,0.65,0.8,1]
    INS_cuteSVTrio_af_res = [[] for x in range(10)]
    INS_HGSVC_af_res = [[] for x in range(10)]
    DEL_cuteSVTrio_af_res = [[] for x in range(10)]
    DEL_HGSVC_af_res = [[] for x in range(10)]
    cuteSVTrio_r2 = []
    HGSVC_r2 = []
    chr_ls = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]
    #chr_ls = ["chr1"]
    for chr in chr_ls :
        cuteSVTrio_file_base = cuteSVTrio_base.replace("chr1",chr) + "/tp-base.vcf.gz"
        INS_af_ls = []
        DEL_af_ls = []
        fh = gzip.open(cuteSVTrio_file_base, 'rt') if cuteSVTrio_file_base.endswith('.gz') else open(cuteSVTrio_file_base, 'rt')
        base_line_ls = []
        for line in fh:
            if line == "" :
                continue
            line = line.strip()
            if line.startswith('#'):
                line_split_ls = line.split("\t")
                if line.startswith('#CHROM'):
                    INS_base_gt_ls = [[] for x in range(len(line_split_ls)-9)]
                    DEL_base_gt_ls = [[] for x in range(len(line_split_ls)-9)]
                continue
            else :
                line_split_ls = line.split("\t")
                base_line_ls.append(line_split_ls)
        fh.close()
                
        comp_line_ls = []
        cuteSVTrio_file_comp = cuteSVTrio_base.replace("chr1",chr) + "/tp-comp.vcf.gz"
        fh = gzip.open(cuteSVTrio_file_comp, 'rt') if cuteSVTrio_file_comp.endswith('.gz') else open(cuteSVTrio_file_comp, 'rt')
        for line in fh:
            if line == "" :
                continue
            line = line.strip()
            if line.startswith('#'):
                line_split_ls = line.split("\t")
                if line.startswith('#CHROM'):
                    INS_comp_gt_ls = [[] for x in range(len(line_split_ls)-9)]
                    DEL_comp_gt_ls = [[] for x in range(len(line_split_ls)-9)]
                continue
            else :
                line_split_ls = line.split("\t")
                comp_line_ls.append(line_split_ls)
        fh.close()

        for i in range(len(base_line_ls)) :
            base_split_ls = base_line_ls[i]
            if "INS" in base_split_ls[7] :
                INS_af_ls.append(float(base_split_ls[7].split("AF=")[1].split(";")[0]))
                for j in range(9,len(base_split_ls)) :
                    INS_base_gt_ls[j-9].append(int(base_split_ls[j][0])+int(base_split_ls[j][2]))
                comp_split_ls = comp_line_ls[i]
                for j in range(9,len(comp_split_ls)) :
                    INS_comp_gt_ls[j-9].append(int(comp_split_ls[j][0])+int(comp_split_ls[j][2]))
            if "DEL" in base_split_ls[7] :
                DEL_af_ls.append(float(base_split_ls[7].split("AF=")[1].split(";")[0]))
                for j in range(9,len(base_split_ls)) :
                    DEL_base_gt_ls[j-9].append(int(base_split_ls[j][0])+int(base_split_ls[j][2]))
                comp_split_ls = comp_line_ls[i]
                for j in range(9,len(comp_split_ls)) :
                    DEL_comp_gt_ls[j-9].append(int(comp_split_ls[j][0])+int(comp_split_ls[j][2]))
                
        INS_base_gt_ls = np.array(INS_base_gt_ls)
        INS_comp_gt_ls = np.array(INS_comp_gt_ls)
        corr_matrix = np.corrcoef(INS_comp_gt_ls, INS_base_gt_ls, rowvar=False)
        n_loci = INS_base_gt_ls.shape[1]
        direct_corr = np.diag(corr_matrix, k=n_loci)  # 提取下三角中偏移 n_loci 的对角线
        for i in range(len(INS_af_ls)) :
            for j in range(1,len(af_range_ls)) :
                if str(direct_corr[i]) == "nan" or direct_corr[i] < 0 :
                    continue
                if INS_af_ls[i] == 0 :
                    INS_cuteSVTrio_af_res[0].append(direct_corr[i])
                    cuteSVTrio_r2.append(direct_corr[i])
                    break
                elif INS_af_ls[i] > af_range_ls[j-1] and INS_af_ls[i] <= af_range_ls[j] :
                    INS_cuteSVTrio_af_res[j-1].append(direct_corr[i])
                    cuteSVTrio_r2.append(direct_corr[i])
                    break
                
        DEL_base_gt_ls = np.array(DEL_base_gt_ls)
        DEL_comp_gt_ls = np.array(DEL_comp_gt_ls)
        corr_matrix = np.corrcoef(DEL_comp_gt_ls, DEL_base_gt_ls, rowvar=False)
        n_loci = DEL_base_gt_ls.shape[1]
        direct_corr = np.diag(corr_matrix, k=n_loci)  # 提取下三角中偏移 n_loci 的对角线
        for i in range(len(DEL_af_ls)) :
            for j in range(1,len(af_range_ls)) :
                if str(direct_corr[i]) == "nan" or direct_corr[i] < 0 :
                    continue
                if DEL_af_ls[i] == 0 :
                    DEL_cuteSVTrio_af_res[0].append(direct_corr[i])
                    cuteSVTrio_r2.append(direct_corr[i])
                    break
                elif DEL_af_ls[i] > af_range_ls[j-1] and DEL_af_ls[i] <= af_range_ls[j] :
                    DEL_cuteSVTrio_af_res[j-1].append(direct_corr[i])
                    cuteSVTrio_r2.append(direct_corr[i])
                    break    
    
        cuteSVTrio_file_fn = cuteSVTrio_base.replace("chr1",chr) + "/fn.vcf.gz"
        fh = gzip.open(cuteSVTrio_file_fn, 'rt') if cuteSVTrio_file_fn.endswith('.gz') else open(cuteSVTrio_file_fn, 'rt')
        for line in fh:
            if line == "" :
                continue
            line = line.strip()
            if line.startswith('#'):
                continue
            else :
                line_split_ls = line.split("\t")
                if "INS" in line_split_ls[7] :
                    af = float(line_split_ls[7].split("AF=")[1].split(";")[0])
                    for j in range(1,len(af_range_ls)) :
                        if af == 0 :
                            INS_cuteSVTrio_af_res[0].append(0)
                            cuteSVTrio_r2.append(0)
                            break
                        elif af > af_range_ls[j-1] and af <= af_range_ls[j] :
                            INS_cuteSVTrio_af_res[j-1].append(0)
                            cuteSVTrio_r2.append(0)
                            break
                if "DEL" in line_split_ls[7] :
                    af = float(line_split_ls[7].split("AF=")[1].split(";")[0])
                    for j in range(1,len(af_range_ls)) :
                        if af == 0 :
                            DEL_cuteSVTrio_af_res[0].append(0)
                            cuteSVTrio_r2.append(0)
                            break
                        elif af > af_range_ls[j-1] and af <= af_range_ls[j] :
                            DEL_cuteSVTrio_af_res[j-1].append(0)
                            cuteSVTrio_r2.append(0)
                            break
        fh.close()

        HGSVC_file_base = HGSVC_base.replace("chr1",chr) + "/tp-base.vcf.gz"
        INS_af_ls = []
        DEL_af_ls = []
        fh = gzip.open(HGSVC_file_base, 'rt') if HGSVC_file_base.endswith('.gz') else open(HGSVC_file_base, 'rt')
        base_line_ls = []
        for line in fh:
            if line == "" :
                continue
            line = line.strip()
            if line.startswith('#'):
                line_split_ls = line.split("\t")
                if line.startswith('#CHROM'):
                    INS_base_gt_ls = [[] for x in range(len(line_split_ls)-9)]
                    DEL_base_gt_ls = [[] for x in range(len(line_split_ls)-9)]
                continue
            else :
                line_split_ls = line.split("\t")
                base_line_ls.append(line_split_ls)
        fh.close()
                
        comp_line_ls = []
        HGSVC_file_comp = HGSVC_base.replace("chr1",chr) + "/tp-comp.vcf.gz"
        fh = gzip.open(HGSVC_file_comp, 'rt') if HGSVC_file_comp.endswith('.gz') else open(HGSVC_file_comp, 'rt')
        for line in fh:
            if line == "" :
                continue
            line = line.strip()
            if line.startswith('#'):
                line_split_ls = line.split("\t")
                if line.startswith('#CHROM'):
                    INS_comp_gt_ls = [[] for x in range(len(line_split_ls)-9)]
                    DEL_comp_gt_ls = [[] for x in range(len(line_split_ls)-9)]
                continue
            else :
                line_split_ls = line.split("\t")
                comp_line_ls.append(line_split_ls)
        fh.close()

        #print(len(base_line_ls))
        #print(len(comp_line_ls))

        for i in range(len(base_line_ls)) :
            base_split_ls = base_line_ls[i]
            if "INS" in base_split_ls[7] :
                INS_af_ls.append(float(base_split_ls[7].split("AF=")[1].split(";")[0]))
                for j in range(9,len(base_split_ls)) :
                    INS_base_gt_ls[j-9].append(int(base_split_ls[j][0])+int(base_split_ls[j][2]))
                comp_split_ls = comp_line_ls[i]
                for j in range(9,len(comp_split_ls)) :
                    INS_comp_gt_ls[j-9].append(int(comp_split_ls[j][0])+int(comp_split_ls[j][2]))
            if "DEL" in base_split_ls[7] :
                DEL_af_ls.append(float(base_split_ls[7].split("AF=")[1].split(";")[0]))
                for j in range(9,len(base_split_ls)) :
                    DEL_base_gt_ls[j-9].append(int(base_split_ls[j][0])+int(base_split_ls[j][2]))
                comp_split_ls = comp_line_ls[i]
                for j in range(9,len(comp_split_ls)) :
                    DEL_comp_gt_ls[j-9].append(int(comp_split_ls[j][0])+int(comp_split_ls[j][2]))
                
        INS_base_gt_ls = np.array(INS_base_gt_ls)
        INS_comp_gt_ls = np.array(INS_comp_gt_ls)
        #print(INS_base_gt_ls)
        #print(INS_comp_gt_ls)
        corr_matrix = np.corrcoef(INS_comp_gt_ls, INS_base_gt_ls, rowvar=False)
        n_loci = INS_base_gt_ls.shape[1]
        direct_corr = np.diag(corr_matrix, k=n_loci)  # 提取下三角中偏移 n_loci 的对角线
        for i in range(len(INS_af_ls)) :
            for j in range(1,len(af_range_ls)) :
                if str(direct_corr[i]) == "nan" or direct_corr[i] < 0 :
                    continue
                if INS_af_ls[i] == 0 :
                    INS_HGSVC_af_res[0].append(direct_corr[i])
                    HGSVC_r2.append(direct_corr[i])
                    break
                elif INS_af_ls[i] > af_range_ls[j-1] and INS_af_ls[i] <= af_range_ls[j] :
                    INS_HGSVC_af_res[j-1].append(direct_corr[i])
                    HGSVC_r2.append(direct_corr[i])
                    break
                
        DEL_base_gt_ls = np.array(DEL_base_gt_ls)
        DEL_comp_gt_ls = np.array(DEL_comp_gt_ls)
        corr_matrix = np.corrcoef(DEL_comp_gt_ls, DEL_base_gt_ls, rowvar=False)
        n_loci = DEL_base_gt_ls.shape[1]
        direct_corr = np.diag(corr_matrix, k=n_loci)  # 提取下三角中偏移 n_loci 的对角线
        for i in range(len(DEL_af_ls)) :
            for j in range(1,len(af_range_ls)) :
                if str(direct_corr[i]) == "nan" or direct_corr[i] < 0 :
                    continue
                if DEL_af_ls[i] == 0 :
                    DEL_HGSVC_af_res[0].append(direct_corr[i])
                    HGSVC_r2.append(direct_corr[i])
                    break
                elif DEL_af_ls[i] > af_range_ls[j-1] and DEL_af_ls[i] <= af_range_ls[j] :
                    DEL_HGSVC_af_res[j-1].append(direct_corr[i])
                    HGSVC_r2.append(direct_corr[i])
                    break    
    
        HGSVC_file_fn = HGSVC_base.replace("chr1",chr) + "/fn.vcf.gz"
        fh = gzip.open(HGSVC_file_fn, 'rt') if HGSVC_file_fn.endswith('.gz') else open(HGSVC_file_fn, 'rt')
        for line in fh:
            if line == "" :
                continue
            line = line.strip()
            if line.startswith('#'):
                continue
            else :
                line_split_ls = line.split("\t")
                if "INS" in line_split_ls[7] :
                    af = float(line_split_ls[7].split("AF=")[1].split(";")[0])
                    for j in range(1,len(af_range_ls)) :
                        if af == 0 :
                            INS_HGSVC_af_res[0].append(0)
                            HGSVC_r2.append(0)
                            break
                        elif af > af_range_ls[j-1] and af <= af_range_ls[j] :
                            INS_HGSVC_af_res[j-1].append(0)
                            HGSVC_r2.append(0)
                            break
                if "DEL" in line_split_ls[7] :
                    af = float(line_split_ls[7].split("AF=")[1].split(";")[0])
                    for j in range(1,len(af_range_ls)) :
                        if af == 0 :
                            DEL_HGSVC_af_res[0].append(0)
                            HGSVC_r2.append(0)
                            break
                        elif af > af_range_ls[j-1] and af <= af_range_ls[j] :
                            DEL_HGSVC_af_res[j-1].append(0)
                            HGSVC_r2.append(0)
                            break
        fh.close()
    
    INS_cuteSVTrio_mean = [0 if len(x) == 0 else sum(x)/len(x) for x in INS_cuteSVTrio_af_res]
    DEL_cuteSVTrio_mean = [0 if len(x) == 0 else sum(x)/len(x) for x in DEL_cuteSVTrio_af_res]
    INS_HGSVC_mean = [0 if len(x) == 0 else sum(x)/len(x) for x in INS_HGSVC_af_res]
    DEL_HGSVC_mean = [0 if len(x) == 0 else sum(x)/len(x) for x in DEL_HGSVC_af_res]
    cuteSVTrio_r2_mean = sum(cuteSVTrio_r2)/len(cuteSVTrio_r2)
    HGSVC_r2_mean = sum(HGSVC_r2)/len(HGSVC_r2)
    print("cuteSVTrio_r2:",cuteSVTrio_r2_mean)
    print("HGSVC_r2:",HGSVC_r2_mean)

    data = []
    data.append(["sv_type","af","cuteSVTrio_r2","HGSVC_r2"])
    for i in range(len(INS_cuteSVTrio_af_res)) :
        data.append(["INS",i+1,INS_cuteSVTrio_mean[i],INS_HGSVC_mean[i]])
    for i in range(len(DEL_cuteSVTrio_af_res)) :
        data.append(["DEL",i+1,DEL_cuteSVTrio_mean[i],DEL_HGSVC_mean[i]])
    with open(write_file, mode='w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerows(data)
        f.close()

    #print([0 if len(x) == 0 else sum(x)/len(x) for x in total_cuteSVTrio_af_res])
    #print([0 if len(x) == 0 else sum(x)/len(x) for x in total_HGSVC_af_res])
