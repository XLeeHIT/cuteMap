# -*- coding: utf-8 -*-
# 该脚本按照染色体筛选sv记录

import sys
import os
import csv
import gzip
import copy
import math
import numpy as np
from scipy.interpolate import CubicSpline
import random

child_ls = ["HG002","HG005","HG00514","HG00733","NA12878","NA19240","NA19836","HG03522","NA12877","HG01891","HG02668","HG03371","HG03453","NA19705"]
father_ls = ["HG003","HG006","HG00512","HG00731","NA12891","NA19238","NA19834","HG03521","NA12889","HG01890","HG02666","","","NA19703"]
mother_ls = ["HG004","HG007","HG00513","HG00732","NA12892","NA19239","NA19835","HG03520","NA12890","","","HG03369","HG03452",""]

child_trio_ls = ["HG002","HG005","HG00514","HG00733","NA12878","NA19240","NA19836","HG03522","NA12877"]
father_trio_ls = ["HG003","HG006","HG00512","HG00731","NA12891","NA19238","NA19834","HG03521","NA12889"]
mother_trio_ls = ["HG004","HG007","HG00513","HG00732","NA12892","NA19239","NA19835","HG03520","NA12890"]

child_duo_ls = ["HG01891","HG02668","HG03371","HG03453","NA19705"]
father_duo_ls = ["HG01890","HG02666","","","NA19703"]
mother_duo_ls = ["","","HG03369","HG03452",""]

# 只有这7个完整家庭同时出现在1kGP的样本集中
child_1kgp_trio_ls = ["HG00405","HG00408","HG00514","HG00733","NA12878","NA19240","NA19836","HG03522","NA12877"]
father_1kgp_trio_ls = ["HG00403","HG00406","HG00512","HG00731","NA12891","NA19238","NA19834","HG03521","NA12889"]
mother_1kgp_trio_ls = ["HG00404","HG00407","HG00513","HG00732","NA12892","NA19239","NA19835","HG03520","NA12890"]

cuteSVTrio_1kgp_member_ls = ["HG00512","HG00513","HG00514","HG00731","HG00732","HG00733",
                             "HG01890","HG01891","HG02666","HG02668","HG03369","HG03371",
                             "HG03452","HG03453","HG03520","HG03521","HG03522","NA12877",
                             "NA12878","NA12889","NA12890","NA12891","NA12892","NA19238",
                             "NA19239","NA19240","NA19703","NA19705","NA19834","NA19835","NA19836"]

HGSVC_1kgp_member_ls = ["HG00171","HG00268","HG00358","HG00512","HG00513","HG00514",
                        "HG00731","HG00732","HG00733","HG01457","HG01890","HG02282","HG02587",
                        "HG02666","HG02769","HG03009","HG03371","HG03452","HG03520",
                        "HG03732","NA18534","NA19036","NA19129","NA19238","NA19239",
                        "NA19240","NA19331","NA19384","NA19705","NA19836","NA19983"]

cuteSVTrio_member_ls = ["HG002","HG003","HG004","HG005","HG006","HG007",
                        "HG00514","HG00733","NA12878","NA19240","NA19836","HG03522","NA12877",
                        "HG00512","HG00731","NA12891","NA19238","NA19834","HG03521","NA12889",
                        "HG00513","HG00732","NA12892","NA19239","NA19835","HG03520","NA12890",
                        "HG01891","HG02668","HG03371","HG03453","NA19705",
                        "HG01890","HG02666","NA19703","HG03369","HG03452",
                        "HG01889","HG02667","HG03370","HG03451","NA19704"]

# 这个样本列表是从1kGP/1kGP_chr1_SV.15.target.vcf.gz中提取出来的
target_member_ls = ["HG00242","HG01148","HG01674","HG02017","HG02628","HG02682","HG02696","NA06994","NA10843","NA10859","NA18862","NA18908","NA19027","NA19360","NA19438","NA20804"]

cuteSVTrio_relative_member_ls = ["HG002","HG003","HG004","HG005","HG006","HG007",
                        "HG00514","HG00733","NA12878","NA19240","NA19836","HG03522","NA12877",
                        "HG00512","HG00731","NA12891","NA19238","NA19834","HG03521","NA12889",
                        "HG00513","HG00732","NA12892","NA19239","NA19835","HG03520","NA12890",
                        "HG01891","HG02668","HG03371","HG03453","NA19705",
                        "HG01890","HG02666","NA19703","HG03369","HG03452",
                        "HG01889","HG02667","HG03370","HG03451","NA19704"]

HGSVC_relative_member_ls = ["HG00171","HG00268","HG00358","HG00512","HG00513","HG00514",
                   "HG00731","HG00732","HG00733","HG01890","HG02282","HG02587",
                   "HG02666","HG02769","HG03009","HG03371","HG03452","HG03520",
                   "HG03732","NA18534","NA19036","NA19129","NA19238","NA19239",
                   "NA19240","NA19331","NA19384","NA19705","NA19836","NA19983",
                   "HG01455","HG01456","HG02281","HG02585","HG02586","HG02770","HG02768","HG02667",
                   "HG02668","HG03727","HG03721","NA19128","NA19127","NA19982","NA19713"]


area_ls = ["AFR","EAS","EAS","AMR","EUR","AFR","AFR","AFR","EUR","AFR","AFR","AFR","AFR","AFR"]

chr_ls = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]

# 将vcf文件中某个染色体的变异记录提取出来
def pick_chr_from_vcf(read_vcf,write_vcf,chr) :
    line_ls = []
    fh = gzip.open(read_vcf, 'rt') if read_vcf.endswith('.gz') else open(read_vcf, 'rt')
    for line in fh:
        if line == "" :
            continue
        line = line.strip()
        if line.startswith('#'):
            line_ls.append(line)
        else :
            line_split_ls = line.split("\t")
            if line_split_ls[0].replace("chr","") == chr.replace("chr","") :
                line_ls.append(line)
    fh.close()
    with open(write_vcf,'w') as f:
        for line in line_ls :
            f.write(line+"\n")
    f.close()

def main(argv):
    pick_chr_from_vcf("cutesv_trio/merge.truvari.rectified.vcf.gz".replace("chr1",argv[0]),"cutesv_trio/merge.truvari.rectified.chr1.vcf".replace("chr1",argv[0]),argv[0].replace("chr",""))
    pick_chr_from_vcf("freeze4/merge.truvari.rectified.vcf.gz".replace("chr1",argv[0]),"freeze4/merge.truvari.rectified.chr1.vcf".replace("chr1",argv[0]),argv[0].replace("chr",""))


if __name__ == "__main__":
    main(sys.argv[1:])