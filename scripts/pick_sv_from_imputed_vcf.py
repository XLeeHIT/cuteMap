# -*- coding: utf-8 -*-
# 从imputed变异中筛选sv

import sys
import os
import csv
import gzip
import copy
import math
import numpy as np
from scipy.interpolate import CubicSpline
import random

def pick_sv_from_imputed_vcf(read_vcf,write_vcf) :
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
            if "IMPUTED" in line_split_ls[7].split(";") and (len(line_split_ls[3]) >= 50 or len(line_split_ls[4]) >= 50) :
                line_ls.append(line)
    fh.close()
    with open(write_vcf,'w') as f:
        for line in line_ls :
            f.write(line+"\n")
    f.close()


def main(argv):
    pick_sv_from_imputed_vcf(argv[0],argv[1])
    

if __name__ == "__main__":
    main(sys.argv[1:])