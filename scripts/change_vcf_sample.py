import sys
import gzip
# 修改9+5个家庭的trio vcf的样本名，并进行初步的筛选
def change_vcf_sample(read_vcf,write_vcf,sample_name) :
    line_ls = []
    chr_ls = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
    fh = gzip.open(read_vcf, 'rt') if read_vcf.endswith('.gz') else open(read_vcf, 'rt')
    for line in fh:
        if line == "" :
            continue
        line = line.strip()
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                line_split_ls = line.split("\t")[:9] + sample_name.split(",")
                line_ls.append("\t".join(line_split_ls))
            elif line.startswith('##contig'):
                if line.split("ID=")[1].split(",")[0] not in chr_ls :
                    continue
                else :
                    line_ls.append(line.replace("chr",""))
            else :
                line_ls.append(line)
        else :
            line_split_ls = line.split("\t")
            if "PASS" not in line_split_ls[7] :
                continue
            line_split_ls[6] = "PASS"
            line_split_ls[0] = line_split_ls[0].replace("chr","")
            if line_split_ls[0] not in chr_ls :
                continue
            INFO_ls = line_split_ls[7].split(";")
            new_INFO_ls = []
            for INFO in INFO_ls :
                if "SVTYPE=" in INFO or "SVLEN=" in INFO or "END=" in INFO:
                    new_INFO_ls.append(INFO)
            line_split_ls[7] = ";".join(new_INFO_ls)
            hg_gt_index = line_split_ls[8].split(":").index("HP_GT")
            is_variant = False
            for i in range(9,len(line_split_ls)) :
                line_split_ls[i] = line_split_ls[i].split(":")[hg_gt_index]
                if "1" in line_split_ls[i][0:3] :
                    #print(line_split_ls[i][0:3])
                    is_variant = True
            line_split_ls[8] = "GT"
            if is_variant :
                #print(line_split_ls)
                if abs(int(line_split_ls[7].split("SVLEN=")[1].split(";")[0])) >= 50 :
                    line_ls.append("\t".join(line_split_ls))
    fh.close()
    with open(write_vcf,'w') as f:
        for line in line_ls :
            f.write(line+"\n")
    f.close()


def main(argv):
    change_vcf_sample(argv[0],argv[1],argv[2])

if __name__ == "__main__":
    main(sys.argv[1:])