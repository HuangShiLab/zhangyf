import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

ref_id_list = []
GeneName = {}
FPKM = {}
GeneSize = {}

with open (input_file, "r") as fr:
    lines = fr.readlines()
    for line in lines:
        s = line.rstrip().split("\t")
        sp = int(s[3])
        ep = int(s[4])
        s = s[-1]
        s = s.replace('"', '')
        s = s.split(";")
        s = [d.lstrip() for d in s]
        s = [d.split(" ") for d in s]
        
        ref_id = fpkm = gene_name = ""
        
        for data in s:
            if data[0] == "reference_id":
                ref_id = data[1]
            elif data[0] == "FPKM":
                fpkm = data[1]
            elif data[0] == "ref_gene_name":
                gene_name = data[1]

            if ref_id[0:4] == "rna-":
                ref_id = ref_id.replace('rna-', 'gene-')
            
        if ref_id not in ref_id_list and fpkm != "":
            ref_id_list.append(ref_id)
            
            FPKM[ref_id] = fpkm
            GeneSize[ref_id] = ep - sp + 1
            
            if gene_name != "":
                GeneName[ref_id] = gene_name
            else:
                GeneName[ref_id] = ref_id.replace('gene-', '')

with open(output_file, "w") as fw:
    fw.write("GeneID\tGeneName\tgene_size\tFPKM\n")
    for ref_id in ref_id_list:
        fw.write(ref_id + "\t" + GeneName[ref_id] + "\t" + str(GeneSize[ref_id]) + "\t" + FPKM[ref_id] + "\n")

