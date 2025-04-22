# process filtered variance from /work/igvf-pm/IGVF-DACC/IGVFFI2231RETG.vcf.gz and make chunks for prediction

import pandas as pd
import vcfpy
import sys
import gzip

n = 0
chunk = 1
with gzip.open("/work/igvf-pm/IGVF-DACC/IGVFFI2231RETG.vcf.gz", "r") as infile:#, open('/work/igvf-pm/A549/extra_GCs/IGVF_var_preds/chunks/VARs'+str(chunk)+'.txt',"wt") as outfile:
    for line in infile:
        line = line.decode("utf-8")
        if line.startswith("#") == False:
            record = line.strip().split("\t")
            chrom = record[0]
            pos = str(int(record[1])-1) # convert 1-based index to 0-based index to fit 2bittofa
            ref = record[3].upper()
            alt = record[4].upper()
            out = chrom+':'+pos+':'+'ref='+ref+':'+ref+','+alt
            if ref == 'N':
                print(out)
            # print(out, file=outfile)
            n+=1
            
            # if n%100000==0:
            #     outfile.close()
            #     chunk+=1
            #     outfile=open('/work/igvf-pm/A549/extra_GCs/IGVF_var_preds/chunks/VARs'+str(chunk)+'.txt',"wt")
print(n)



