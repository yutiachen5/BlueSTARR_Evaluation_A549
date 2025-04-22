# split
# biased downsampling for training + unbiased downsampling for testing and validation

import pandas as pd
import numpy as np
from FastaWriter import FastaWriter
from tqdm import tqdm
import sys
import gzip
import os
from FastaWriter import FastaWriter

fastaWriter=FastaWriter()

def loadFasta(fasta_path, as_dict=False,uppercase=False, stop_at=None,
              revcomp=False):
    fastas = []
    seq = None
    header = None
    for r in (gzip.open(fasta_path) if fasta_path.endswith(".gz") else open(fasta_path)):
        if type(r) is bytes: r = r.decode("utf-8")
        r = r.strip()
        if r.startswith(">"):
            if seq != None and header != None:
                fastas.append([header, seq])
                if stop_at != None and len(fastas) >= stop_at:
                    break
            seq = ""
            header = r[1:]
        else:
            if seq != None:
                seq += r.upper() if uppercase else r
            else:
                seq = r.upper() if uppercase else r
    if stop_at != None and len(fastas) < stop_at:
        fastas.append([header, seq])
    elif stop_at == None:
        fastas.append([header, seq])
    if as_dict:
        return {h: s for h, s in fastas}
    if(revcomp):
        for rec in fastas:
            rc=generate_complementary_sequence(rec[1])
            rec[1]=rec[1]+"NNNNNNNNNNNNNNNNNNNN"+rc
    return pd.DataFrame({'location': [e[0] for e in fastas],
                         'sequence': [e[1] for e in fastas]})

def write(fasta, count, label, name, header):
    fasta.reset_index(drop=True, inplace=True)
    count.reset_index(drop=True, inplace=True)
    assert len(fasta) == len(count)

    os.makedirs(output_dir, exist_ok=True)
    out_fasta = open(output_dir+name+'.fasta',"wt")
    out_count = gzip.open(output_dir+name+'-counts.txt.gz',"wt")

    ls_loc = fasta['location']
    ls_seq = fasta['sequence']

    # header = "DNA=5\tRNA=4\n"
    print(header,end="",file=out_count)
    count.to_csv(out_count, sep='\t', index=False, header=None)

    for i in range(len(fasta)):
        fastaWriter.addToFasta(ls_loc[i],ls_seq[i],out_fasta)

    out_fasta.close()
    out_count.close()

def main(label, n_train, n_val, n_test):
    fasta = loadFasta(input_dir+label+'/all.fasta.gz')
    count = pd.read_csv(input_dir+label+'/all-counts.txt.gz', sep = '\t', compression="gzip", skiprows=1, header=None)

    header = gzip.open(input_dir+label+'/all-counts.txt.gz',"rt").readline()
    
    all_idx = np.arange(len(fasta))
    train_idx = np.random.choice(all_idx, size=int(n_train), replace=False)
    remaining_idx = np.setdiff1d(all_idx, train_idx)
    val_idx = np.random.choice(remaining_idx, size=int(n_val), replace=False)
    remaining_idx = np.setdiff1d(remaining_idx, val_idx)
    test_idx = np.random.choice(remaining_idx, size=int(n_test), replace=False)


    fasta_train = fasta.iloc[train_idx]
    count_train = count.iloc[train_idx] 
    fasta_val = fasta.iloc[val_idx]
    count_val = count.iloc[val_idx]
    fasta_test = fasta.iloc[test_idx]
    count_test = count.iloc[test_idx]

    write(fasta_train, count_train, label, 'train', header)
    write(fasta_val, count_val, label, 'validation', header)
    write(fasta_test, count_test, label, 'test', header)
    

if(len(sys.argv)!=7):
    exit(ProgramName.get()+" <input-dir:/.../> <output-dir:/.../> <label:A,B,C> <n-train> <n-val> <n-test>\n")
(input_dir,output_dir,labels,n_train,n_val,n_test)=sys.argv[1:]

labels = labels.split(',') 
for label in labels:
    main(label,int(n_train),int(n_val),int(n_test))

# command line:
# python /work/igvf-pm/A549/extra_GCs/downsample_unbiased.py /work/igvf-pm/A549/extra_GCs/ /work/igvf-pm/A549/extra_GCs/ AZD2906,AZD9567,CORT108297,CpdA,GW870086,Hydrocortisone,Mapracorat,RU486,ZK216348 1600000 500000 500000
# python /work/igvf-pm/A549/extra_GCs/downsample_unbiased.py /work/igvf-pm/A549/extra_GCs/data-normalized/ /work/igvf-pm/A549/extra_GCs/data-normalized/ AZD2906,AZD9567,CORT108297,CpdA,GW870086,Hydrocortisone,Mapracorat,RU486,ZK216348 1600000 500000 500000
# python /work/igvf-pm/A549/extra_GCs/prepare_input_unbiased.py /work/igvf-pm/K562/full-set/600bp/data-normalized/ /work/igvf-pm/K562/full-set/600bp/data-normalized/split-unbiased/ '' 700000 200000 200000