import pandas as pd
from pathlib import Path
import math
from tqdm import tqdm
import sys

spdi_dic = {'chr1':'NC_000001.11', 'chr2':'NC_000002.12', 'chr3':'NC_000003.12', 'chr4':'NC_000004.12', 'chr5':'NC_000005.10', 'chr6':'NC_000006.12',
        'chr7':'NC_000007.14', 'chr8':'NC_000008.11', 'chr9':'NC_000009.12', 'chr10':'NC_000010.11', 'chr11':'NC_000011.10', 'chr12':'NC_000012.12',
        'chr13':'NC_000013.11', 'chr14':'NC_000014.9', 'chr15':'NC_000015.10', 'chr16':'NC_000016.10', 'chr17':'NC_000017.11', 'chr18':'NC_000018.10',
        'chr19':'NC_000019.10', 'chr20':'NC_000020.11', 'chr21':'NC_000021.9', 'chr22':'NC_000022.11', 'chrX':'NC_000023.11', 'chrY':'NC_000024.10'}
        
ROOT_DIR = Path('/hpc/group/igvf/A549/extra_GCs/IGVF_var_preds/preds/')

def combine(label):
    headers = ['chrom','pos','ref','alt','log2_pred','spdi']

    with open(ROOT_DIR/label/'all.txt', 'w') as outfile:
        outfile.write("\t".join(headers) + "\n")
        for i in tqdm(range(1,268)):
            file = ROOT_DIR/label/('pred'+str(i)+'.txt')
            with open(file,'r') as infile:
                for line in infile:
                    line = line.strip().split('\t')
                    chrom = line[0].split(':')[0]
                    pos = line[1].split('=')[-1]
                    ref = line[2].split('=')[-1]
                    alt = line[3]
                    log2_pred = str(float(line[4])*math.log2(math.e))
                    spdi = spdi_dic[chrom]+':'+ pos+':'+ref+':'+alt # NC_000001.11:826009:T:A
                    out_line = [chrom, pos, ref, alt, log2_pred, spdi]
                    
                    outfile.write("\t".join(out_line) + "\n")

def process(label):
    all_df = pd.read_table(ROOT_DIR/label/'all.txt', sep = '\t')
    print(len(all_df))
    all_df = all_df.drop_duplicates(subset=['chrom','pos','ref','alt'])
    print(len(all_df))
    ref = all_df.loc[all_df['ref'] == all_df['alt'], ]
    ref = ref.drop(columns=['alt'])
    ref = ref.rename(columns={'log2_pred':'log2_ref_pred','spdi':'spdi_ref'})
    all_merged = pd.merge(all_df, ref, on=['chrom','pos','ref'], how='left')
    all_merged['effect'] = all_merged['log2_pred'] - all_merged['log2_ref_pred'] # log2FC
    all_merged = all_merged.loc[all_merged['ref'] != all_merged['alt'], ]
    print(len(all_merged))
    all_merged = all_merged[['chrom','pos','effect','spdi']]
    all_merged.to_csv(ROOT_DIR/label/'all.txt', sep = '\t', index=False)

def main(label):
    all_df = combine(label)
    process(label)

if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <label:A,B,C>\n")
(labels)=sys.argv[1].split(',')

for label in labels:
    print(label)
    main(label)
                