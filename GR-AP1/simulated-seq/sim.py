import numpy as np
import math
import pandas as pd
import sys
import gzip
from tqdm import tqdm
from FastaWriter import FastaWriter

forward_ap1 = 'TGAGTCAT'
forward_gr  = 'GAACATTATGTTC'
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

def calculate_activation(dist):
    if dist <= 200:
        score = 200 - dist
    else:
        score = dist - 200
    return min(score, 200)

def insert_motifs(fasta, label, seq_len=300):
    seq = fasta['sequence']
    loc = fasta['location']
    ls_score = []
    ls_distance = []

    output_fasta = open(outdir + label + ".fasta","wt") 
    for i in tqdm(range(len(fasta))):
        sequence = seq[i]
        max_dist = seq_len - len(forward_gr)
        dist = np.random.randint(len(forward_ap1), max_dist + 1)  # dist: from begining of a motif to the begining of another motif
        ls_distance.append(dist)
        
        if dist == max_dist: pos_ap1 = 0
        else: pos_ap1 = np.random.randint(0, seq_len - len(forward_gr) - dist)
        pos_gr = pos_ap1 + dist

        sequence = sequence[:pos_ap1] + forward_ap1 + sequence[pos_ap1 + len(forward_ap1):]
        sequence = sequence[:pos_gr] + forward_gr + sequence[pos_gr + len(forward_gr):]
        ls_score.append(calculate_activation(dist))

        fastaWriter.addToFasta(loc[i],sequence,output_fasta)
        
    score = pd.DataFrame({'DNA':[1]*len(fasta), 'RNA':ls_score})
    score_dist = pd.DataFrame({'score':ls_score, 'distance':ls_distance})

    output_score = outdir + label + '-counts.txt'
    with open(output_score, 'w') as file:
        file.write("DNA=1\tRNA=1\n")  
        score.to_csv(file, sep='\t', index=False, header=None) 
    
    output_score_dist = outdir + label + '-scores-dist.txt'
    with open(output_score_dist, 'w') as file:
        score_dist.to_csv(file, sep='\t', index=False, header=None) 

    return fasta, score


def main(subdir, outdir):
    train = loadFasta(subdir+'train.fasta.gz')
    test = loadFasta(subdir+'test.fasta.gz')
    validation = loadFasta(subdir+'validation.fasta.gz')

    train_fasta, train_score = insert_motifs(train, "train")
    test_fasta, test_score = insert_motifs(test, "test")
    validation_fasta, validation_score = insert_motifs(validation, "validation")



if(len(sys.argv)!=3):
    exit(ProgramName.get()+"<data-subdir> <output-dir> \n")
(subdir,outdir)=sys.argv[1:]
main(subdir,outdir)