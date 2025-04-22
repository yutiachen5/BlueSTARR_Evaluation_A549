import numpy as np
import pandas as pd
import gzip
from FastaWriter import FastaWriter
from scipy.stats import lognorm

AP1 = 'TGAGTCAT'
GR = 'GAACATTATGTTC'

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
    return pd.DataFrame({'location': [e[0].split('=')[-1] for e in fastas],
                         'sequence': [e[1] for e in fastas]})

def get_rdm_seq(rdm_seq_path, p): # sample random sequences and generate random scores
    df = pd.read_csv(rdm_seq_path, sep='\t', header=None)
    df_sampled = df.sample(frac=1-p, random_state=42)
    score = lognorm(s=1, scale=np.exp(0)).rvs(size=len(df_sampled))
    df_sampled.columns = ['location','sequence']
    df_sampled['score'] = score
    return df_sampled
    
def calculate_activation(dist):
    if dist <= 200:
        score = 200 - dist
    else:
        score = dist - 200
    return min(score, 200)
    
def insert_motifs(low_act_seq_path, p, seq_len=300): # insert GR and AP1 to low activation sequences at varying distances
    df = pd.read_csv(low_act_seq_path, sep='\t', header=None)
    df.columns = ['location','sequence']
    seq = list(df['sequence'])
    loc = list(df['location'])
    ls_seq = []
    ls_loc = []
    ls_score = []

    n = int(170_000*p)
    for i in range(n):
        sampled = np.random.choice(df.index, size=1, replace=True)
        sequence = seq[sampled[0]]
        max_dist = seq_len - len(GR)
        dist = np.random.randint(len(AP1), max_dist + 1)  # dist: from begining of a motif to the begining of another motif
        
        if dist == max_dist: pos_ap1 = 0
        else: pos_ap1 = np.random.randint(0, seq_len - len(GR) - dist)
        pos_gr = pos_ap1 + dist

        sequence = sequence[:pos_ap1] + AP1 + sequence[pos_ap1 + len(AP1):]
        sequence = sequence[:pos_gr] + GR + sequence[pos_gr + len(GR):]
        
        ls_seq.append(sequence)
        ls_loc.append(loc[sampled[0]])
        ls_score.append(calculate_activation(dist))
        
    installed = pd.DataFrame({'location':ls_loc,'sequence':ls_seq, 'score':ls_score})

    return installed

def write(df, outdir):
    fastaWriter=FastaWriter()
    output_fasta = open(outdir+".fasta","wt") 
    for i in range(len(df)):
        fastaWriter.addToFasta(df.loc[i,'location'], df.loc[i,'sequence'], output_fasta)
    
    output_counts = outdir+"-counts.txt"
    with open(output_counts, 'w') as file:
        file.write("DNA=1\tRNA=1\n")  
        score = pd.DataFrame({'DNA':[1]*len(df), 'RNA':list(df['score'])})
        score.to_csv(file, sep='\t', index=False, header=None)   

p = 0.02
rdm_seq_score = get_rdm_seq('/datacommons/igvf-pm/A549/GR-AP1/simulated-seq/data/Dex-200/diluting/rdm_seq.txt', p) 
low_act_installed_score = insert_motifs('/datacommons/igvf-pm/A549/GR-AP1/simulated-seq/data/Dex-200/diluting/low_act_5000_train.txt', p)

# concat random seq and low-activation seq with scores 
train = pd.concat([rdm_seq_score, low_act_installed_score], axis=0, ignore_index=True)
train.columns = ['location','sequence','score']

# generate random score for testing and validation set
# test = loadFasta('/datacommons/igvf-pm/A549/GR-AP1/simulated-seq/data/Dex-200/test.fasta.gz')
# test['score'] = lognorm(s=1, scale=np.exp(0)).rvs(size=len(test))
# validation = loadFasta('/datacommons/igvf-pm/A549/GR-AP1/simulated-seq/data/Dex-200/validation.fasta.gz')
# validation['score'] = lognorm(s=1, scale=np.exp(0)).rvs(size=len(validation))

# write data to .fasta and .counts for training
write(train, '/datacommons/igvf-pm/A549/GR-AP1/simulated-seq/data/Dex-200/diluting/train-'+str(1-p))
# write(test, '/datacommons/igvf-pm/A549/GR-AP1/simulated-seq/data/Dex-200/diluting/test')
# write(validation, '/datacommons/igvf-pm/A549/GR-AP1/simulated-seq/data/Dex-200/diluting/validation')
