import pandas as pd
import gzip

## check if 2 set overlap with each other

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
                         'idx': [e[0].split(' /coord')[0].split('>')[-1] for e in fastas],  # >7485398 /coord=chr10:47800-48100
                         'sequence': [e[1] for e in fastas]})

df1 = loadFasta('/datacommons/igvf-pm/A549/full-set/DMSO-200/600-bases/data-normalized/train.fasta.gz')
df2 = loadFasta('/datacommons/igvf-pm/A549/full-set/DMSO-200/600-bases/data-normalized/test.fasta.gz')
df3 = loadFasta('/datacommons/igvf-pm/A549/full-set/DMSO-200/600-bases/data-normalized/validation.fasta.gz')
same_idx_train_test = list(set(list(df1['idx'])).intersection(set(list(df2['idx']))))
same_idx_train_validation = list(set(list(df1['idx'])).intersection(set(list(df3['idx']))))
same_idx_test_validation = list(set(list(df2['idx'])).intersection(set(list(df3['idx']))))
print(same_idx_train_test)
print(same_idx_train_validation)
print(same_idx_test_validation)

df1 = loadFasta('/datacommons/igvf-pm/A549/full-set/Dex-200/600-bases/data-normalized/train.fasta.gz')
df2 = loadFasta('/datacommons/igvf-pm/A549/full-set/Dex-200/600-bases/data-normalized/test.fasta.gz')
df3 = loadFasta('/datacommons/igvf-pm/A549/full-set/Dex-200/600-bases/data-normalized/validation.fasta.gz')

same_idx_train_test = list(set(list(df1['idx'])).intersection(set(list(df2['idx']))))
same_idx_train_validation = list(set(list(df1['idx'])).intersection(set(list(df3['idx']))))
same_idx_test_validation = list(set(list(df2['idx'])).intersection(set(list(df3['idx']))))
print(same_idx_train_test)
print(same_idx_train_validation)
print(same_idx_test_validation)