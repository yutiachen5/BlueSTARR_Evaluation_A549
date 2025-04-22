# adapted from: https://gitlab.oit.duke.edu/igvf/bluestarr-viz/-/blob/main/sim/Mixture-biased-sampling.ipynb
# for simulated sequences

import numpy as np
import pandas as pd
import sys
import gzip
import random
from FastaWriter import FastaWriter
import ProgramName
from Rex import Rex
from scipy.stats import lognorm
from scipy.stats import gaussian_kde
from scipy.interpolate import PchipInterpolator
rex=Rex()

RANGE=0.99

class HistKDE:
    def __init__(self, samples, bins=1000):
        (densities, bin_edges) = np.histogram(samples, bins=bins, density=True)
        self.bins = bin_edges
        self.mids = bin_edges[1:] - np.diff(bin_edges)/2
        self.probs = densities * np.diff(bin_edges)
        self.cum_probs = np.cumsum(self.probs)
        self.interp_pdf = PchipInterpolator(self.mids, self.probs, extrapolate=True)
        self.interp_cdf = PchipInterpolator(self.mids, self.cum_probs, extrapolate=True)

    def __call__(self, x):
        return self.pdf(x)

    def pdf(self, x, interpolate=True):
        if interpolate:
            return self.interp_pdf(x)
        else:
            bin_indices = np.digitize(x, self.bins[1:], right=False).clip(0, len(self.probs)-1)
            return self.probs[bin_indices]

    def cdf(self, x, interpolate=True):
        if interpolate:
            return self.interp_cdf(x)
        else:
            bin_indices = np.digitize(x, self.bins[1:], right=False).clip(0, len(self.cum_probs)-1)
            return self.cum_probs[bin_indices]

def scores_and_unif(scores: np.ndarray):
    if scores.ndim > 1:
        scores = scores[:,0]
    (min, max) = (int(scores.min()), int(scores.max()))
    unif = (1 / (max - min + 1))
    return scores, unif

# def accept_hist_pdf(scores, kde=None, hist_bins=1000, interpolate=True):
#     # scores = np.array(scores)
#     scores, _ = scores_and_unif(scores)
#     if not kde:
#         kde = HistKDE(scores, bins=hist_bins)
#     return 1 - kde.pdf(scores, interpolate=interpolate)

def accept_hist_pdf(scores, kde=None, hist_bins=10000, interpolate=True):
    scores, _ = scores_and_unif(scores)
    if not kde:
        kde = HistKDE(scores, bins=hist_bins)
    acceptance_probs = kde.pdf(scores, interpolate=interpolate)
    biased_probs = acceptance_probs**2
    biased_probs /= biased_probs.max()
    return 1 - biased_probs

def accept_hist_cdf(scores, kde=None, hist_bins=1000, interpolate=True):
    scores, _ = scores_and_unif(scores)
    if not kde:
        kde = HistKDE(scores, bins=hist_bins)
    return kde.cdf(scores, interpolate=interpolate)

def accept_cdf(scores, distclass=lognorm, postfunc=lambda x: x, **cdf_params):
    scores, _ = scores_and_unif(scores)
    return postfunc(distclass.cdf(scores, **cdf_params))

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
    
def sample_scores_accept(scores, accept_probs, frac=0.1):
    rng = np.random.default_rng()
    accept_probs = accept_probs / accept_probs.sum()
    indices = np.arange(len(scores))
    chosen_indices = rng.choice(indices, size=int(len(scores) * frac), replace=False, p=accept_probs)
    chosen_scores = scores[chosen_indices]
    return chosen_indices, chosen_scores

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=8):
    exit(ProgramName.get()+" <in:fasta.gz> <in:counts-scores.txt.gz> <#bins> <desired-sample-size> <out-dir> <label> <acceptance func>")
(fastaFile,countsFile,numBins,N,outDir,fileLabel,accFunc)=sys.argv[1:]

numBins=int(numBins)
N=int(N)
outCountFile=outDir+"/"+fileLabel+"-counts.txt.gz"
outFastaFile=outDir+"/"+fileLabel+".fasta.gz"

# Create output files
OUT_COUNTS=gzip.open(outCountFile,"wt")
OUT_FASTA=gzip.open(outFastaFile,"wt")

# Downsample the data
IN_COUNTS=gzip.open(countsFile,"rt")
header=IN_COUNTS.readline()
if(not rex.find("DNA=(\d+)\s+RNA=(\d+)",header)):
    raise Exception("Can't parse header in counts file: "+header)
print(header,end="",file=OUT_COUNTS)
DNA_REPS=int(rex[1]); RNA_REPS=int(rex[2])
fastaWriter=FastaWriter()

fasta = loadFasta(fastaFile)
counts = pd.read_csv(countsFile, sep='\t', header=None, skiprows=1, compression="gzip")
counts.columns = ['DNA','RNA']
score = np.array(counts['RNA'])
histkde = HistKDE(score, bins=numBins)
fraction = N/len(counts)

if accFunc == 'BlueSTARR-like':
    chosen_indices, chosen_scores = sample_scores_accept(score, frac=fraction, accept_probs=accept_hist_pdf(score, histkde, interpolate=False))  # BlueSTARR like
elif accFunc == 'ECDF':
    chosen_indices, chosen_scores = sample_scores_accept(score, frac=fraction, accept_probs=accept_hist_cdf(score, histkde, interpolate=True))  # ECDF
elif accFunc == 'lognormal':
    chosen_indices, chosen_scores = sample_scores_accept(score, frac=fraction, accept_probs=accept_cdf(score, s=1, scale=1))  # lognormal CDF
elif accFunc == 'lognormal^2':
    chosen_indices, chosen_scores = sample_scores_accept(score, frac=fraction, accept_probs=accept_cdf(score, s=1, scale=1, postfunc=np.square)) # lognormal CDF^2
elif accFunc == 'lognormal^10':
    chosen_indices, chosen_scores = sample_scores_accept(score, frac=fraction, accept_probs=accept_cdf(score, s=1, scale=1, postfunc=lambda x: np.power(x, 10))) # lognormal CDF^10
else:
    raise Exception('choice of acceptance function: BlueSTARR-like, ECDF, lognormal, lognormal^2, lognormal^10')


chosen_fasta = fasta.iloc[chosen_indices]
chosen_fasta = chosen_fasta.reset_index(drop=True)
# print how many percentage of positive samples were selected (top p% were positive)
postfrac = np.sum(chosen_indices<(len(counts)*0.005))/N*100
print(f'{postfrac:.2f}% positive samples after biased down-sampling by acceptance function {accFunc}')

loc = list(chosen_fasta['location'])
seq = list(chosen_fasta['sequence'])
for i in range(len(chosen_fasta)):
    fastaWriter.addToFasta(loc[i],seq[i],OUT_FASTA)

chosen_scores = pd.DataFrame({'DNA':[1]*len(chosen_fasta), 'RNA':chosen_scores})
chosen_scores.to_csv(OUT_COUNTS, sep='\t', header=None, index=False)

IN_COUNTS.close(); OUT_COUNTS.close(); OUT_FASTA.close()

# python /hpc/group/igvf/A549/full-set/BlueSTARR/downsample-biased-sim.py /hpc/group/igvf/A549/GR-AP1/simulated-seq/data/mound/all-train-0.005.fasta.gz /hpc/group/igvf/A549/GR-AP1/simulated-seq/data/mound/all-train-0.005-counts.txt.gz 100 1700000 /hpc/group/igvf/A549/GR-AP1/simulated-seq/data/mound train-lognormal10 lognormal^10
