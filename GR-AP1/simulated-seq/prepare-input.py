# prepare the input file for model training in simulation analysis

import pandas as pd
import numpy as np
from scipy.stats import lognorm
from FastaWriter import FastaWriter
import gzip
import sys

MOTIF1 = 'TGAGTCAT' # AP1
MOTIF2 = 'GAACATTATGTTC' # GR
START = 10 # minimum position of MOTIF1

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

def v_shape(dist):
    score = np.where(dist <= 200, 200 - dist, dist - 200)
    return np.minimum(score, 200)

def exp_fit_1(dist):
    a, b, c = 3.83942702e-03, -9.39553441e-01, 4.61986242e+01
    return a * dist**2 + b * dist + c

def exp_fit_2(dist):
    a1, b1, c1 = 3.39285714e-03, -4.32500000e-01, 2.74669643e+01
    a2, b2, c2 = 1.44791666e-02, -5.25383925e+00, 4.75250250e+02

    part1 = a1 * dist**2 + b1 * dist + c1
    part2 = a2 * dist**2 + b2 * dist + c2
    return np.where(dist < 123, part1, part2)

def hill_activation(dist: np.ndarray,
                    y_max: float = 10,
                    y_min: float = 1, 
                    x_min: float = 10,
                    x_max: float = 300):

    if not (y_min < y_max):
        raise ValueError(f"y_min must be less than {y_max}. Got {y_min}")
        
    k = (y_max - y_min)/(x_min - x_max)
    b = y_max - x_min * k

    return dist * k + b

def cliff_activation(dist: np.ndarray,
                     y_max: float = 10,
                     y_min: float = 1,
                     x_intersection: float = 150,
                     y_intersection: float = 1.0,
                     x_min: float = 10,
                     x_max: float = 300):

    if not (x_min < x_intersection < x_max):
        raise ValueError(f"x_intersection must be in ({x_offset}, 300). Got {x_intersection}")
    if not (y_min < y_intersection < y_max):
        raise ValueError(f"y_intersection must be in ({y_min}, {y_max}). Got {y_intersection}")
        
    x1, y1 = x_min, y_max
    x2, y2 = x_max, y_min

    k1 = (y_intersection - y1) / (x_intersection - x1)  
    k2 = (y2 - y_intersection) / (x2 - x_intersection)

    y = np.zeros_like(dist)

    left_mask = dist <= x_intersection
    right_mask = dist > x_intersection

    y[left_mask] = y1 + k1 * (dist[left_mask] - x1)
    y[right_mask] = y_intersection + k2 * (dist[right_mask] - x_intersection)

    k1 = (y_max - y_intersection)/(x_min - x_intersection)
    b1 = y_max - x_min * k1
    k2 = (y_intersection - y_min)/(x_intersection - x_max)
    b2 = y_intersection - x_intersection * k2

    part1 = k1 * dist + b1
    part2 = k2 * dist + b2

    return np.where(dist < x_intersection, part1, part2)

def mound_activation(dists: np.ndarray,
                     d0: float,
                     k: float,
                     decay_rate: float, 
                     height: float,
                     x_offset: float=10) -> np.ndarray:
    """
    Generate a mound-shaped activation signal with a rapid rise to a peak 
    followed by a gentle decay.

    Parameters
    ----------
    dists : np.ndarray
        Array of distances (in basepairs) over which the activation is computed.
    d0 : float
        Distance at which the rapid rise reaches its midpoint.
    k : float
        Steepness of the sigmoid function controlling the rapid rise.
    decay_rate : float
        Exponential decay rate controlling the gentle decline after the peak.
    height : float
        Maximum height of the activation signal.
    x_offset : float, optional
        Offset applied to the distances, by default 10. Normally this is the minimum distance.

    Returns
    -------
    np.ndarray
        Array of activation values corresponding to the input distances.
    """
    sigmoid = 1 / (1 + np.exp(-k * (dists - d0)))
    decay = np.exp(-decay_rate * (dists-x_offset))
    return height * sigmoid * decay

def bowl_activation(dists: np.ndarray,
                    y_max0: float,
                    y_max1: float,
                    bowl_depth: float,
                    bowl_width: int=None,
                    x_offset: float=10) -> np.ndarray:
    """
    Generate a bowl-shaped activation signal with high values at short and long 
    distances and a parabolic decay and rise in between.

    Parameters
    ----------
    dists : np.ndarray
        Array of distances (in basepairs) over which the activation is computed.
    y_max0 : float
        Maximum activation value at the start of the distance range.
    y_max1 : float
        Maximum activation value at the end of the distance range.
    bowl_depth : float
        Depth of the bowl, representing the minimum activation value at the center.
    bowl_width : int, optional
        Width of the bowl, representing the distance range over which the parabolic 
        decay occurs. If None, it defaults to the range of `dists` minus `x_offset`.
    x_offset : float, optional
        Offset applied to the distances, by default 10. Normally, this is the
        minimum distance in the input array.

    Returns
    -------
    np.ndarray
        Array of activation values corresponding to the input distances.
    """
    if bowl_width is None:
        bowl_width = dists[-1] - x_offset
    
    xm = (2 * x_offset + bowl_width) / 2
    ym = min(y_max0, y_max1) - abs(bowl_depth)  # The minimum y at the center (bowl bottom)

    # Fit a parabola: y = a(x - xm)^2 + ym
    # Make sure it passes through (x0, y0) and (x1, y1)
    a = (y_max0 - ym) / ((x_offset - xm) ** 2)

    return a * (dists - xm) ** 2 + ym

def periodicity_activation(dists: np.ndarray,
                           period: float,
                           amplitude: float,
                           ampl_decay_rate,
                           y_offset: float,
                           offset_decay_rate: float,
                           x_offset: int=10) -> np.ndarray:
    """
    Generate a periodic activation signal with an amplitude and offset that decay 
    exponentially with distance.

    Parameters
    ----------
    dists : np.ndarray
        Array of distances (in basepairs) over which the activation is computed.
    period : float
        Period (in basepairs) of the sine wave representing the periodicity.
    amplitude : float
        Initial amplitude of the sine wave.
    ampl_decay_rate : float
        Exponential decay rate of the amplitude with distance.
    y_offset : float
        Initial vertical offset of the sine wave.
    offset_decay_rate : float
        Exponential decay rate of the vertical offset with distance.
    x_offset : int, optional
        Offset applied to the distances, by default 10. Normally, this is the
        distance at which the activation starts to be non-zero.

    Returns
    -------
    np.ndarray
        Array of activation values corresponding to the input distances.
    """
    offs_dists = dists - x_offset
    envelope = amplitude * np.exp(-ampl_decay_rate * offs_dists)
    offset_envelope = y_offset * np.exp(-offset_decay_rate * offs_dists)
    freq = 1 / period
    # set distances below x_offset to zero
    envelope[offs_dists < 0] = 0
    offset_envelope[offs_dists < 0] = 0
    return envelope * np.sin(2 * np.pi * freq * offs_dists) + offset_envelope
    
def insert_motifs(df, act_func, seq_len=300): 
    seq = list(df['sequence'])
    loc = list(df['location'])
    ls_seq = []
    ls_loc = []
    ls_dist = []

    for i in range(len(df)):
        sequence = seq[i]
        max_dist = seq_len - len(MOTIF2) - START
        dist = np.random.randint(len(MOTIF1), max_dist + 1)  # dist: from begining of a motif to the begining of another motif
        
        if dist == max_dist: 
            pos_motif1 = START
        else: 
            pos_motif1 = np.random.randint(START, seq_len - len(MOTIF2) - dist)
        pos_motif2 = pos_motif1 + dist

        sequence = sequence[:pos_motif1] + MOTIF1 + sequence[pos_motif1 + len(MOTIF1):]
        sequence = sequence[:pos_motif2] + MOTIF2 + sequence[pos_motif2 + len(MOTIF2):]
        
        ls_seq.append(sequence)
        ls_loc.append(loc[i])
        ls_dist.append(dist)
    
    if act_func == 'v-shape':
        score = v_shape(np.asarray(ls_dist))
    elif act_func == 'exp-fit-1':
        score = exp_fit_1(np.asarray(ls_dist))
    elif act_func == 'exp-fit-2':
        score = exp_fit_2(np.asarray(ls_dist))
    elif act_func == 'mound':
        score = mound_activation(np.asarray(ls_dist),
                                d0=40,
                                k=0.2,
                                decay_rate=0.002,
                                height=11,
                                x_offset=START)
    elif act_func == 'bowl':
        score = bowl_activation(np.asarray(ls_dist), 
                                y_max0=10,
                                y_max1=10,
                                bowl_depth=2,
                                x_offset=START)
    elif act_func == 'periodicity':
        score = periodicity_activation(np.asarray(ls_dist),
                                        period=10.4,
                                        amplitude=3,
                                        ampl_decay_rate=0.03,
                                        y_offset=8,
                                        offset_decay_rate=0.04,
                                        x_offset=START)
    elif act_func == 'cliff':
        score = cliff_activation(np.asarray(ls_dist),
                                y_max=10,
                                y_min=4,
                                x_intersection=25,
                                y_intersection=6,
                                x_min=10,
                                x_max=300) 
    elif act_func == 'hill':
        score = hill_activation(np.asarray(ls_dist),
                                y_max=10,
                                y_min=7,
                                x_min=10,
                                x_max=300) 
    else: 
        raise Exception('wrong activation function.')

    installed = pd.DataFrame({'location':ls_loc,'sequence':ls_seq, 'score':score})
    return installed
    
def write(df, outdir):
    fastaWriter = FastaWriter() 

    with gzip.open(outdir + ".fasta.gz", "wt") as OUT_FASTA, \
         gzip.open(outdir + "-counts.txt.gz", "wt") as OUT_COUNTS:

        print("DNA=1\tRNA=1\n", end="", file=OUT_COUNTS)

        for location, sequence, score in zip(df['location'], df['sequence'], df['score']):
            fastaWriter.addToFasta(location, sequence, OUT_FASTA)
            print(f"1\t{score}", file=OUT_COUNTS)


def main(indir, outdir, p, act_func):
    # select certain percentage of low-activation sequence from training set and install motifs

    # counts = pd.read_csv(indir+'/all-train-blacklisted-counts.txt.gz', compression='gzip', sep='\t', header=None, skiprows=1)
    # fasta = loadFasta(indir+'/all-train-blacklisted.fasta.gz')
    counts = pd.read_csv(indir+'/all-train-counts.txt.gz', compression='gzip', sep='\t', header=None, skiprows=1)
    fasta = loadFasta(indir+'/all-train.fasta.gz')
    n = int(len(counts)*p)

    counts.columns = ['DNA1','DNA2','DNA3','DNA4','DNA5','RNA1','RNA2','RNA3','RNA4']
    counts['DNA_aver'] = (counts['DNA1']+counts['DNA1']+counts['DNA1']+counts['DNA1']+counts['DNA1'])/5
    counts['RNA_aver'] = (counts['RNA1']+counts['RNA2']+counts['RNA3']+counts['RNA4'])/4
    counts['theta'] = np.log(counts['RNA_aver']/counts['DNA_aver'])
    counts['activation'] = abs(counts['theta'].median()-counts['theta'])
    counts.sort_values(by='activation',inplace=True)

    low_act_seq = fasta.loc[fasta.index.isin(counts[:n].index)]
    rdm_seq = fasta.loc[~fasta.index.isin(counts[:n].index)]
    low_act_seq = low_act_seq.reset_index(drop=True)
    rdm_seq = rdm_seq.reset_index(drop=True)

    # concat p% low-activation sequence and (1-p)% random sequence to be training set for downsampling
    low_act_seq = insert_motifs(low_act_seq, act_func)
    rdm_seq.loc[:, 'score'] = list(lognorm(s=1, scale=np.exp(0)).rvs(size=len(rdm_seq)))
    all_train = pd.concat([low_act_seq, rdm_seq], axis=0, ignore_index=True)

    # generate random score for downsampled testing and validation set
    test = loadFasta(indir+'/test.fasta.gz')
    test.loc[:,'score'] = lognorm(s=1, scale=np.exp(0)).rvs(size=len(test))
    validation = loadFasta(indir+'/validation.fasta.gz')
    validation.loc[:,'score'] = lognorm(s=1, scale=np.exp(0)).rvs(size=len(validation))

    write(test, outdir+'/test')
    write(validation, outdir+'/validation')
    write(all_train, outdir+'/all-train-'+str(round(p*100, 1)))


if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <indir> <outdir> <p:percent of positive samples> <act-func: v-shape, exp-fit-1,  exp-fit-2, mound, bowl, periodicity, cliff, hill>\n")
(indir, outdir, p, act_func)=sys.argv[1:]
main(indir, outdir, float(p)/100, act_func)

# python /hpc/group/igvf/A549/GR-AP1/simulated-seq/prepare-input.py /hpc/group/igvf/A549/GR-AP1/simulated-seq/data/Dex-200 /hpc/group/igvf/A549/GR-AP1/simulated-seq/data/mound 0.5 mound
# python /hpc/group/igvf/A549/GR-AP1/simulated-seq/prepare-input.py /hpc/group/igvf/A549/GR-AP1/simulated-seq/data/Dex-200 /hpc/group/igvf/A549/GR-AP1/simulated-seq/data/cliff 0.5 cliff
# python /hpc/group/igvf/A549/GR-AP1/simulated-seq/prepare-input.py /hpc/group/igvf/A549/GR-AP1/simulated-seq/data/Dex-200 /hpc/group/igvf/A549/GR-AP1/simulated-seq/data/bowl 0.5 bowl
# python /hpc/group/igvf/A549/GR-AP1/simulated-seq/prepare-input.py /hpc/group/igvf/A549/GR-AP1/simulated-seq/data/Dex-200 /hpc/group/igvf/A549/GR-AP1/simulated-seq/data/periodicity 0.5 periodicity