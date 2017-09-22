import pandas as pd

def readThreeColumnTruth(fn, suffix=""):
    df = pd.read_csv(fn, sep=' ', skiprows=1,
                     names=['Name', 'Gene{}'.format(suffix),
                            'TPM{}'.format(suffix)], engine='c')
    df.set_index("Name", inplace=True)
    pd.to_numeric(df["TPM{}".format(suffix)], errors='ignore')
    return df

def readRSEMTruth(fn, suffix=""):
    df = pd.read_csv(fn, sep='\t', skiprows=1,
                     names=['Name', 'Gene{}'.format(suffix),
                            'Length{}'.format(suffix),
                            'EffectiveLength{}'.format(suffix),
                            'NumReads{}'.format(suffix),
                            'TPM{}'.format(suffix),
                            'FPKM{}'.format(suffix),
                            'IsoPct{}'.format(suffix)], engine='c').set_index('Name')
    for col in ["TPM", "Length", "EffectiveLength", "NumReads"]:
        pd.to_numeric(df["{}{}".format(col, suffix)], errors='ignore')
    return df

def readRSEM(fn, suffix=""):
    df = pd.read_csv(fn, sep='\t', skiprows=1,
                     names=['Name', 'Gene{}'.format(suffix),
                            'Length{}'.format(suffix),
                            'EffectiveLength{}'.format(suffix),
                            'NumReads{}'.format(suffix),
                            'TPM{}'.format(suffix),
                            'FPKM{}'.format(suffix),
                            'IsoPct{}'.format(suffix)], engine='c').set_index('Name')
    for col in ["TPM", "Length", "EffectiveLength", "NumReads"]:
        pd.to_numeric(df["{}{}".format(col, suffix)], errors='ignore')
    return df

def readStringTie(fn, suffix=""):
    """
    Not yet tested
    """
    df = pd.read_csv(fn, sep="\t", skiprows=1,
                     names=["tid{}".format(suffix),
                            "chr{}".format(suffix),
                            "strand{}".format(suffix),
                            "start{}".format(suffix),
                            "end{}".format(suffix),
                            "Name",
                            "num_exons{}".format(suffix),
                            "Length{}".format(suffix),
                            "gene_id{}".format(suffix),
                            "gene_name{}".format(suffix),
                            "cov{}".format(suffix),
                            "FPKM{}".format(suffix)])
    df.set_index('Name', inplace=True)
    pd.to_numeric(df)
    return df

def readExpress(fn, suffix=""):
    df = pd.read_csv(fn, sep="\t", skiprows=1,
                     names=["bundle_id{}".format(suffix),
                            "Name",
                            "Length{}".format(suffix),
                            "EffectiveLength{}".format(suffix),
                            "tot_counts{}".format(suffix),
                            "uniq_counts{}".format(suffix),
                            "NumReads{}".format(suffix),
                            "NumReadsEffective{}".format(suffix),
                            "ambig_distr_alpha{}".format(suffix),
                            "ambig_distr_beta{}".format(suffix),
                            "fpkm{}".format(suffix),
                            "fpkm_conf_low{}".format(suffix),
                            "fpkm_conf_high{}".format(suffix),
                            "solvable{}".format(suffix),
                            "TPM{}".format(suffix)]).set_index('Name')
    for col in ["TPM", "Length", "EffectiveLength", "NumReads"]:
        pd.to_numeric(df["{}{}".format(col, suffix)], errors='ignore')
    return df

def readKallistoBoot(fn, suffix=""):
    import h5py
    import os
    import numpy as np
    h5file = os.path.sep.join([fn,"abundance.h5"])
    f = h5py.File(h5file)
    names = map(str, f['aux']['ids'].value)
    nboot = len(f['bootstrap'])
    boots = []
    for i in xrange(nboot):
        boots.append(f['bootstrap']['bs{}'.format(i)].value)
    y = pd.DataFrame(data=boots, dtype=np.float64).T
    y = y.assign(Name=names).set_index('Name').sort_index()
    y = y.apply(np.sort, axis=1)
    return y

def readSalmonBoot(fn, suffix=""):
    import os
    import gzip
    import pandas as pd
    import numpy as np
    import struct
    import json
    auxDir = "aux"
    # Check for a custom auxDir
    with open(os.path.sep.join([fn, "cmd_info.json"])) as cmdFile:
        dat = json.load(cmdFile)
        if 'auxDir' in dat:
            auxDir = dat['auxDir']

    bootstrapFile = os.path.sep.join([fn, auxDir, "bootstrap", "bootstraps.gz"])
    nameFile = os.path.sep.join([fn, auxDir, "bootstrap", "names.tsv.gz"])
    bootstrapFile = os.path.sep.join([fn, auxDir, "bootstrap", "bootstraps.gz"])
    nameFile = os.path.sep.join([fn, auxDir, "bootstrap", "names.tsv.gz"])
    if not os.path.isfile(bootstrapFile):
       print("The required bootstrap file {} doesn't appear to exist".format(bootstrapFile)) 
       sys.exit(1)
    if not os.path.isfile(nameFile):
       print("The required transcript name file {} doesn't appear to exist".format(nameFile)) 
       sys.exit(1)

    txpNames = None
    with gzip.open(nameFile) as nf:
        txpNames = nf.read().strip().split('\t')
    
    ntxp = len(txpNames)
    print("Expecting bootstrap info for {} transcripts".format(ntxp))
    
    with open(os.path.sep.join([fn, auxDir, "meta_info.json"])) as fh:
        meta_info = json.load(fh)
        
    stype = None
    if meta_info['samp_type'] == 'gibbs':
        s = struct.Struct('@' + 'd' * ntxp)
        stype = 'g'
    elif meta_info['samp_type'] == 'bootstrap':
        s = struct.Struct('@' + 'd' * ntxp)
        stype = 'b'
    else:
        print("Unknown sampling method: {}".format(meta_info['samp_type']))
        sys.exit(1)
        
    numBoot = 0
    samps = []
    convert = float
    # Now, iterate over the bootstrap samples and write each
    with gzip.open(bootstrapFile) as bf:
        while True:
            try:
                x = s.unpack_from(bf.read(s.size))
                xs = map(convert, x)
                samps.append(xs)
                numBoot += 1
            except:
                print("read all posterior values")
                break

    print("wrote {} bootstrap samples".format(numBoot))
    print("converted bootstraps successfully.")
    y = pd.DataFrame(data=samps, dtype=np.float64).T
    y = y.assign(Name=txpNames).set_index('Name').sort_index()
    y = y.apply(np.sort, axis=1)
    return y

def readSailfish(fn, suffix=""):
    df = pd.read_table(fn, engine='c').set_index('Name')
    df.columns = [ "{}{}".format(cn, suffix) for cn in df.columns.tolist()]
    for col in ["TPM", "Length", "EffectiveLength", "NumReads"]:
        pd.to_numeric(df["{}{}".format(col, suffix)], errors='ignore')
    return df

def readSalmon(fn, suffix=""):
    return readSailfish(fn, suffix)

def readSailfishDeprecated(fn, suffix=""):
    df = pd.read_csv(fn, sep='\t', comment='#',
                     names=['Name',
                            'Length{}'.format(suffix),
                            'TPM{}'.format(suffix),
                            'RPKM{}'.format(suffix),
                            'KPKM{}'.format(suffix),
                            'NumKmers{}'.format(suffix),
                            'NumReads{}'.format(suffix)])
    df.dropna(how='all', inplace=True)
    df.convert_objects(convert_numeric=True)
    df.set_index('Name', inplace=True)
    return df

def readKallisto(fn, suffix=""):
    df = pd.read_csv(fn, sep='\t', skiprows=1,
                     names=['Name',
                            'Length{}'.format(suffix),
                            'EffectiveLength{}'.format(suffix),
                            'NumReads{}'.format(suffix),
                            'TPM{}'.format(suffix)], engine='c').set_index('Name')
    for col in ["TPM", "Length", "EffectiveLength", "NumReads"]:
        pd.to_numeric(df["{}{}".format(col, suffix)], errors='ignore')
    return df

def readProFile(fn, suffix=""):
    df = pd.read_csv(fn, sep='\t',
                     names=['Locus{}'.format(suffix),
                            'Name',
                            'Coding{}'.format(suffix),
                            'Length{}'.format(suffix),
                            'ExpFrac{}'.format(suffix),
                            'ExpNum{}'.format(suffix),
                            'LibFrac{}'.format(suffix),
                            'LibNum{}'.format(suffix),
                            'SeqFrac{}'.format(suffix),
                            'SeqNum{}'.format(suffix),
                            'CovFrac{}'.format(suffix),
                            'ChiSquare{}'.format(suffix),
                            'CV{}'.format(suffix)]).set_index('Name')
    for col in ["Length", "ExpFrac", "ExpNum", "LibFrac", "LibNum", "SeqFrac", "SeqNum", "CovFrac", "ChiSquare", "CV"]:
        pd.to_numeric(df["{}{}".format(col, suffix)], errors='ignore')
    return df

def readResFile(fn, suffix=""):
    df = pd.read_csv(fn, sep='\t',
                     names=['Name',
                            'Length{}'.format(suffix),
                            'Abund{}'.format(suffix)]).set_index('Name')
    for col in ["Length", "Abund"]:
        pd.to_numeric(df["{}{}".format(col, suffix)], errors='ignore')
    return df
