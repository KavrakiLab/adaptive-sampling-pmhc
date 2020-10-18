import numpy as np
import os
import sys
import mdtraj as md
import pyemma

import pyemma.coordinates as coor
#import matplotlib
#matplotlib.use('Agg')
#from matplotlib.pyplot import *
import glob
from pyemma import config 
import time
from pyemma import plots
import math
from subprocess import call
from subprocess import check_output

print("pyemma version:", pyemma.__version__) 


def get_raw_features(traj_name, featurizer):

    inp = coor.source(traj_name, featurizer, chunksize=5000)
    raw_Y = inp.get_output()
    return raw_Y

def get_pMHC_featurizer(feat_type, top):

    featurizer = coor.featurizer(top)

    peptide_residues = []
    system_residues = np.arange(top.n_residues)
    for resi in system_residues:
        if len(top.top.select("chainid == 1 and resi == " + str(resi))) > 0: peptide_residues.append(resi)

    if feat_type == 'pep_to_MHC':
        residue_pairs = []
        for peptide_residue in peptide_residues:
            for residue in system_residues:
                if peptide_residue == residue: continue
                residue_pairs.append([peptide_residue, residue])
        featurizer.add_residue_mindist(residue_pairs=np.array(residue_pairs), scheme='closest-heavy')

    else:
        print("Featurizer type not recognized")
        sys.exit(0)

    print("Number of atoms:", top.n_atoms)
    print("Number of residues:", top.n_residues)
    print("Number of features:", featurizer.dimension())

    return featurizer


def main():


    traj_folder = sys.argv[1]
    feature_type = 'pep_to_MHC'

    print(traj_folder, feature_type)

    traj_feature_prefix = traj_folder + "/" + feature_type

    traj = traj_folder + "/output.dcd"
    topfile = glob.glob(traj_folder + "/aln*.pdb")[0]
    top = md.load(topfile)
    featurizer = get_pMHC_featurizer(feature_type, top)
    Y = get_raw_features(traj, featurizer)
    call(["mkdir -p " + traj_feature_prefix], shell=True)
    np.save(traj_feature_prefix + "/Y.npy", Y[0])
    
    # should be done only once, since all trajectories have the same features
    if not os.path.exists("feats_des.npz"): np.savez_compressed("feats_des.npz", feat=featurizer.describe()) 

if __name__ == "__main__":
    main()

