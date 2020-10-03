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

# specific to this file
def verifyParam(param):
    # check if has peptide start/end and system start/end
    return True

def get_raw_features(traj_name, featurizer):

    inp = coor.source(traj_name, featurizer, chunksize=5000)
    raw_Y = inp.get_output()
    return raw_Y

# TODO: try implementing a featurizer that looks at the interface!! (using electrostatic/SASA values?)
def get_pMHC_featurizer(feat_type, top): #, peptide_residues, system_residues):

    featurizer = coor.featurizer(top)

    peptide_residues = []
    system_residues = np.arange(top.n_residues)
    #new_system_residues = []
    for resi in system_residues:
        if len(top.top.select("chainid == 1 and resi == " + str(resi))) > 0: peptide_residues.append(resi)
        #elif len(top.top.select("chainid == 0 and resi == " + str(resi))) > 0 and (resi < 45 or (resi >= 95 and resi <= 120) ): new_system_residues.append(resi)
        #elif len(top.top.select("chainid == 0 and resi == " + str(resi))) > 0: new_system_residues.append(resi)
    #system_residues = new_system_residues

    if feat_type == 'pep_to_MHC':
        residue_pairs = []
        for peptide_residue in peptide_residues:
            for residue in system_residues:
                if peptide_residue == residue: continue
                residue_pairs.append([peptide_residue, residue])
        featurizer.add_residue_mindist(residue_pairs=np.array(residue_pairs), scheme='closest-heavy')

    elif feat_type == 'pep_to_MHC_ca':
        residue_pairs = []
        for peptide_residue in peptide_residues:
            for residue in system_residues:
                if peptide_residue == residue: continue
                residue_pairs.append([peptide_residue, residue])
        featurizer.add_residue_mindist(residue_pairs=np.array(residue_pairs), scheme='ca')

    elif feat_type == 'pep_bb_ca_torsions':
        resi_str = "resi " + str(peptide_residues[0]) + " to " + str(peptide_residues[-1])
        featurizer.add_backbone_torsions(selstr=resi_str, cossin=True)
        featurizer.add_sidechain_torsions(selstr=resi_str, cossin=True)

    elif feat_type == 'pep_bb_torsions':
        resi_str = "resi " + str(peptide_residues[0]) + " to " + str(peptide_residues[-1])
        featurizer.add_backbone_torsions(selstr=resi_str, cossin=True)

    elif feat_type == 'pep_bb_ca':
        resi_str = "resi " + str(peptide_residues[0]) + " to " + str(peptide_residues[-1])
        bb_ca_str = resi_str + " and backbone and name == 'CA'"
        bb_ca_indices = top.top.select(bb_ca_str)
        featurizer.add_distances(indices=featurizer.pairs(bb_ca_indices))

    elif feat_type == 'pep_bb_torsions_and_ca':
        resi_str = "resi " + str(peptide_residues[0]) + " to " + str(peptide_residues[-1])
        featurizer.add_backbone_torsions(selstr=resi_str, cossin=True)
        bb_ca_str = resi_str + " and backbone and name == 'CA'"
        bb_ca_indices = top.top.select(bb_ca_str)
        featurizer.add_distances(indices=featurizer.pairs(bb_ca_indices))

    elif feat_type == 'sasa':
        featurizer.add_custom_func(get_sasa, len(system_residues))

    else:
        print("Featurizer type not recognized")
        sys.exit(0)

    print("Number of atoms:", top.n_atoms)
    print("Number of residues:", top.n_residues)
    print("Number of features:", featurizer.dimension())

    return featurizer

def get_sasa(traj):
    return md.shrake_rupley(traj, mode='residue')

def one_over_d(traj):
    dd = 1./(md.compute_contacts(traj,contacts='all',scheme='closest-heavy')[0])
    return dd

def main():


    traj_folder = sys.argv[1]
    feature_type = 'pep_to_MHC' #sys.argv[2]

    print(traj_folder, feature_type)

    #traj_folder = str(traj_folder).zfill(4)
    traj_feature_prefix = traj_folder + "/" + feature_type
    #alreadyFeaturized = os.path.exists(traj_feature_prefix + "/raw_features.npz")
    #if not alreadyFeaturized:

    traj = traj_folder + "/output.dcd"
    topfile = glob.glob(traj_folder + "/*.pdb")[0] #traj_folder + "/eq.pdb"
    top = md.load(topfile)
    featurizer = get_pMHC_featurizer(feature_type, top) #, peptide_residues, all_residues)
    Y = get_raw_features(traj, featurizer)
    call(["mkdir -p " + traj_feature_prefix], shell=True)
    #np.savez_compressed(feature_type + "/feats_des.npz", feat=featurizer.describe()) 
    np.save(traj_feature_prefix + "/Y.npy", Y[0])
    
    if not os.path.exists("feats_des.npz"): np.savez_compressed("feats_des.npz", feat=featurizer.describe()) 

if __name__ == "__main__":
    main()

