import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pyemma
import nglview
import glob
import matplotlib as mpl
import msmtools
import os
from pyemma import config
import sys

config.mute = True
config.show_progress_bars = False

mode = "D4P"
num_trajs = 0
for i in range(1000):
    if os.path.exists(str(i).zfill(4) + "/output.dcd"): num_trajs += 1
    else: break
print(num_trajs)  

traj_to_exclude = []
unbound_trajs = [] #[e for e in range(707, 723)
global_traj_indices = np.array([gi for gi in np.arange(num_trajs) if gi not in traj_to_exclude])
local_traj_indices = [li for li, gi in enumerate(global_traj_indices)]

bound_trajs = [gi for gi in global_traj_indices if gi not in unbound_trajs]
local_bound_trajs = [li for li in local_traj_indices if global_traj_indices[li] not in unbound_trajs]
local_unbound_trajs = [li for li in local_traj_indices if global_traj_indices[li] in unbound_trajs]

f = np.load("discretization2.npz", allow_pickle=True)
dtrajs = list(f['dtrajs'])
cc_all = f['cc_all']
index_clusters = f['index_clusters']
n_clusters = len(index_clusters)
f.close()

kT = 2.479
umbrella_sampling_g_indices = []
umbrella_sampling_l_indices = []
us_trajs = []
md_trajs = []
us_centers = []
us_force_constants = []
unbiased_g_trajs = []
unbiased_l_trajs = []

for progress, i in enumerate(global_traj_indices):
    if os.path.exists(str(i).zfill(4) + "/dCOM.npz"):
        dCOM = np.load(str(i).zfill(4) + "/dCOM.npz", allow_pickle=True)["dCOM"]
    else:
        print(progress)
        file_prefix = str(i).zfill(4) + "/"
        topfile = glob.glob(file_prefix + "aln*.pdb")[0]
        f = md.load(file_prefix + "output_every1ns_fix.dcd", top=topfile)
        peptide = f.top.select("chainid == 1")
        mhc = f.top.select("chainid != 1 and name == 'CA' and (resi < 45 or (resi >= 95 and resi <= 120))")

        peptide_frame = f.atom_slice(atom_indices=peptide, inplace=False)
        mhc_frame = f.atom_slice(atom_indices=mhc, inplace=False)

        pep_com = md.compute_center_of_mass(peptide_frame)
        mhc_com = md.compute_center_of_mass(mhc_frame)

        dCOM = np.abs(pep_com - mhc_com)[:,2]
        np.savez_compressed(str(i).zfill(4) + "/dCOM.npz", dCOM=dCOM)
        
    if os.path.exists(str(i).zfill(4) + "/us_info.npz"):
        us_info_file = np.load(str(i).zfill(4) + "/us_info.npz", allow_pickle=True)
        us_centers.append(us_info_file["center"])
        us_force_constants.append(us_info_file["force_constant"] / kT)
        us_trajs.append(dCOM)
        umbrella_sampling_g_indices.append(i)
        umbrella_sampling_l_indices.append(progress)
        
    else:
        md_trajs.append(dCOM)
        unbiased_g_trajs.append(i)
        unbiased_l_trajs.append(progress)
        
        
#us_centers = np.array(us_centers)
#us_force_constants = np.array(us_force_constants)
print(len(us_trajs), len(md_trajs))
print(us_centers)
print(us_force_constants)

us_dtrajs = [dtrajs[i] for i in umbrella_sampling_l_indices]
md_dtrajs = [dtrajs[i] for i in unbiased_l_trajs]

indices_bootstrap = []
for i in range(len(md_dtrajs)):
    ind_b_i = [(i,j,0) for j in np.arange(len(md_dtrajs[i]))]
    indices_bootstrap.append(ind_b_i)
for i in range(len(us_dtrajs)):
    ind_b_i = [(i,j,1) for j in np.arange(len(us_dtrajs[i]))]
    indices_bootstrap.append(ind_b_i)

lags = 25*np.array([1,10,20,25,50,75,100,250,500])
lags = [lags[int(sys.argv[1])]]
for lag in lags:
    print(lag)
    for b_ind in range(100):
        new_ind = msmtools.estimation.bootstrap_trajectories(indices_bootstrap, -1)

        new_md_trajs = []
        new_md_dtrajs = []
        new_us_trajs = []
        new_us_dtrajs = []
        new_us_centers = []
        new_us_force_constants = []
        for t in new_ind:
            traj_index = t[0][0]
            beg = t[0][1]
            end = t[-1][1]
            traj_type = t[0][2]
            if traj_type == 0:
                new_md_trajs.append(md_trajs[traj_index][beg:end+1])
                new_md_dtrajs.append(md_dtrajs[traj_index][beg:end+1])
            elif traj_type == 1:
                new_us_trajs.append(us_trajs[traj_index][beg:end+1])
                new_us_dtrajs.append(us_dtrajs[traj_index][beg:end+1])
                new_us_centers.append(us_centers[traj_index])
                new_us_force_constants.append(us_force_constants[traj_index])

        #us = pyemma.thermo.estimate_umbrella_sampling(us_trajs=us_trajs, us_dtrajs=us_dtrajs, 
        #                                              us_centers=us_centers, us_force_constants=us_force_constants,
        #                                              md_trajs=new_md_trajs, md_dtrajs=new_md_dtrajs, estimator="dtram", lag=250*25, 
        #                                              width=None, save_convergence_info=1, maxiter=10000)#, maxerr=1E-6)
        us = pyemma.thermo.estimate_umbrella_sampling(us_trajs=new_us_trajs, us_dtrajs=new_us_dtrajs,
                                                      us_centers=new_us_centers, us_force_constants=new_us_force_constants,
                                                      md_trajs=new_md_trajs, md_dtrajs=new_md_dtrajs, estimator="dtram", lag=lag,
                                                      width=None, save_convergence_info=1, maxiter=10000)#, maxerr=1E-6)


        MSM = us.msm
        us.save("bootstrap/us_obj_"+str(lag)+"_"+str(b_ind),"us_obj",overwrite=True)
        MSM.save("bootstrap/msm_obj_"+str(lag)+"_"+str(b_ind),"msm_obj",overwrite=True)



