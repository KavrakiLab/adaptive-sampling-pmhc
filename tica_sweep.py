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

#global_traj_indices = np.concatenate( (np.arange(386), np.arange(387,542), np.arange(543,646), np.arange(647,723)) )#np.arange(723) 542, 646
# 542 and 646 have errors from the trajectory saving process
# 386 produces an error in the TICA analysis
if mode == "wild": traj_to_exclude = [71,89,125,121] # 71,89,125 sink - edge moves to peptide, 121 totally unbound
elif mode == "D4A": traj_to_exclude = [75,95,96,97,98,99,100] # D4A weird starting states
elif mode == "D4P": traj_to_exclude = [50, 85,102, 86,99, 90, 327] 
unbound_trajs = [] #[e for e in range(707, 723)
global_traj_indices = np.array([gi for gi in np.arange(num_trajs) if gi not in traj_to_exclude])
local_traj_indices = [li for li, gi in enumerate(global_traj_indices)]

bound_trajs = [gi for gi in global_traj_indices if gi not in unbound_trajs]
local_bound_trajs = [li for li in local_traj_indices if global_traj_indices[li] not in unbound_trajs]
local_unbound_trajs = [li for li in local_traj_indices if global_traj_indices[li] in unbound_trajs]


traj_filenames = [str(i).zfill(4) for i in global_traj_indices]
inp_str = []
for traj_folder in traj_filenames:
    input_prefix = traj_folder + "/pep_to_MHC"
    inp_i = input_prefix + "/Y.npy"
    inp_str.append(inp_i)

inp = pyemma.coordinates.source(inp_str)

TICA_lagtimes = 25*np.array([1,10,25,50,75,100,250,500])
if os.path.exists("tica_timescales.npz"):
    f = np.load("tica_timescales.npz")
    lags = f["lags"]
    all_tica_timescales = f["timescales"]
else:
    all_tica_timescales = []
    for l in TICA_lagtimes:
        print("lag:",l)
        tica_obj_temp = pyemma.coordinates.tica(inp, lag=l, commute_map=True, kinetic_map=False, dim=50)
        tica_timescales = tica_obj_temp.timescales
        all_tica_timescales.append(tica_timescales)
    np.savez_compressed("tica_timescales.npz", lags=TICA_lagtimes, timescales=all_tica_timescales)


