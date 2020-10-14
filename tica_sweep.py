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

# Because there is no bounding box in the implicit solvent simulations, 
# trajectories must be removed that unbind and travel far from the binding site
traj_to_exclude = []
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


