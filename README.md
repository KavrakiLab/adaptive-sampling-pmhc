# How to build a MSM for pMHCs

## 1) System preparation (`system_prep/`)

### a) If a crystal structure of the system exists, use `get_pMHC_pdb.py` from `APE-Gen` to download PDB

```
python get_pMHC_pdb.py 3I6L
```

### b) The `mutate.py` script can be used to do simple point mutations on structures

```
pymol -qc mutate.py 3I6L.pdb C/4/ PRO D4P.pdb
```

### c) If no structure exists, you can use a docking program like `APE-Gen` to model a structure of your pMHC. For example,

```
python APE_Gen.py QFKPNVILL HLA-A*24:02
```

### d) Truncate and minimize structure as preparation for equilibration

```
python truncate_and_minimize.py D4P.pdb
rm temp.pdb fix.pdb
mv eq.pdb D4P_pre_eq.pdb
```

Temporary files `temp.pdb` and `fix.pdb` are created and can be safely removed. The main output is `eq.pdb` which holds the truncated and minimized conformation.

### e) Equilibration

Note that if the system you are running consists of a known unstable peptide, this step may end up detaching the peptide. 
Check the PDB after or use a shorter equilibration time.

```
python -u run_md_equilibration.py D4P_pre_eq.pdb output.dcd > eq.log &
```

### f) Get last frame for the umbrella sampling stage

Note that the files should be prepended with `aln-` which is used throughout other scripts to identify the topology files for a given trajectory.

```
mdconvert -t D4P_pre_eq.pdb -i -1 -o aln-D4P-start.pdb output.dcd
```

## Important note on file organization

There is a specific file structure that the scripts assume. 
Each trajectory is saved into its own folder, ordered by the folder name. 
Folders begin with at `0000/` and are padded with zeros until the folder name is 4 characters long.
Inside each folder, there should be one PDB file that begins with `aln-*` which represents the topology file, and a `output.dcd` trajectory file.
If the trajectory originates from umbrella sampling, there should be a `us_info.npz` file, which contains the equilibrium distance and the force constant of the umbrella.
Note that the existence of an npz file is required by the scripts in order to distinguish between biased (from umbrella sampling) and unbiased (from adaptive sampling) trajectories.

All of the main analysis is found in `msm_analysis.ipynb`. For running the script remotely, run `ssh -N -L localhost:8888:localhost:8888 COMPUTER` on the local machine, and `jupyter notebook --no-browser` on the remote machine.


## 2) Umbrella sampling stage (`us/`)

A set of umbrella sampling simulations are run to explore the major conformational states of the pMHC.
The reaction coordinate used requires that the pMHC is aligned such that the longest axis of the system (roughly from the binding site to the alpha_3 portion).
After doing so, the z coordinate is approximately aligned with the vector normal of the beta sheet floor, which characterizes the unbinding direction that a peptide must take to escape the binding site.
The reaction coordinate is measured as the z-coordinate distance between the center of mass between the beta sheet floor and the center of mass of the peptide.
In the paper, we refer to this as the z-dist value.
For the exact implementation, refer to `run_md_us.py`

We run umbrella sampling simulations centered from 1.0 nm (representing the native state z-dist value) to 3.0 nm (representing z-dist value of the unbound state) in increments of 0.1 nm using a force constant of 100 kJ/mol/nm^2.
From experience, simulations centered from the 2.0nm to 3.0nm range have the trajectories which contain the most binding and unbinding events, which the DTRAM reweighting algorithm needs to accurate approximate the unbiased MSM.
If the system is particularly unstable, a looser force constant may be required (like 10 kJ/mol/nm^2) to get smoother unbinding trajectories.

### a) Setting up simulations

HPC with GPUs is needed for this step.
We provide an example SLURM script for our runs.
With the file organization in mind, use the `setup_us_sim.py` to setup new folders for umbrella sampling simulations.

```
python setup_us_sim.py FOLDERNAME JOBNAME PDBNAME EQ_VAL US_FORCE_CONSTANT
```

This script creates a folder and copies over the file to run the OpenMM-based simulations, the SLURM script used to run on the computing cluster, and the starting conformation from the previous equilibration step (must have the PDB in this folder).
We also provide `setup_all_us.py` which contains the recommended initial set of umbrella sampling simulations to run.

### b) Run simulations

To submit a batch from `0000` to `0022` at once, it should be something like,

```
for i in $(seq -f "%04g" 0 22); do cd $i; sbatch run_us.sbatch; cd ..; done
```

Open the `output.csv` to check estimates of how long the simulation will take to finish. Check the queue with `squeue -u USERID`

### c) Visualizing trajectories

#### Watching movies

Use the cells under the `Visualize Trajectory` section to use nglview to visualize a particular simulation.

#### Getting basic statistics

All the cells up to the `Compute MSM` section contain intermediate statistics that may be helpful to visualize the sampling.


## 3) Adaptive sampling stage (`as/`)

Batches of unbiased simulations are run at a time. 
The selection of the starting conformations is based on the distribution of the trajectories in a 2D TICA space.
Conformations in the less sampled regions of the TICA space are more likely to be chosen as starting conformations.

### Initialization

The first round of simulations is recommended to be started from the native state.
From our experience, a batch of about 20 trajectories should suffice native state sampling.
Use the following script to create the trajectories (must have the native conformation in the same folder)

```
python setup_sim.py FOLDERNAME JOBNAME PDBNAME
```

### Choosing restarting conformations

Run each cell of the script up until the creation of the MSM.
Along the way, the number of connected microstates is computed, which can be used as an initial measure for convergence.
Additionally, there will be a cell which randomly chooses restarting conformations and plots on the TICA space.
Finally, there will be a cell which generates new folders with the restarting conformations extracted.
Note the starting folder for this round of simulations must be specified.
Batches of 20 simulations are run at a time.

## 4) Building the MSM

Run each cell of the script in the `Compute MSM` section.
Along the way, various plots will be generated related to the model, including for example, contact probabilities.
The subsequent section called `Network Analysis` will plot the percentage of unbinding trajectories going from the native state to the unbound state.

## 5) Mutational Analysis

In order for this section to be run, the energies of the conformations must be recomputed upon alanine mutation.
Subsample the trajectories to be every 1ns (stride=25) and follow the example script in `mut/array.sbatch`.

## 6) Bootstrapping

The last part of the notebook relates to incorporating bootstrapping.
Run `bootstrap_us.py` before executing those cells.








