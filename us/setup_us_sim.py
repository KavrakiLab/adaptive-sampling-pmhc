from subprocess import call
import sys
import os

foldername = sys.argv[1]
jobname = sys.argv[2]
pdbname = sys.argv[3]
eq_val = sys.argv[4]
us_force_constant = sys.argv[5]

if os.path.exists(foldername):
    print("Folder exists ... exiting.")
    sys.exit(0)

call(["mkdir " + foldername], shell=True)
os.chdir(foldername)

call(["cp ../run_md_us.py ."], shell=True)
call(["cp ../run_us.sbatch ."], shell=True)
call(["cp ../aln*.pdb ."], shell=True)

call(["sed -i 's/JOBNAME/" + jobname + "/g' run_us.sbatch"], shell=True) 
call(["sed -i 's/PDBNAME/" + pdbname + "/g' run_us.sbatch"], shell=True)
call(["sed -i 's/EQ_VAL/" + eq_val + "/g' run_us.sbatch"], shell=True)
call(["sed -i 's/US_FORCE_CONSTANT/" + us_force_constant + "/g' run_us.sbatch"], shell=True)

