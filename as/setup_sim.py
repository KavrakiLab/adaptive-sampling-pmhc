from subprocess import call
import sys
import os

foldername = sys.argv[1]
jobname = sys.argv[2]
pdbname = sys.argv[3]

if os.path.exists(foldername):
    print("Folder exists ... exiting.")
    sys.exit(0)

call(["mkdir " + foldername], shell=True)
os.chdir(foldername)

call(["cp ../run_md.py ."], shell=True)
call(["cp ../run.sbatch ."], shell=True)
call(["cp ../" + pdbname + " ."], shell=True)

call(["sed -i 's/JOBNAME/" + jobname + "/g' run.sbatch"], shell=True) 
call(["sed -i 's/PDBNAME/" + pdbname + "/g' run.sbatch"], shell=True)

