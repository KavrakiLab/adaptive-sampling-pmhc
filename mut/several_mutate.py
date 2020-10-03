import sys,os
from subprocess import call
import glob


first, last = int(sys.argv[1]), int(sys.argv[2])

for i in range(first, last+1):

	folder_name = str(i).zfill(4)

	#call(["tar xopf " + folder_name + ".tar && rm " + folder_name + ".tar"], shell=True)

	os.chdir(folder_name)

	topfile = glob.glob("aln-*.pdb")[0]
	call(["cp ../*.py ."], shell=True)
	call(["python mutate_and_rescore.py " + topfile], shell=True)

	os.chdir("..")


