import glob
from subprocess import call
import os
import sys


filename = sys.argv[1]

# truncate
f = open(filename, 'r')
f_out = open("temp.pdb", 'w')

for line in f:
	line_arr = line.split()
	if line_arr[0] == "ATOM": 
		if line_arr[4] == "A":
			if int(line_arr[5]) > 180:
				continue
		elif line_arr[4] == "B": continue
	elif line_arr[0] == "TER": continue

	f_out.write(line)

f.close()
f_out.close()

# prep for minimize
call(["pdbfixer temp.pdb --output=fix.pdb"], shell=True)

# minimize
call(["python minimize.py fix.pdb eq.pdb"], shell=True)

