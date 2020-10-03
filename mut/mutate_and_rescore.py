import numpy as np
import mdtraj as md
from subprocess import call
from subprocess import check_output
import sys

#call(["mkdir wild"], shell=True)

topfile = sys.argv[1]

#call(["mdconvert output.dcd -s 25 -t " + topfile + " -o output_every1ns.dcd"], shell=True)

f = md.load("output_every1ns_fix.dcd", top=topfile)

f_out = open("energies.txt", 'w')
for i, conf in enumerate(f):

    conf_name = str(i).zfill(5) + ".pdb"
    conf.save_pdb(conf_name)
    for mutant in ["wild", "Q1", "F2", "K3", "D4", "N5", "V6", "I7", "L8", "L9", "Q155", "Y159"]:
        if mutant not in ["Q155", "Y159"]: continue
        new_conf_name = mutant + "-" + conf_name
        if mutant == "wild": call(["cp " + conf_name + " " + new_conf_name], shell=True)
        elif mutant in ["V6","I7"]: call(["pymol -qc mutate_V6I7.py " + conf_name + " B/" + mutant[1] + "/ ALA " + new_conf_name], shell=True)
        #else: call(["pymol -qc mutate.py " + conf_name + " B/" + mutant[1] + "/ ALA " + new_conf_name], shell=True)
        else: call(["pymol -qc mutate.py " + conf_name + " A/" + mutant[1:] + "/ ALA " + new_conf_name], shell=True)        

        if mutant in ["Q1"]: 
            call(["sed '/B   /d' " + new_conf_name + " | sed '/TER/d' | sed '/END/d' > temp1.pdb"], shell=True)
            call(["grep 'H2  GLN B' " + conf_name + " | sed 's/H2  GLN B/H2  ALA B/g' >> temp2.pdb"], shell=True)
            call(["grep 'H3  GLN B' " + conf_name + " | sed 's/H3  GLN B/H3  ALA B/g' >> temp2.pdb"], shell=True)
            call(["grep 'B   ' " + new_conf_name + " > temp3.pdb"], shell=True)
            call(["cat temp1.pdb temp2.pdb temp3.pdb > " + new_conf_name], shell=True)
            call(["rm temp1.pdb temp2.pdb temp3.pdb"], shell=True)
            #call(["pdbfixer " + new_conf_name + " --output=" + new_conf_name], shell=True)
        elif mutant in ["L9"]:
            call(["sed '/TER/d' " + new_conf_name + " | sed '/END/d' > temp1.pdb"], shell=True)
            call(["grep 'OXT LEU B' " + conf_name + " | sed 's/LEU/ALA/g' >> temp1.pdb"], shell=True)
            call(["cp temp1.pdb " + new_conf_name], shell=True)
            call(["rm temp1.pdb"], shell=True)
            #call(["pdbfixer " + new_conf_name + " --output=" + new_conf_name], shell=True)
        energy = float(check_output(["python energy.py " + new_conf_name + " | awk '{ print $1 }'"], shell=True))
        print(mutant, i, energy)
        f_out.write(mutant + " " + str(energy) + "\n")
    #break

f_out.close()

