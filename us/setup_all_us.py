from subprocess import call

all_centers = [1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,
               2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,
               1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,
               2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,
               2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,
               2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0]

us_force_constant = 100

for i in range(len(all_centers)):
    c = all_centers[i]
    inputs = [str(i).zfill(4), str(i).zfill(4), "aln-D4P-start.pdb", str(c), str(us_force_constant)]
    call(["python setup_us_sim.py " + " ".join(inputs)], shell=True)


