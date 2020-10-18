from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import sys
import numpy as np
import time
import mdtraj as md
#from pdbfixer import PDBFixer
#import matplotlib.pyplot as plt

import glob


force_constant = 100 # force constant for restraining the beta sheet floor
num_steps = 250000000 # 4 fs timestep, 1 us production
printInterval = 10000 # at 4fs timestep is saving every (1/25) nanoseconds

def main():

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    dz0 = float(sys.argv[3])
    us_force_constant = float(sys.argv[4])
    temp = 300

    np.savez_compressed("us_info.npz", center=dz0, force_constant=us_force_constant)


    ref = input_filename    

    pdb = PDBFile(input_filename)

    reference_structure = md.load(ref)
    peptide_particles = reference_structure.top.select("chainid == 1")
    peptide_top = reference_structure.top.subset(peptide_particles)
    peptide_traj = md.Trajectory(reference_structure.xyz[:, peptide_particles, :], peptide_top)


    platform = Platform.getPlatformByName('OpenCL')

    forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')
    modeller = Modeller(pdb.topology, pdb.positions)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic, constraints=app.HBonds, hydrogenMass=4*amu) # changed from AllBonds

    # constrain alpha carbon atoms of beta sheet floor of MHC
    force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k", force_constant)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    protein_particles = md.load(ref).top.select("chainid != 1 and name == 'CA' and (resi < 45 or (resi >= 95 and resi <= 120))")

    particle_indices = []
    for protein_particle in protein_particles:
        particle_indices.append(force.addParticle(int(protein_particle), modeller.positions[protein_particle]) )
    system.addForce(force)

    # umbrella sampling force
    restraint_force = CustomCentroidBondForce(2, "k*(sqrt( (z1-z2)^2 ) - dz0)^2")
    restraint_force.addGlobalParameter("dz0", dz0)
    restraint_force.addPerBondParameter("k")
    
    peptide_particles = md.load(ref).top.select("chainid == 1")
    new_pep = [int(p) for p in peptide_particles]
    new_mhc = [int(p) for p in protein_particles]
    
    restraint_force.addGroup(new_pep)
    restraint_force.addGroup(new_mhc)
    restraint_force.addBond([int(0), int(1)], [us_force_constant])
    system.addForce(restraint_force)


    forces = system.getForces()
    i = 0
    for f in forces:
        f.setForceGroup(i)
        i = i + 1

    integrator = LangevinIntegrator(temp*kelvin, 1/picosecond, 0.004*picoseconds)
    simulation = Simulation(modeller.topology, system, integrator, platform)

    simulation.context.setPositions(modeller.positions)

    totalSimulationTime = num_steps
    simulation.reporters.append(app.StateDataReporter(stdout, printInterval, step=True, 
    potentialEnergy=True, temperature=True, progress=False, remainingTime=True, 
    speed=True, totalSteps=totalSimulationTime, separator='\t'))

    r = DCDReporter(output_filename, printInterval, enforcePeriodicBox=False)
    r.report(simulation, simulation.context.getState(getPositions=True, getEnergy=True)) # get starting state
    simulation.reporters.append(r)

    simulation.reporters.append(app.StateDataReporter(output_filename[:-4] + ".csv", printInterval, step=True, 
    potentialEnergy=True, temperature=True, progress=False, remainingTime=True, 
    speed=True, totalSteps=totalSimulationTime, separator='\t'))

    print("Running MD with k=" + str(force_constant) + " at T=" + str(temp) + " for " + str(num_steps) + " steps.")
    simulation.step(num_steps)



if __name__ == "__main__":
    main()



