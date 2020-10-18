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


force_constant = 100
num_steps = 25000000 # 4 fs timestep, 100 ns equilibration
printInterval = 10000

def main():

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    temp = 300

    ref = input_filename    

    pdb = PDBFile(input_filename)

    reference_structure = md.load(ref)
    peptide_particles = reference_structure.top.select("chainid == 1")
    peptide_top = reference_structure.top.subset(peptide_particles)
    peptide_traj = md.Trajectory(reference_structure.xyz[:, peptide_particles, :], peptide_top)

    platform = Platform.getPlatformByName('OpenCL')

    forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')
    modeller = Modeller(pdb.topology, pdb.positions)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic, constraints=app.HBonds, hydrogenMass=4*amu)

    # constrain alpha carbon atoms of MHC
    force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k", force_constant)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")


    protein_particles = md.load(ref).top.select("chainid != 1 and name == 'CA'")
    particle_indices = []
    for protein_particle in protein_particles:
        particle_indices.append(force.addParticle(int(protein_particle), modeller.positions[protein_particle]) )
    system.addForce(force)

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


def printForces(simulation):

    for i in range(simulation.system.getNumForces()):
        f_name = simulation.system.getForce(i).__class__.__name__
        s0 = simulation.context.getState(getEnergy=True, groups=2**i)
        print(f_name + ": " + str(s0.getPotentialEnergy()))

    print("total:", simulation.context.getState(getEnergy=True).getPotentialEnergy())


if __name__ == "__main__":
    main()



