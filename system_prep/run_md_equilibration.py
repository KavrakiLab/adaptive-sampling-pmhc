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

    #ref = sys.argv[1]
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    temp = 300

    ref = input_filename    

    #print "Opening:", input_filename
    pdb = PDBFile(input_filename)

    reference_structure = md.load(ref)
    peptide_particles = reference_structure.top.select("chainid == 1")
    peptide_top = reference_structure.top.subset(peptide_particles)
    peptide_traj = md.Trajectory(reference_structure.xyz[:, peptide_particles, :], peptide_top)
    #reference_com = md.compute_center_of_mass(peptide_traj)[0]
    #print reference_com

    platform = Platform.getPlatformByName('OpenCL')

    forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')
    #forcefield = app.ForceField('amber03.xml', 'amber03_obc.xml')
    #forcefield = app.ForceField('amber10.xml', 'amber10_obc.xml')
    modeller = Modeller(pdb.topology, pdb.positions)
    #system = forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic, constraints=app.HBonds)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic, constraints=app.HBonds, hydrogenMass=4*amu) # changed from AllBonds

    force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k", force_constant)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")


    # constrain all MHC backbone
    protein_particles = md.load(ref).top.select("chainid != 1 and name == 'CA'") # alpha carbon of MHC
    #protein_particles = md.load(ref).top.select("backbone")

    # constrain only MHC backbone atoms involved in beta sheets
    """
    resi_str = ""
    dssp_assignment = md.compute_dssp(md.load(ref))[0]
    for resi, ss in enumerate(dssp_assignment):
        if ss == "E": resi_str += "resi " + str(resi) + " or "
    resi_str = resi_str[:-4]
    protein_particles = md.load(ref).top.select("chainid != 1 and name == 'CA' and (" + resi_str + ")")
    """
    #protein_particles = md.load(ref).top.select("chainid != 1 and name == 'CA' and (resi < 45 or (resi >= 95 and resi <= 120))")

    particle_indices = []
    for protein_particle in protein_particles:
        particle_indices.append(force.addParticle(int(protein_particle), modeller.positions[protein_particle]) )
    system.addForce(force)

    forces = system.getForces()
    i = 0
    for f in forces:
        f.setForceGroup(i)
        i = i + 1

    #integrator = LangevinIntegrator(temp*kelvin, 1/picosecond, 0.002*picoseconds)
    integrator = LangevinIntegrator(temp*kelvin, 1/picosecond, 0.004*picoseconds)
    simulation = Simulation(modeller.topology, system, integrator, platform)

    simulation.context.setPositions(modeller.positions)

    totalSimulationTime = num_steps
    simulation.reporters.append(app.StateDataReporter(stdout, printInterval, step=True, 
    potentialEnergy=True, temperature=True, progress=False, remainingTime=True, 
    speed=True, totalSteps=totalSimulationTime, separator='\t'))

    #r = md.reporters.HDF5Reporter(output_filename, printInterval) #, enforcePeriodicBox=False)
    r = DCDReporter(output_filename, printInterval, enforcePeriodicBox=False)
    r.report(simulation, simulation.context.getState(getPositions=True, getEnergy=True)) # get starting state
    simulation.reporters.append(r)

    simulation.reporters.append(app.StateDataReporter(output_filename[:-4] + ".csv", printInterval, step=True, 
    potentialEnergy=True, temperature=True, progress=False, remainingTime=True, 
    speed=True, totalSteps=totalSimulationTime, separator='\t'))
    
    """
    procedure_intervals = [(5000, 350*1000/2), (1000, 10*1000/2), (100, 10*1000/2), (10, 10*1000/2), (0, 120*1000/2)]
    #procedure_intervals = [(5000, 500), (1000, 100), (100, 100), (10, 100), (0, 200)]
    for interval in procedure_intervals:
        force_constant = interval[0]
        num_steps_i = interval[1]
        print "Running MD with k=" + str(force_constant) + " for " + str(num_steps_i) + " steps."
        force.setGlobalParameterDefaultValue(0, force_constant)
        force.updateParametersInContext(simulation.context)
        simulation.step(num_steps_i)
    """

    print("Running MD with k=" + str(force_constant) + " at T=" + str(temp) + " for " + str(num_steps) + " steps.")
    simulation.step(num_steps)

    # currentSimTime = 0
    #while currentSimTime < num_steps:

    #    simulation.step(1000)
    #    currentSimTime += 1000

        #cur_com = get_current_com(simulation, peptide_top, peptide_particles)
        #d_com = np.linalg.norm(reference_com - cur_com)
        #print d_com


def printForces(simulation):

    for i in range(simulation.system.getNumForces()):
        f_name = simulation.system.getForce(i).__class__.__name__
        s0 = simulation.context.getState(getEnergy=True, groups=2**i)
        print(f_name + ": " + str(s0.getPotentialEnergy()))

    print("total:", simulation.context.getState(getEnergy=True).getPotentialEnergy())


def get_current_com(simulation, top, peptide_particles):
    peptide_top = top
    positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
    new_positions = np.zeros((1, len(peptide_particles), 3))
    new_positions[0, :, :] = positions[peptide_particles, :]
    #print positions.shape, peptide_particles.shape, np.array([ positions[peptide_particles, :] ])
    peptide_traj = md.Trajectory(new_positions, peptide_top)
    return md.compute_center_of_mass(peptide_traj)[0]

if __name__ == "__main__":
    main()



