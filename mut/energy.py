from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import sys
import numpy as np
import time
import mdtraj

# USAGE: python energy.py <pdb file>

def main():

    filename = sys.argv[1]

    #print "Opening:", filename

    pdb = PDBFile(filename)
    top = pdb.getTopology()
    positions = np.array(pdb.positions) #pdb.getPositions(asNumpy=True)
    numAtoms = len(positions)

    #print "Number of atoms:", numAtoms
    #print "Number of residues:", top.getNumResidues()

    positions = np.reshape(positions, (3*numAtoms,1))

    # run file through pdb fixer first

    #forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
    forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')
    #forcefield = app.ForceField('amber03.xml', 'amber03_obc.xml')
    #forcefield = app.ForceField('amber10.xml', 'amber10_obc.xml')
    modeller = Modeller(pdb.topology, pdb.positions)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic, constraints=None)


    forces = system.getForces()
    i = 0
    for f in forces:
        f.setForceGroup(i)
        i = i + 1


    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    #integrator = VerletIntegrator(0.002*picoseconds)
    platform = Platform.getPlatformByName('CPU')
    simulation = Simulation(modeller.topology, system, integrator, platform)

    simulation.context.setPositions(modeller.positions)
    #print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

    #printForces(simulation)
    print(simulation.context.getState(getEnergy=True, groups=2**3+2**4).getPotentialEnergy())

def printForces(simulation):

    
    for i in range(simulation.system.getNumForces()):
        f_name = simulation.system.getForce(i).__class__.__name__
        s0 = simulation.context.getState(getEnergy=True, groups=2**i)
        print(f_name + ": " + str(s0.getPotentialEnergy()))
    

    print("total:", simulation.context.getState(getEnergy=True).getPotentialEnergy())


if __name__ == "__main__":
    main()



