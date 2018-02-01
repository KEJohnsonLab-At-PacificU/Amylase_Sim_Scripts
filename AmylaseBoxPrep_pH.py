#prep_enzyme_box.py
# Kevin E. Johnson  October 2017, December 2017

#  Script to take a bare enzyme pdb structure (missing H atoms), add H atoms at specific pH, create water box with 10Angstrom buffer between molecules, and set ionic strength with NaCl.  Returns starting point for simulatuons as pdb file.

# The process uses two openmm tools:  First PDBFixer to add missing atoms. Second Modeller to add hydrogens, and create a padded solvent box. It then minimizes the bos using AMBER parameters.



import sys
import os.path
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from pdbfixer import *

	
# Get file name and pH from input line
	
if len(sys.argv)!=3:
	print("Two input arguments required: base_file_name & pH\nScript Terminated")
	sys.exit(1)
	
# check for presence of input file
BaseFileName=sys.argv[1]
pdb_ext=".pdb"
if os.path.isfile(BaseFileName+pdb_ext):
	print ("Base file name: ",BaseFileName)
else:
	print("Base_file_name not found\nScript Terminated")
	sys.exit(1)
# Check for valid pH
try:
	sim_pH=float(sys.argv[2])
except:
	print('pH not recognized\nScript Terminated')
	sys.exit(1)
	
#	First use pdbfixer to repair issueswith the pdb file
	
print('loading file into pdbfixer')
fixer=PDBFixer(filename=BaseFileName+pdb_ext)

print('Assume no missing residues.')
fixer.missingResidues = {}

print('Replace non-standard residues')	
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()

print('Find and replace missing heavyatoms')
fixer.findMissingAtoms()
fixer.addMissingAtoms()

#Determine the size of the molecule in 3 dimensions.  strip nm units.

print('Determining the molecule size')
molec_size=[]
for i in range(3):
	molec_size.append(max((pos[i] for pos in fixer.positions))-min((pos[i] for pos in fixer.positions)))
# 1.0 nm pad in each direction.  UNITLESS!?!
pad = 1.0*nanometer
print('Solvent box padded {} in each dimension'.format(pad))
simsize = Vec3(molec_size[0]+pad,molec_size[1]+pad,molec_size[2]+pad)/nanometer


print('Modeling with Amber')
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
modeller = Modeller(fixer.topology, fixer.positions)

print('Adding hydrogens...')
modeller.addHydrogens(forcefield,pH=sim_pH)
print('Adding solvent...')
modeller.addSolvent(forcefield, boxSize=simsize, ionicStrength=0.050*molar, positiveIon='Na+', negativeIon='Cl-')
print('Minimizing with AMBER...')
platform = Platform.getPlatformByName('CUDA')
properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)
integrator = VerletIntegrator(0.001*picoseconds)
simulation = Simulation(modeller.topology, system, integrator, platform, properties)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy(maxIterations=500)

print('Saving...')
positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
with open(BaseFileName+"_box"+pdb_ext, 'w') as f:
	PDBFile.writeFile(simulation.topology, positions, f)


print('Done')
