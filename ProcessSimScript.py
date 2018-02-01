#!/usr/bin/python
###################
# This Script reads a text file with lines instructing the temperature simulations.
# Script is formatted with name interval_lentgh(ns) interval_steps and a list of temperatures 
# Expects folders for system to exist in the same directory as the execution location
# Each system folder contains a pdb with either name_box.pdb or name_min.pdb
# Creates if needed sub folders for each temperature within system folder
# Processing for each system: Initialize, followed by T equlib and trajectory simulation
#
#  K.E. Johnson, Pacific University
#  Dec 2017
#
##################
import sys
import traceback
import os.path
import time
# openMM Imports
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from openmmtools import integrators
# Garbage Collector!
import gc

def main():
#main will open and parse the simulation script file


# These two globalvariables are used to count simulation errors that should be trapped by 'try:'
# Used for each time simulation steps are run
	global errorcount
	errorcount = 0
	global maxerrors
	maxerrors = 5

# Check for presence of the script file
	if len(sys.argv)!=2:
		print("A single input argument is expected: The script file for the simulations.")
		sys.exit(1)
	infilepath = sys.argv[1]
	
	if not os.path.isfile(infilepath):
		print("File {} does not exist.  Exiting!".format(infilepath))
		sys.exit(1)

# read  and process one line at a time.
	with open(infilepath,'r') as scriptfile:
		
		for simcount, line in enumerate(scriptfile):
#read and parse each line of the script file
			list=line.split()	
			if len(list) < 4:
				print("Simulation Script line {} does not include minimum 4 items. Exiting!".format(simcount))
				sys.exit(1)
#Check for the input directory
			baseSimFile = list[0]
			if os.path.exists(baseSimFile):
				print("Base file name: {}".format(baseSimFile))
			else:
				print("Directory {} not found.  Exiting!".format(baseSimFile))
				sys.exit(1)
				
#Check for a valid simulation interval length
			try:
				intervalLength=float(list[1])
			except:
				print("{} is not a valid interval length.  Exiting!".format(list[1]))
				sys.exit(1)
			if intervalLength<0.01 or intervalLength>5.0:
				print("Simulation interval length must be between 0.1 and 5.0ns. Exiting!")
				sys.exit(1)
			print("Simulation interval: {}ns".format(intervalLength))
            
#Check for valid count of simulationintervals
			try:
				intervalCount=int(list[2])
			except:
				print("Interval Count must be an integer. Exiting!")
				sys.exit(1)
			if intervalCount<1 or intervalCount>20:
				print("Interval count be between 1 and 20.  Exiting!")
				sys.exit(1)

			print("Data saved at intervals:")
			for i in range(intervalCount):
				print("  {}ns".format((i+1)*intervalLength))
			print("")
#Read simulation temperatures into list.
			Temp_List=[]
			for i in range(3,len(list)):
				try:
					readTemp=float(list[i])
				except:
					print("{} is not a valid starting temperature.  Exiting!".format(list[1]))
					sys.exit(1)
				if readTemp<50 or readTemp>500:
					print("Starting Temperature must be between 50 and 500K. Exiting!")
					sys.exit(1)
				Temp_List.append(readTemp)
# List theTemperatures				
			print("Temperatures to simulate:")
			for t in Temp_List:
				print("  {}K".format(t))
			print("")

# Run simulations in the preexisting directories
			os.chdir(baseSimFile)

			Simulate(baseSimFile,Temp_List,intervalLength,intervalCount)
			
			os.chdir('../')
			
		print("Simulations complete.  Congratulations!")

#  END main


#########################################################	
def Simulate(FileBase,T_List,intervalTime,intervalCount):
	print("Starting simulation of {}".format(FileBase))
#Set up a new simulatiom
	
	mySimulation = InitSimulation(FileBase, T_List[0])
	
#Iterate over the Temperatures.  (Time iteration built in)
	
	for T_Sim in T_List:
		TempStr="{:.0f}".format(T_Sim)+"K"
		print("Starting simulations at T = "+TempStr)

# All tempeaature simulations in named directory.  Make if necessary
		if not os.path.exists(TempStr):
			os.mkdir(TempStr)
		os.chdir(TempStr)
		
		FileNameBase=FileBase+TempStr
		
# Check for .chk file at this temperature
		ChkFile=FileNameBase+".chk"
		if os.path.isfile(ChkFile):
			print('  Loading checkpoint...')
			with open(ChkFile, 'rb') as chk_File:
			    mySimulation.context.loadCheckpoint(chk_File.read())

# Now it's time to start the simu;ation!
		SetSimTemp(mySimulation,T_Sim,FileNameBase)
		
		result = SimAtT(mySimulation,T_Sim,intervalTime,intervalCount,FileNameBase)

#Move back out of temperature directory to ain system directory
		os.chdir('../')

#result is False when the error capture exceeds the maximum simulation faults at one temperature 
#Abort simulations of that molecule
		if not result:
			break


# Clean Up		
	del mySimulation.context
	del mySimulation.integrator
	del mySimulation
	gc.collect()
		
	print("Simulation of {} complete.\n\n".format(FileBase))

# END Simulate


#########################################################	
def InitSimulation(FileBase,TStart):
#Initialize and minimize simulation

	global errorcount
	errorcount=0
	global maxerrors

	print("  Initializing simulations.")	
	
# Read in pre-minimized pdb if present, otherwise the box.
	pdb_in=FileBase+"_box.pdb"
	pdb_min=FileBase+"_min.pdb"
	if os.path.isfile(pdb_min):
		print ("    Loading minimized file")
		pdb = PDBFile(pdb_min)
	elif os.path.isfile(pdb_in):
		print ("    Loading initial pdb file")
		pdb = PDBFile(pdb_in)
	else:
		print(pdb_in,": file not found\nScript Terminated")
		sys.exit(1)
		
# Set up for system, NoCutof for Amoeba, Constrain only x-H bond lengths
	force_Field_File = 'amoeba2013.xml'
	print ('    Loading Force Field')
	forcefield = ForceField(force_Field_File)
	
	print ('    Setting up system')
	system = forcefield.createSystem(	pdb.topology, 
										nonbondedMethod = PME,
										nonbondedCutoff=1.0*nanometer,
										polarization='extrapolated', 
										constraints = HBonds)

# Set up integrator
	print ('    Setting up integrator')
	Sim_Temp=TStart * kelvin
	integrator_Friction = 1.0/picoseconds	# friction is 1.0 per ps, Future possibliltuy more frequent collisions
	step_size = 1.0*femtoseconds        # Short time steps unconstrained water.  May need to go to 0.5ps 
	
#Experiment with new integrator.  Use 1.0fsstep??  7% Slower than Langevin but tollerates longer step?
	integrator = integrators.GeodesicBAOABIntegrator(2,Sim_Temp,
					integrator_Friction, 
					step_size)

# setup to run CUDA at mixed precision. Single precision for forces, integrate @ double precision.
	print ('Setting Platform to CUDA')
	platform = Platform.getPlatformByName('CUDA')
	properties = {'DeviceIndex': '0', 'CudaPrecision': 'mixed'}

# simulation is the object that accomplishes the MD work
	print ('Creating Simulation')
	mySimulation = Simulation(pdb.topology, 
							system, 
							integrator, 
							platform, 
							properties)

	mySimulation.context.setPositions(pdb.positions)
	
# Check if necessary and minimize.  Save result
	if not os.path.isfile(pdb_min):
		
		print('Minimizing System Energy' )
		done=False	
		while errorcount <= maxerrors and not done:
			try:
				mySimulation.minimizeEnergy()
				done=True
			except:
				print('Simulation error in minimization call')
				traceback.print_exc()
				errorcount += 1
				if errorcount > maxerrors:
					print('Maximum number of errors reached.  Terminating!')
					sys.exit(1)
				else:
					print('Resetting simulation positions at ', time.ctime())
					mySimulation.context.setPositions(pdb.positions)
					mySimulation.context.reinitialize()
					print('Restarting minimization')
			
		print ('Minimized Complete')
			
		print ('Saving PDB File')
		positions = mySimulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
		with open(pdb_min, 'w') as f:
			PDBFile.writeFile(mySimulation.topology, positions, f)
			
	print("  Initializatiom Complete\n")

	return mySimulation	
# Emd InitSimulation

	
#########################################################	
def SetSimTemp(mySimulation,Temp,F_name):
# Bring simulation to desired Temperature

	global errorcount
	errorcount=0
	global maxerrors

	print("Bringing Simulation to set Temperature: {} K".format(Temp))
	
	Sim_Temp = Temp*kelvin
	Equlib_Steps = 500
	
# Calculate Temperature Conversion. Used to calculate Temperature
	BOLTZ = MOLAR_GAS_CONSTANT_R
# calculate degrees of freedom in the system.  This used to calculate simulation temperature
	particles = mySimulation.system.getNumParticles()
	constraints = mySimulation.system.getNumConstraints()
	forces = [mySimulation.system.getForce(i) for i in range(mySimulation.system.getNumForces())]
	usesCMMR = any(isinstance(f, CMMotionRemover) for f in forces)
	DOF = 3*particles-constraints 
	if usesCMMR: DOF -= 3

# Set minimum error in sim temperature error reltive to set temperature
	Temp_tolerance = 0.0075

# Staring new simulation
	mySimulation.context.setTime(0)
	mySimulation.currentStep = 0
	mySimulation.integrator.setTemperature(Sim_Temp)

# Begin Cycle of testing simulation temperature to set point temperature
	Equlib_Cycle = 0
	Current_Temp=mySimulation.context.getState(getEnergy=True).getKineticEnergy()/(0.5*DOF*BOLTZ)
	Temp_Err = Current_Temp/Sim_Temp - 1

	while (abs(Temp_Err) > Temp_tolerance):
		Equlib_Cycle +=1
		print ("     Starting Thermal Equlib Step: {}".format(Equlib_Cycle))
		done=False
		while errorcount <= maxerrors and not done:
			try:
				mySimulation.step(Equlib_Steps)
				done=True
			except:
				print('Simulation error in Thermal Equilibration call')
				traceback.print_exc()
				errorcount += 1
				if errorcount > maxerrors:
					print('Maximum number of errors reached.  Terminating!')
					sys.exit(1)
				else:
					print('Resetting Simulation Temperature at ', time.ctime())
					statePositions = mySimulation.context.getState(getPositions=True, enforcePeriodicBox=True)
					mySimulation.context.reinitialize()
					mySimulation.context.setState(statePositions)
					mySimulation.context.setVelocitiesToTemperature(Sim_Temp,int(time.time()))
					print('Restarting Thermal Equilibration')
		
		Current_Temp=mySimulation.context.getState(getEnergy=True).getKineticEnergy()/(0.5*DOF*BOLTZ)
		Temp_Err = Current_Temp/Sim_Temp - 1
		print("         Simulation T=","{:.1f}".format(Current_Temp/kelvin),"K.  Temp Error = ", "{:+.2f}".format(Temp_Err*100),"%")
	print ('    Equilibrated Temp = {:.1f}K'.format(Current_Temp/kelvin))
	print ('	Thermal Equilibrated for ', mySimulation.currentStep, ' steps')

#Write .pdb & .chk file

	print ('Saving PDB File')
	positions = mySimulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
	with open(F_name + "_00.pdb","w") as f:
			PDBFile.writeFile(mySimulation.topology, positions, f)

	print('Saving Checkpoint File')
	with open(F_name+".chk",'wb') as f:
		f.write(mySimulation.context.createCheckpoint())

	return
# End SetSimTemp


#########################################################	
def SimAtT(mySimulation,Temp,intervalTime,intervalCount,FileNameBase):
#Run the requested number of simulations at set temperature.  Save incremental .pdb files.  Replace .chk file

	global errorcount
	errorcount = 0
	global maxerrors
	
	simTemp=Temp*kelvin

# Simulation steps reads integrator step length set in InitSim function
	Simulation_steps = int(intervalTime*nanoseconds/mySimulation.integrator.getStepSize())

# Reporting Step Definitions
	Trajectory_Step = int(Simulation_steps/200)  # 200 per 1 ns simulation
	Log_Step= int(Simulation_steps/100)          #100 per 1 ns simulatuon
# Initialize Reporters
	mySimulation.reporters.append(DCDReporter(FileNameBase + '_traj.dcd', Trajectory_Step))
	mySimulation.reporters.append(StateDataReporter(FileNameBase + '_sim.log', Log_Step, step=True, 
		time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, 
		temperature=True,  progress=True, 
		remainingTime=True, speed=True, totalSteps=Simulation_steps*intervalCount, separator='\t'))

# Reset simulation
	mySimulation.context.setTime(0)
	mySimulation.currentStep = 0
	
# Start the iterated simulation!
	print('Running Temperature Trajectory: ',FileNameBase)
	print('Startingtime: ',time.ctime())
			  
	for i in range(intervalCount):
		print("   Start Simulation Interval #: ", i+1)
		done=False
		while errorcount <= maxerrors and not done:
			try:
				mySimulation.step(Simulation_steps)
				done=True
			except:
				print('Simulation error in Trajectory Simuation call at ', time.ctime())
				traceback.print_exc()
				errorcount += 1
				if errorcount > maxerrors:
					print('Maximum number of errors reached.  Terminating!')
					break
#					sys.exit(1)
				else:
					print('  ReLoading checkpoint & resetting simulation...')
# Trial and error indicates the context should be reset.
# But first get the positions (translated into the box).
# Then reset, set positons, set velocities to new randomized temperature/ Boltzman distribution
					with open(FileNameBase+'.chk', 'rb') as chk_File:
					    mySimulation.context.loadCheckpoint(chk_File.read())
					statePositions = mySimulation.context.getState(getPositions=True, enforcePeriodicBox=True)
					mySimulation.context.reinitialize()
					mySimulation.context.setState(statePositions)
					mySimulation.context.setVelocitiesToTemperature(simTemp,int(time.time()))
					print('  RESTARTING Simulation')
		if errorcount > maxerrors:
			break
 
					
		print("      Saving Data")
		with open(FileNameBase+'.chk','wb') as f:
			f.write(mySimulation.context.createCheckpoint())
		positions = mySimulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
		with open(FileNameBase+"_" +"{:0>2d}".format(i+1)+".pdb","w") as f:
			PDBFile.writeFile(mySimulation.topology, positions, f)
		
#stop reporting progress and trajectory
	mySimulation.reporters.pop()
	mySimulation.reporters.pop()
	
	if errorcount<=maxerrors:
		print('Simulation Complete!')	
	else:
		print('Simulation '+FileNameBase+ 'Failed.')
	
	return (errorcount<=maxerrors)
# End SimAtT


			
if __name__ == '__main__':  
	main()						
									
				
