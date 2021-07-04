import pickle as pk
from sys import stdout
from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *
from input_parameters import *
parm = "system_TIP3P.prmtop"
rst = "system_TIP3P.inpcrd"
annealing_output_pdb = "system_annealing_output.pdb"
annealing_steps = annealing_step
pdb_freq = pbd_freq_annealing
starting_temp = 0
target_temp = 300
temp_incr = 3
prmtop = AmberPrmtopFile(parm)
inpcrd = AmberInpcrdFile(rst)
annealing_system = prmtop.createSystem(nonbondedMethod=PME,nonbondedCutoff=1*nanometer,constraints=HBonds)
annealing_integrator = LangevinIntegrator(0*kelvin, 1/picosecond, 2*femtoseconds)
total_steps = ((target_temp / temp_incr) + 1) * annealing_steps
annealing_temp_range = int((target_temp / temp_incr) + 1)
annealing_platform = Platform.getPlatformByName('CUDA')
annealing_properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}
annealing_simulation = Simulation(prmtop.topology, annealing_system, annealing_integrator, annealing_platform,annealing_properties)
annealing_simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    annealing_simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
annealing_simulation.minimizeEnergy()
annealing_simulation.reporters.append(PDBReporter(annealing_output_pdb, pdb_freq))
simulated_annealing_last_frame = annealing_output_pdb[:-4] + '_last_frame.pdb'
annealing_simulation.reporters.append(PDBReporter(simulated_annealing_last_frame, total_steps))
annealing_simulation.reporters.append(StateDataReporter(stdout,pdb_freq,step=True, time=True,potentialEnergy=True, totalSteps=total_steps,temperature=True,progress=True,remainingTime=True,speed=True, separator='\t'))
temp = starting_temp
while temp <= target_temp:
    annealing_integrator.setTemperature(temp*kelvin)
    if temp == starting_temp:
        annealing_simulation.step(annealing_steps)
        annealing_simulation.saveState('annealing.state')
    else:
        annealing_simulation.loadState('annealing.state')
        annealing_simulation.step(annealing_steps)
    temp += temp_incr
state = annealing_simulation.context.getState()
print(state.getPeriodicBoxVectors())
annealing_simulation_box_vectors = state.getPeriodicBoxVectors()
print(annealing_simulation_box_vectors)
with open('annealing_simulation_box_vectors.pkl', 'wb') as f:
    pk.dump(annealing_simulation_box_vectors, f)       
print("Finshed NVT Simulated Annealing Simulation")
