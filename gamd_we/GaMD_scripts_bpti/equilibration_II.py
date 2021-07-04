import re
import pickle as pk
from sys import stdout
from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *
from input_parameters import *
def create_vectors(x):
    x = str(x)
    x = x.replace("Vec3", "")
    x = re.findall('\d*\.?\d+',x)
    for i in range(0, len(x)): 
        x[i] = float(x[i]) 
    x = tuple(x)
    n = int(len(x)/3)
    x = [x[i * n:(i + 1) * n] for i in range((len(x) + n - 1) // n )]  
    return (x) 
parm = "system_TIP4P.prmtop"
npt_output_pdb = "system_npt_output.pdb" 
pdb_freq = pbd_freq_npt
npt_steps = npt_step
target_temp = 300
npt_pdb = "system_annealing_output_last_frame.pdb"
npt_init_pdb = PDBFile(npt_pdb)
prmtop = AmberPrmtopFile(parm)
npt_system = prmtop.createSystem(nonbondedMethod=PME,nonbondedCutoff=1*nanometer,constraints=HBonds)
barostat = MonteCarloBarostat(25.0*bar,target_temp*kelvin, 25)
npt_system.addForce(barostat)
npt_integrator = LangevinIntegrator(target_temp*kelvin, 1/picosecond, 2*femtoseconds)
npt_platform = Platform.getPlatformByName('CUDA')
npt_properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}
npt_simulation = Simulation(prmtop.topology, npt_system, npt_integrator, npt_platform, npt_properties)
npt_simulation.context.setPositions(npt_init_pdb.positions)
npt_simulation.context.setVelocitiesToTemperature(target_temp*kelvin)
with open('annealing_simulation_box_vectors.pkl', 'rb') as f:
    annealing_simulation_box_vectors = pk.load(f)
annealing_simulation_box_vectors = create_vectors(annealing_simulation_box_vectors)
npt_simulation.context.setPeriodicBoxVectors(annealing_simulation_box_vectors[0],annealing_simulation_box_vectors[1],annealing_simulation_box_vectors[2])
npt_last_frame = npt_output_pdb[:-4] + '_last_frame.pdb' 
npt_simulation.reporters.append(PDBReporter(npt_output_pdb,pdb_freq))
npt_simulation.reporters.append(PDBReporter(npt_last_frame,npt_steps))
npt_simulation.reporters.append(StateDataReporter(stdout,pdb_freq,step=True, time=True, potentialEnergy=True,totalSteps=npt_steps,temperature=True,progress=True,  remainingTime=True,speed=True, separator='\t'))
npt_simulation.minimizeEnergy()
npt_simulation.step(npt_steps)
npt_simulation.saveState('npt_simulation.state')
state = npt_simulation.context.getState()
print(state.getPeriodicBoxVectors())
npt_simulation_box_vectors = state.getPeriodicBoxVectors()
print(npt_simulation_box_vectors)
with open('npt_simulation_box_vectors.pkl', 'wb') as f:
    pk.dump(npt_simulation_box_vectors, f)       
print("Finished NPT Simulation")
