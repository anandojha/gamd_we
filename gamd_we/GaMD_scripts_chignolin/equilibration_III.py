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
parm = "system_TIP3P.prmtop"
nvt_output_pdb = "system_nvt_output.pdb"
pdb_freq = pbd_freq_nvt
nvt_steps = nvt_step
target_temp = 300
nvt_pdb = "system_npt_output_last_frame.pdb"       
nvt_init_pdb = PDBFile(nvt_pdb)
prmtop = AmberPrmtopFile(parm)
nvt_system = prmtop.createSystem(nonbondedMethod=PME,nonbondedCutoff=1*nanometer,constraints=HBonds)
nvt_integrator = LangevinIntegrator(target_temp*kelvin, 1/picosecond, 2*femtoseconds)
nvt_platform = Platform.getPlatformByName('CUDA')
nvt_properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}
nvt_simulation = Simulation(prmtop.topology, nvt_system, nvt_integrator,nvt_platform,nvt_properties)
nvt_simulation.context.setPositions(nvt_init_pdb.positions)
nvt_simulation.context.setVelocitiesToTemperature(target_temp*kelvin)
with open('npt_simulation_box_vectors.pkl', 'rb') as f:
    npt_simulation_box_vectors = pk.load(f)
npt_simulation_box_vectors = create_vectors(npt_simulation_box_vectors)
nvt_simulation.context.setPeriodicBoxVectors(npt_simulation_box_vectors[0], npt_simulation_box_vectors[1], npt_simulation_box_vectors[2])        
nvt_last_frame = nvt_output_pdb[:-4] + '_last_frame.pdb'
nvt_simulation.reporters.append(PDBReporter(nvt_output_pdb,pdb_freq))
nvt_simulation.reporters.append(PDBReporter(nvt_last_frame,nvt_steps))
nvt_simulation.reporters.append(StateDataReporter(stdout,pdb_freq,step=True, time=True,potentialEnergy=True, totalSteps=nvt_steps,temperature=True,progress=True,remainingTime=True,speed=True, separator='\t'))
nvt_simulation.minimizeEnergy()
nvt_simulation.step(nvt_steps)
nvt_simulation.saveState('nvt_simulation.state')
state = nvt_simulation.context.getState()
print(state.getPeriodicBoxVectors())
nvt_simulation_box_vectors = state.getPeriodicBoxVectors()
print(nvt_simulation_box_vectors)
with open('nvt_simulation_box_vectors.pkl', 'wb') as f:
    pk.dump(nvt_simulation_box_vectors, f)       
print("Finished NVT Simulation")
