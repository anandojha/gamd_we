import matplotlib.image as mpimg
import matplotlib.style as style
import matplotlib.pyplot as plt
from matplotlib import rcParams
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import seaborn as sns
from math import exp
import pandas as pd
import mdtraj as md
import pickle as pk
import numpy as np
import statistics
import itertools
import fileinput
import fnmatch
import shutil
import random
import math
import os
import re


def create_vectors(x):

    """
    Extracts peridic box information from the
    given line.

    """
    x = str(x)
    x = x.replace("Vec3", "")
    x = re.findall("\d*\.?\d+", x)
    for i in range(0, len(x)):
        x[i] = float(x[i])
    x = tuple(x)
    n = int(len(x) / 3)
    x = [x[i * n : (i + 1) * n] for i in range((len(x) + n - 1) // n)]
    return x


def prepare_bpti():

    """

    Prepares the Bovine pancreatic trypsin inhibitor
    system for Gaussian Accelerated Molecular Dynamics
    (GaMD) simulations. Downloads the pdb structure from
    http://ambermd.org/tutorials/advanced/tutorial22/files/5PTI-DtoH-dry.pdb
    and parameterizes it using General Amber Force Field
    (GAFF).

    """

    os.system(
        "curl -O http://ambermd.org/tutorials/advanced/tutorial22/files/5PTI-DtoH-dry.pdb"
    )
    os.system("rm -rf system_inputs")
    # Removes any existing directory named system_inputs
    os.system("mkdir system_inputs")
    # Creates a directory named system_inputs
    cwd = os.getcwd()
    target_dir = cwd + "/" + "system_inputs"
    # save the tleap script to file
    with open("input_TIP4P.leap", "w") as f:
        f.write(
            """
    source leaprc.protein.ff14SB
    loadOff solvents.lib
    loadOff tip4pbox.off
    loadOff tip4pewbox.off
    source leaprc.water.tip4pew
    HOH = TP4
    pdb = loadpdb 5PTI-DtoH-dry.pdb
    bond pdb.55.SG pdb.5.SG
    bond pdb.30.SG pdb.51.SG
    bond pdb.14.SG pdb.38.SG
    charge pdb
    addions2 pdb Cl- 6
    charge pdb
    solvatebox pdb TIP4PEWBOX 12.0
    saveamberparm pdb system_TIP4P.prmtop system_TIP4P.inpcrd
    saveamberparm pdb system_TIP4P.parm7 system_TIP4P.rst7
    savepdb pdb system_TIP4P.pdb
    quit
    """
        )
    os.system("tleap -f input_TIP4P.leap")
    os.system("rm -rf leap.log")
    shutil.copy(
        cwd + "/" + "system_TIP4P.inpcrd",
        target_dir + "/" + "system_TIP4P.inpcrd",
    )
    shutil.copy(
        cwd + "/" + "system_TIP4P.parm7",
        target_dir + "/" + "system_TIP4P.parm7",
    )
    shutil.copy(
        cwd + "/" + "system_TIP4P.pdb", target_dir + "/" + "system_TIP4P.pdb"
    )
    shutil.copy(
        cwd + "/" + "system_TIP4P.prmtop",
        target_dir + "/" + "system_TIP4P.prmtop",
    )
    shutil.copy(
        cwd + "/" + "system_TIP4P.rst7", target_dir + "/" + "system_TIP4P.rst7"
    )
    shutil.copy(
        cwd + "/" + "input_TIP4P.leap", target_dir + "/" + "input_TIP4P.leap"
    )
    shutil.copy(
        cwd + "/" + "5PTI-DtoH-dry.pdb", target_dir + "/" + "5PTI-DtoH-dry.pdb"
    )
    os.system("rm -rf system_TIP4P.inpcrd")
    os.system("rm -rf system_TIP4P.parm7")
    os.system("rm -rf system_TIP4P.pdb")
    os.system("rm -rf system_TIP4P.rst7")
    os.system("rm -rf system_TIP4P.prmtop")
    os.system("rm -rf input_TIP4P.leap")
    os.system("rm -rf 5PTI-DtoH-dry.pdb")


def simulated_annealing(
    parm="system_TIP4P.prmtop",
    rst="system_TIP4P.inpcrd",
    annealing_output_pdb="system_annealing_output.pdb",
    annealing_steps=100000,
    pdb_freq=100000,
    starting_temp=0,
    target_temp=300,
    temp_incr=3,
):

    """

    Performs simulated annealing of the system from
    0K to 300 K (default) using OpenMM MD engine and
    saves the last frame of the simulation to be
    accessed by the next simulation.

    Parameters
    ----------
    parm: str
        System's topology file

    rst: str
        System's coordinate file

    annealing_output_pdb: str
        System's output trajectory file

    annealing_steps: int
        Aneealing steps at each temperatrure jump

    pdb_freq: int
        Trajectory to be saved after every pdb_freq steps

    starting_temp: int
        Initial temperature of Simulated Annealing

    target_temp: int
        Final temperature of Simulated Annealing

    temp_incr: int
        Temmperature increase for every step

    """

    prmtop = AmberPrmtopFile(parm)
    inpcrd = AmberInpcrdFile(rst)
    annealing_system = prmtop.createSystem(
        nonbondedMethod=PME, nonbondedCutoff=1 * nanometer, constraints=HBonds
    )
    annealing_integrator = LangevinIntegrator(
        0 * kelvin, 1 / picosecond, 2 * femtoseconds
    )
    total_steps = ((target_temp / temp_incr) + 1) * annealing_steps
    annealing_temp_range = int((target_temp / temp_incr) + 1)
    annealing_platform = Platform.getPlatformByName("CUDA")
    annealing_properties = {"CudaDeviceIndex": "0", "CudaPrecision": "mixed"}
    annealing_simulation = Simulation(
        prmtop.topology,
        annealing_system,
        annealing_integrator,
        annealing_platform,
        annealing_properties,
    )
    annealing_simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        annealing_simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    annealing_simulation.minimizeEnergy()
    annealing_simulation.reporters.append(
        PDBReporter(annealing_output_pdb, pdb_freq)
    )
    simulated_annealing_last_frame = (
        annealing_output_pdb[:-4] + "_last_frame.pdb"
    )
    annealing_simulation.reporters.append(
        PDBReporter(simulated_annealing_last_frame, total_steps)
    )
    annealing_simulation.reporters.append(
        StateDataReporter(
            stdout,
            pdb_freq,
            step=True,
            time=True,
            potentialEnergy=True,
            totalSteps=total_steps,
            temperature=True,
            progress=True,
            remainingTime=True,
            speed=True,
            separator="\t",
        )
    )
    temp = starting_temp
    while temp <= target_temp:
        annealing_integrator.setTemperature(temp * kelvin)
        if temp == starting_temp:
            annealing_simulation.step(annealing_steps)
            annealing_simulation.saveState("annealing.state")
        else:
            annealing_simulation.loadState("annealing.state")
            annealing_simulation.step(annealing_steps)
        temp += temp_incr
    state = annealing_simulation.context.getState()
    print(state.getPeriodicBoxVectors())
    annealing_simulation_box_vectors = state.getPeriodicBoxVectors()
    print(annealing_simulation_box_vectors)
    with open("annealing_simulation_box_vectors.pkl", "wb") as f:
        pk.dump(annealing_simulation_box_vectors, f)
    print("Finshed NVT Simulated Annealing Simulation")


def npt_equilibration(
    parm="system_TIP4P.prmtop",
    npt_output_pdb="system_npt_output.pdb",
    pdb_freq=500000,
    npt_steps=5000000,
    target_temp=300,
    npt_pdb="system_annealing_output_last_frame.pdb",
):

    """

    Performs NPT equilibration MD of the system
    using OpenMM MD engine and saves the last
    frame of the simulation to be accessed by
    the next simulation.

    Parameters
    ----------
    parm: str
        System's topology file

    npt_output_pdb: str
        System's output trajectory file

    pdb_freq: int
        Trajectory to be saved after every pdb_freq steps

    npt_steps: int
        NPT simulation steps

    target_temp: int
        Temperature for MD simulation

    npt_pdb: str
        Last frame of the simulation

    """

    npt_init_pdb = PDBFile(npt_pdb)
    prmtop = AmberPrmtopFile(parm)
    npt_system = prmtop.createSystem(
        nonbondedMethod=PME, nonbondedCutoff=1 * nanometer, constraints=HBonds
    )
    barostat = MonteCarloBarostat(25.0 * bar, target_temp * kelvin, 25)
    npt_system.addForce(barostat)
    npt_integrator = LangevinIntegrator(
        target_temp * kelvin, 1 / picosecond, 2 * femtoseconds
    )
    npt_platform = Platform.getPlatformByName("CUDA")
    npt_properties = {"CudaDeviceIndex": "0", "CudaPrecision": "mixed"}
    npt_simulation = Simulation(
        prmtop.topology,
        npt_system,
        npt_integrator,
        npt_platform,
        npt_properties,
    )
    npt_simulation.context.setPositions(npt_init_pdb.positions)
    npt_simulation.context.setVelocitiesToTemperature(target_temp * kelvin)
    with open("annealing_simulation_box_vectors.pkl", "rb") as f:
        annealing_simulation_box_vectors = pk.load(f)
    annealing_simulation_box_vectors = create_vectors(
        annealing_simulation_box_vectors
    )
    npt_simulation.context.setPeriodicBoxVectors(
        annealing_simulation_box_vectors[0],
        annealing_simulation_box_vectors[1],
        annealing_simulation_box_vectors[2],
    )
    npt_last_frame = npt_output_pdb[:-4] + "_last_frame.pdb"
    npt_simulation.reporters.append(PDBReporter(npt_output_pdb, pdb_freq))
    npt_simulation.reporters.append(PDBReporter(npt_last_frame, npt_steps))
    npt_simulation.reporters.append(
        StateDataReporter(
            stdout,
            pdb_freq,
            step=True,
            time=True,
            potentialEnergy=True,
            totalSteps=npt_steps,
            temperature=True,
            progress=True,
            remainingTime=True,
            speed=True,
            separator="\t",
        )
    )
    npt_simulation.minimizeEnergy()
    npt_simulation.step(npt_steps)
    npt_simulation.saveState("npt_simulation.state")
    state = npt_simulation.context.getState()
    print(state.getPeriodicBoxVectors())
    npt_simulation_box_vectors = state.getPeriodicBoxVectors()
    print(npt_simulation_box_vectors)
    with open("npt_simulation_box_vectors.pkl", "wb") as f:
        pk.dump(npt_simulation_box_vectors, f)
    print("Finished NPT Simulation")


def nvt_equilibration(
    parm="system_TIP4P.prmtop",
    nvt_output_pdb="system_nvt_output.pdb",
    pdb_freq=500000,
    nvt_steps=5000000,
    target_temp=300,
    nvt_pdb="system_npt_output_last_frame.pdb",
):

    """

    Performs NVT equilibration MD of the system
    using OpenMM MD engine  saves the last
    frame of the simulation to be accessed by
    the next simulation.

    Parameters
    ----------
    parm: str
        System's topology file

    nvt_output_pdb: str
        System's output trajectory file

    pdb_freq: int
        Trajectory to be saved after every pdb_freq steps

    nvt_steps: int
        NVT simulation steps

    target_temp: int
        Temperature for MD simulation

    nvt_pdb: str
        Last frame of the simulation

    """

    nvt_init_pdb = PDBFile(nvt_pdb)
    prmtop = AmberPrmtopFile(parm)
    nvt_system = prmtop.createSystem(
        nonbondedMethod=PME, nonbondedCutoff=1 * nanometer, constraints=HBonds
    )
    nvt_integrator = LangevinIntegrator(
        target_temp * kelvin, 1 / picosecond, 2 * femtoseconds
    )
    nvt_platform = Platform.getPlatformByName("CUDA")
    nvt_properties = {"CudaDeviceIndex": "0", "CudaPrecision": "mixed"}
    nvt_simulation = Simulation(
        prmtop.topology,
        nvt_system,
        nvt_integrator,
        nvt_platform,
        nvt_properties,
    )
    nvt_simulation.context.setPositions(nvt_init_pdb.positions)
    nvt_simulation.context.setVelocitiesToTemperature(target_temp * kelvin)
    with open("npt_simulation_box_vectors.pkl", "rb") as f:
        npt_simulation_box_vectors = pk.load(f)
    npt_simulation_box_vectors = create_vectors(npt_simulation_box_vectors)
    nvt_simulation.context.setPeriodicBoxVectors(
        npt_simulation_box_vectors[0],
        npt_simulation_box_vectors[1],
        npt_simulation_box_vectors[2],
    )
    nvt_last_frame = nvt_output_pdb[:-4] + "_last_frame.pdb"
    nvt_simulation.reporters.append(PDBReporter(nvt_output_pdb, pdb_freq))
    nvt_simulation.reporters.append(PDBReporter(nvt_last_frame, nvt_steps))
    nvt_simulation.reporters.append(
        StateDataReporter(
            stdout,
            pdb_freq,
            step=True,
            time=True,
            potentialEnergy=True,
            totalSteps=nvt_steps,
            temperature=True,
            progress=True,
            remainingTime=True,
            speed=True,
            separator="\t",
        )
    )
    nvt_simulation.minimizeEnergy()
    nvt_simulation.step(nvt_steps)
    nvt_simulation.saveState("nvt_simulation.state")
    state = nvt_simulation.context.getState()
    print(state.getPeriodicBoxVectors())
    nvt_simulation_box_vectors = state.getPeriodicBoxVectors()
    print(nvt_simulation_box_vectors)
    with open("nvt_simulation_box_vectors.pkl", "wb") as f:
        pk.dump(nvt_simulation_box_vectors, f)
    print("Finished NVT Simulation")


def run_equilibration():

    """

    Runs systematic simulated annealing followed by
    NPT and NVT equilibration MD simulation.

    """

    cwd = os.getcwd()
    target_dir = cwd + "/" + "equilibration"
    os.system("rm -rf equilibration")
    os.system("mkdir equilibration")
    shutil.copy(
        cwd + "/" + "system_inputs" + "/" + "system_TIP4P.inpcrd",
        target_dir + "/" + "system_TIP4P.inpcrd",
    )
    shutil.copy(
        cwd + "/" + "system_inputs" + "/" + "system_TIP4P.parm7",
        target_dir + "/" + "system_TIP4P.parm7",
    )
    shutil.copy(
        cwd + "/" + "system_inputs" + "/" + "system_TIP4P.pdb",
        target_dir + "/" + "system_TIP4P.pdb",
    )
    shutil.copy(
        cwd + "/" + "system_inputs" + "/" + "system_TIP4P.prmtop",
        target_dir + "/" + "system_TIP4P.prmtop",
    )
    shutil.copy(
        cwd + "/" + "system_inputs" + "/" + "system_TIP4P.rst7",
        target_dir + "/" + "system_TIP4P.rst7",
    )
    shutil.copy(
        cwd + "/" + "system_inputs" + "/" + "5PTI-DtoH-dry.pdb",
        target_dir + "/" + "5PTI-DtoH-dry.pdb",
    )
    shutil.copy(
        cwd + "/" + "system_inputs" + "/" + "input_TIP4P.leap",
        target_dir + "/" + "input_TIP4P.leap",
    )
    os.chdir(target_dir)
    simulated_annealing()
    npt_equilibration()
    nvt_equilibration()
    os.system("rm -rf system_TIP4P.inpcrd")
    os.system("rm -rf system_TIP4P.parm7")
    os.system("rm -rf system_TIP4P.pdb")
    os.system("rm -rf system_TIP4P.rst7")
    os.system("rm -rf system_TIP4P.prmtop")
    os.system("rm -rf 5PTI-DtoH-dry.pdb")
    os.system("rm -rf input_TIP4P.leap")
    os.system("rm -rf alanine-dipeptide-nowater.pdb")
    os.system("rm -rf input_TIP3P.leap")
    os.chdir(cwd)


def create_starting_structures():

    """
    Prepares starting structures for Amber GaMD simulations.
    All input files required to run Amber GaMD simulations are
    placed in the starting_structures directory.

    """

    cwd = os.getcwd()
    target_dir = cwd + "/" + "starting_structures"
    os.system("rm -rf starting_structures")
    os.system("mkdir starting_structures")
    shutil.copy(
        cwd + "/" + "equilibration" + "/" + "system_nvt_output_last_frame.pdb",
        target_dir + "/" + "system_nvt_output_last_frame.pdb",
    )
    os.chdir(target_dir)
    os.system(
        "pdb4amber -i system_nvt_output_last_frame.pdb -o intermediate_temp.pdb"
    )
    os.system("rm -rf intermediate_temp_renum.txt")
    os.system("rm -rf intermediate_temp_sslink")
    os.system("rm -rf intermediate_temp_nonprot.pdb")
    remove_words = ["H   ARG A   1"]
    with open("intermediate_temp.pdb") as oldfile, open(
        "intermediate.pdb", "w"
    ) as newfile:
        for line in oldfile:
            if not any(word in line for word in remove_words):
                newfile.write(line)
    # Save the tleap script to file
    with open("final_input_TIP4P.leap", "w") as f:
        f.write(
            """
    source leaprc.protein.ff14SB
    source leaprc.water.tip4pew
    pdb = loadpdb intermediate.pdb
    charge pdb
    saveamberparm pdb system_final.prmtop system_final.inpcrd
    saveamberparm pdb system_final.parm7 system_final.rst7
    savepdb pdb system_final.pdb
    quit
    """
        )
    os.system("tleap -f final_input_TIP4P.leap")
    os.system("rm -rf leap.log")
    os.system("rm -rf leap.log")
    os.system("rm -rf intermediate.pdb")
    os.system("rm -rf intermediate_temp.pdb")
    os.system("rm -rf system_nvt_output_last_frame.pdb")
    os.chdir(cwd)


def add_vec_inpcrd():

    """

    Adds box dimensions captured from the last saved
    frame of the NVT simulations to the inpcrd file.
    Only to be used when the box dimensions are not
    present in the inpcrd file.

    """

    cwd = os.getcwd()
    target_dir = cwd + "/" + "starting_structures"
    shutil.copy(
        cwd + "/" + "equilibration" + "/" + "nvt_simulation_box_vectors.pkl",
        target_dir + "/" + "nvt_simulation_box_vectors.pkl",
    )

    os.chdir(target_dir)
    with open("nvt_simulation_box_vectors.pkl", "rb") as f:
        nvt_simulation_box_vectors = pk.load(f)
    nvt_simulation_box_vectors = create_vectors(nvt_simulation_box_vectors)
    vectors = (
        (nvt_simulation_box_vectors[0][0]) * 10,
        (nvt_simulation_box_vectors[1][1]) * 10,
        (nvt_simulation_box_vectors[2][2]) * 10,
    )
    vectors = (
        round(vectors[0], 7),
        round(vectors[1], 7),
        round(vectors[2], 7),
    )
    last_line = (
        "  "
        + str(vectors[0])
        + "  "
        + str(vectors[1])
        + "  "
        + str(vectors[2])
        + "  90.0000000"
        + "  90.0000000"
        + "  90.0000000"
    )
    with open("system_final.inpcrd", "a+") as f:
        f.write(last_line)
    os.system("rm -rf nvt_simulation_box_vectors.pkl")
    os.chdir(cwd)


def add_vec_prmtop():

    """

    Adds box dimensions captured from the last saved
    frame of the NVT simulations to the prmtop file.
    Only to be used when the box dimensions are not
    present in the prmtop file.

    """

    cwd = os.getcwd()
    target_dir = cwd + "/" + "starting_structures"
    shutil.copy(
        cwd + "/" + "equilibration" + "/" + "nvt_simulation_box_vectors.pkl",
        target_dir + "/" + "nvt_simulation_box_vectors.pkl",
    )

    os.chdir(target_dir)
    with open("nvt_simulation_box_vectors.pkl", "rb") as f:
        nvt_simulation_box_vectors = pk.load(f)
    nvt_simulation_box_vectors = create_vectors(nvt_simulation_box_vectors)
    vectors = (
        nvt_simulation_box_vectors[0][0],
        nvt_simulation_box_vectors[1][1],
        nvt_simulation_box_vectors[2][2],
    )
    vectors = round(vectors[0], 7), round(vectors[1], 7), round(vectors[2], 7)
    oldbeta = "9.00000000E+01"
    x = str(vectors[0]) + str(0) + "E+" + "01"
    y = str(vectors[1]) + str(0) + "E+" + "01"
    z = str(vectors[2]) + str(0) + "E+" + "01"
    line1 = "%FLAG BOX_DIMENSIONS"
    line2 = "%FORMAT(5E16.8)"
    line3 = "  " + oldbeta + "  " + x + "  " + y + "  " + z
    with open("system_final.prmtop") as i, open(
        "system_intermediate_final.prmtop", "w"
    ) as f:
        for line in i:
            if line.startswith("%FLAG RADIUS_SET"):
                line = line1 + "\n" + line2 + "\n" + line3 + "\n" + line
            f.write(line)
    os.system("rm -rf system_final.prmtop")
    os.system("mv system_intermediate_final.prmtop system_final.prmtop")
    os.system("rm -rf nvt_simulation_box_vectors.pkl")
    os.chdir(cwd)


def create_filetree(
    nst_lim=251000000,
    ntw_x=1000,
    nt_cmd=1000000,
    n_teb=1000000,
    n_tave=50000,
    ntcmd_prep=200000,
    nteb_prep=200000,
):

    """

    Creates a directory named gamd_simulations. Inside
    this directory, there are subdirectories for dihedral,
    dual and total potential-boosted GaMD with upper and
    lower threshold boosts separately.

    Parameters
    ----------

    nst_lim: int
        Total simulation time including preparatory simulation.
        For example, if nst_lim = 251000000, then, we may have
        2 ns of preparatory simulation i.e. 1000000 preparation steps
        and 500 ns of GaMD simulation i.e. 250000000 simulation steps

    ntw_x: int
        Saving coordinates of the simulation every ntw_x
        timesteps. For example, 2 ps implies 1000 timesteps

    nt_cmd: int
        Number of initial MD simulation step, 2 ns of
        preparatory simulation requires 1000000 preparation
        timesteps

    n_teb: int
        Number of biasing MD simulation steps

    n_tave: int
        Number of simulation steps used to calculate the
        average and standard deviation of potential energies

    ntcmd_prep: int
        Number of preparation conventional molecular dynamics
        steps.This is used for system equilibration and
        potential energies are not collected for statistics

    nteb_prep: int
        Number of preparation biasing molecular dynamics
        simulation steps. This is used for system
        equilibration

    """

    cwd = os.getcwd()
    os.system("rm -rf gamd_simulations")
    os.system("mkdir gamd_simulations")
    os.chdir(cwd + "/" + "gamd_simulations")
    source_dir = cwd + "/" + "starting_structures"
    target_dir = cwd + "/" + "gamd_simulations"
    dir_list = [
        "dihedral_threshold_lower",
        "dihedral_threshold_upper",
        "dual_threshold_lower",
        "dual_threshold_upper",
        "total_threshold_lower",
        "total_threshold_upper",
    ]
    for i in range(len(dir_list)):
        os.mkdir(dir_list[i])
        os.chdir(target_dir + "/" + dir_list[i])
        shutil.copy(
            source_dir + "/" + "system_final.inpcrd",
            target_dir + "/" + dir_list[i] + "/" + "system_final.inpcrd",
        )
        shutil.copy(
            source_dir + "/" + "system_final.prmtop",
            target_dir + "/" + dir_list[i] + "/" + "system_final.prmtop",
        )
        if "lower" in dir_list[i]:
            i_E = 1
        if "upper" in dir_list[i]:
            i_E = 2
        if "total" in dir_list[i]:
            i_gamd = 1
        if "dihedral" in dir_list[i]:
            i_gamd = 2
        if "dual" in dir_list[i]:
            i_gamd = 3
        with open("md.in", "w") as f:
            f.write("&cntrl" + "\n")
            f.write("  imin = 0, irest = 0, ntx = 1," + "\n")
            f.write("  nstlim = " + str(nst_lim) + ", dt = 0.002," + "\n")
            f.write("  ntc = 2, ntf = 2, tol = 0.000001," + "\n")
            f.write("  iwrap = 1, ntb = 1, cut = 8.0," + "\n")
            f.write("  ntt = 3, temp0 = 300.0, gamma_ln = 1.0, " + "\n")
            f.write(
                "  ntpr = 500, ntwx = " + str(ntw_x) + ", ntwr = 500," + "\n"
            )
            f.write("  ntxo = 2, ioutfm = 1, ig = -1, ntwprt = 0," + "\n")
            f.write(
                "  igamd = "
                + str(i_gamd)
                + ", iE = "
                + str(i_E)
                + ", irest_gamd = 0,"
                + "\n"
            )
            f.write(
                "  ntcmd = "
                + str(nt_cmd)
                + ", nteb = "
                + str(n_teb)
                + ", ntave = "
                + str(n_tave)
                + ","
                + "\n"
            )
            f.write(
                "  ntcmdprep = "
                + str(ntcmd_prep)
                + ", ntebprep = "
                + str(nteb_prep)
                + ","
                + "\n"
            )
            f.write("  sigma0D = 6.0, sigma0P = 6.0" + " \n")
            f.write("&end" + "\n")
            os.chdir(target_dir)
    os.chdir(cwd)


def run_simulations():

    """

    Runs GaMD simulations for each of the dihedral, dual and total
    potential boosts for both thresholds i.e. upper and lower potential
    thresholds. (Remember to check md.in files for further details and
    flag information).

    """

    cwd = os.getcwd()
    os.chdir(cwd + "/" + "gamd_simulations")
    os.chdir(cwd + "/" + "gamd_simulations" + "/" + "dihedral_threshold_lower")
    os.system(
        "pmemd.cuda -O -i md.in -o system_final.out -p system_final.prmtop -c system_final.inpcrd -r system_final.rst -x system_final.nc"
    )
    os.chdir(cwd + "/" + "gamd_simulations" + "/" + "dihedral_threshold_upper")
    os.system(
        "pmemd.cuda -O -i md.in -o system_final.out -p system_final.prmtop -c system_final.inpcrd -r system_final.rst -x system_final.nc"
    )
    os.chdir(cwd + "/" + "gamd_simulations" + "/" + "dual_threshold_lower")
    os.system(
        "pmemd.cuda -O -i md.in -o system_final.out -p system_final.prmtop -c system_final.inpcrd -r system_final.rst -x system_final.nc"
    )
    os.chdir(cwd + "/" + "gamd_simulations" + "/" + "dual_threshold_upper")
    os.system(
        "pmemd.cuda -O -i md.in -o system_final.out -p system_final.prmtop -c system_final.inpcrd -r system_final.rst -x system_final.nc"
    )
    os.chdir(cwd + "/" + "gamd_simulations" + "/" + "total_threshold_lower")
    os.system(
        "pmemd.cuda -O -i md.in -o system_final.out -p system_final.prmtop -c system_final.inpcrd -r system_final.rst -x system_final.nc"
    )
    os.chdir(cwd + "/" + "gamd_simulations" + "/" + "total_threshold_upper")
    os.system(
        "pmemd.cuda -O -i md.in -o system_final.out -p system_final.prmtop -c system_final.inpcrd -r system_final.rst -x system_final.nc"
    )
    os.chdir(cwd + "/" + "gamd_simulations")
    os.chdir(cwd)


"""
prepare_bpti()
run_equilibration()
create_starting_structures()
add_vec_inpcrd()
add_vec_prmtop()
create_filetree()
run_simulations()
"""
