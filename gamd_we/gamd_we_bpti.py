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


def create_data_files(
    index_phi=[200, 202, 204, 207],
    index_psi=[584, 586, 588, 591],
    jump=10,
    traj="system_final.nc",
    topology="system_final.prmtop",
    T=300,
):

    """

    Extracts data from GaMD log files and saves them as
    weights.dat, Psi.dat and Phi_Psi.dat. gamd.log file
    contains data excluding the initial equilibration MD
    simulation steps but trajectory output file has all
    the trajectories including the initial equilibration
    MD steps.  This part has ben taken care to make the
    data consistent.

    Parameters
    ----------

    index_phi : list
        List of indices of the phi dihedral to
        be calculated

    index_phi : list
        List of indices of the phi dihedral to
        be calculated

    jump: int
        Every nth frame to be considered for
        reweighting

    traj: str
        System's trajectory file

    topology: str
        System's topology file

    T: int
        MD simulation temperature

    """

    # To make data consistent with gamd.log and .nc file
    factor = 0.001987 * T
    with open("md.in") as f:
        lines = f.readlines()
    for i in lines:
        if "nstlim =" in i:
            nstlim_line = i
        if "ntcmd =" in i:
            ntcmd_line = i
        if "ntwx =" in i:
            ntwx_line = i
    x = re.findall(r"\b\d+\b", ntcmd_line)
    ntcmd = int(x[0])
    x = re.findall(r"\b\d+\b", nstlim_line)
    nstlim = int(x[0])
    x = re.findall(r"\b\d+\b", ntwx_line)
    ntwx = int(x[1])
    # From the .nc trajectory files, we will not consider ntcmd trajectories
    leave_frames = int(ntcmd / ntwx)
    no_frames = int(nstlim / ntwx)
    frame_indices = []
    for i in range(leave_frames, no_frames, jump):
        frame_indices.append(i)
    # Recheck conditions
    file = open("gamd.log", "r")
    number_of_lines = 0
    for line in file:
        line = line.strip("\n")
        number_of_lines += 1
    file.close()
    f = open("gamd.log")
    fourth_line = f.readlines()[3]
    if str(ntcmd) in fourth_line:
        datapoints = number_of_lines - 4
    if not str(ntcmd) in fourth_line:
        datapoints = number_of_lines - 3
    print(datapoints == int((nstlim - ntcmd) / ntwx))
    # Creating Psi.dat and Phi_Psi.dat
    traj = md.load(traj, top=topology)
    traj = traj[leave_frames:no_frames:jump]
    phi = md.compute_dihedrals(
        traj,
        indices=np.array(
            [[index_phi[0], index_phi[1], index_phi[2], index_phi[3]]]
        ),
    )
    phi = np.array([math.degrees(i) for i in phi])  # radians to degrees
    psi = md.compute_dihedrals(
        traj,
        indices=np.array(
            [[index_psi[0], index_psi[1], index_psi[2], index_psi[3]]]
        ),
    )
    psi = np.array([math.degrees(i) for i in psi])  # radians to degrees
    phi = np.array([math.degrees(i) for i in phi])  # radians to degrees
    psi = np.array([math.degrees(i) for i in psi])  # radians to degrees
    df_psi = pd.DataFrame(phi, columns=["Psi"])
    df_psi = df_psi.tail(int(datapoints))
    df_psi.to_csv("Psi.dat", sep="\t", index=False, header=False)
    df_phi = pd.DataFrame(psi, columns=["Phi"])
    df_phi = df_phi.tail(int(datapoints))
    df_phi_psi = pd.concat([df_phi, df_psi], axis=1)
    df_phi_psi.to_csv("Phi_Psi.dat", sep="\t", index=False, header=False)
    # Creating weights.dat
    with open("gamd.log") as f:
        lines = f.readlines()
    column_names = lines[2]
    column_names = column_names.replace("#", "")
    column_names = column_names.replace("\n", "")
    column_names = column_names.replace(" ", "")
    column_names = column_names.split(",")
    list_words = ["#"]
    with open("gamd.log") as oldfile, open("data.log", "w") as newfile:
        for line in oldfile:
            if not any(word in line for word in list_words):
                newfile.write(line)
    df = pd.read_csv("data.log", delim_whitespace=True, header=None)
    df.columns = column_names
    df["dV(kcal/mol)"] = (
        df["Boost-Energy-Potential"] + df["Boost-Energy-Dihedral"]
    )
    df["dV(kbT)"] = df["dV(kcal/mol)"] / factor
    df_ = df[["dV(kbT)", "total_nstep", "dV(kcal/mol)"]]
    df_ = df_[::jump]
    df_.to_csv("weights.dat", sep="\t", index=False, header=False)
    os.system("rm -rf data.log")
    print(df_phi_psi.shape)
    print(df_phi.shape)
    print(df_.shape)


def create_data_files(
    jump=10,
    traj="system_final.nc",
    topology="system_final.prmtop",
    T=300,
):

    """

    Extracts data from GaMD log files and saves them as
    weights.dat, Psi.dat and Phi_Psi.dat. gamd.log file
    contains data excluding the initial equilibration MD
    simulation steps but trajectory output file has all
    the trajectories including the initial equilibration
    MD steps.  This part has ben taken care to make the
    data consistent.

    Parameters
    ----------

    jump: int
        Every nth frame to be considered for reweighting

    traj: str
        System's trajectory file

    topology: str
        System's topology file

    T: int
        MD simulation temperature

    """

    # To make data consistent with gamd.log and .nc file
    factor = 0.001987 * T
    with open("md.in") as f:
        lines = f.readlines()
    for i in lines:
        if "nstlim =" in i:
            nstlim_line = i
        if "ntcmd =" in i:
            ntcmd_line = i
        if "ntwx =" in i:
            ntwx_line = i
    x = re.findall(r"\b\d+\b", ntcmd_line)
    ntcmd = int(x[0])
    x = re.findall(r"\b\d+\b", nstlim_line)
    nstlim = int(x[0])
    x = re.findall(r"\b\d+\b", ntwx_line)
    ntwx = int(x[1])
    # From the .nc trajectory files, we will not consider ntcmd trajectories
    leave_frames = int(ntcmd / ntwx)
    no_frames = int(nstlim / ntwx)
    # Recheck conditions
    file = open("gamd.log", "r")
    number_of_lines = 0
    for line in file:
        line = line.strip("\n")
        number_of_lines += 1
    file.close()
    f = open("gamd.log")
    fourth_line = f.readlines()[3]
    if str(ntcmd) in fourth_line:
        datapoints = number_of_lines - 4
    if not str(ntcmd) in fourth_line:
        datapoints = number_of_lines - 3
    print(datapoints == int((nstlim - ntcmd) / ntwx))
    # Creating Psi.dat and Phi_Psi.dat
    traj = md.load(traj, top=topology)
    traj = traj[leave_frames:no_frames:jump]
    phi = md.compute_phi(traj)
    phi = phi[1]  # 0:indices, 1:phi angles
    phi = np.array([math.degrees(i) for i in phi])  # radians to degrees
    psi = md.compute_psi(traj)
    psi = psi[1]  # 0:indices, 1:psi angles
    psi = np.array([math.degrees(i) for i in psi])  # radians to degrees
    df_psi = pd.DataFrame(phi, columns=["Psi"])
    df_psi = df_psi.tail(int(datapoints))
    df_psi.to_csv("Psi.dat", sep="\t", index=False, header=False)
    df_phi = pd.DataFrame(psi, columns=["Phi"])
    df_phi = df_phi.tail(int(datapoints))
    df_phi_psi = pd.concat([df_phi, df_psi], axis=1)
    df_phi_psi.to_csv("Phi_Psi.dat", sep="\t", index=False, header=False)
    # Creating weights.dat
    with open("gamd.log") as f:
        lines = f.readlines()
    column_names = lines[2]
    column_names = column_names.replace("#", "")
    column_names = column_names.replace("\n", "")
    column_names = column_names.replace(" ", "")
    column_names = column_names.split(",")
    list_words = ["#"]
    with open("gamd.log") as oldfile, open("data.log", "w") as newfile:
        for line in oldfile:
            if not any(word in line for word in list_words):
                newfile.write(line)
    df = pd.read_csv("data.log", delim_whitespace=True, header=None)
    df.columns = column_names
    df["dV(kcal/mol)"] = (
        df["Boost-Energy-Potential"] + df["Boost-Energy-Dihedral"]
    )
    df["dV(kbT)"] = df["dV(kcal/mol)"] / factor
    df_ = df[["dV(kbT)", "total_nstep", "dV(kcal/mol)"]]
    df_ = df_[::jump]
    df_.to_csv("weights.dat", sep="\t", index=False, header=False)
    os.system("rm -rf data.log")
    print(df_phi_psi.shape)
    print(df_phi.shape)
    print(df_.shape)


def create_bins(lower_bound, width, upper_bound):

    """

    Creates bin if given the lower and upper bound
    with the wirdth information.

    """

    bins = []
    for low in range(lower_bound, upper_bound, width):
        bins.append([low, low + width])
    return bins


def find_bin(value, bins):

    """

    Finds which value belongs to which bin.

    """

    for i in range(0, len(bins)):
        if bins[i][0] <= value < bins[i][1]:
            return i
    return -1


def reweight_1d(
    binspace=10, n_structures=4, Xdim=[-180, 180], T=300.0, min_prob=0.000001
):

    """

    Reweights boosted potential energies in one-dimension based on
    Maclaurin series expansion to one, two and three degrees.

    Parameters
    ----------

    binspace: int
        Spacing between the bins

    n_structures: int
        Number of structures per bin chosen
        for Weighted Ensemble (WE) simulations

    Xdim: list
        Range of dihedral angles

    T: float
        MD simulation temperature

    min_prob: float
        minimum probability threshold

    """

    beta = 1.0 / (0.001987 * float(T))
    df_Psi = pd.read_csv("Psi.dat", delim_whitespace=True, header=None)
    df_Psi.columns = ["Psi"]
    df_weight = pd.read_csv("weights.dat", delim_whitespace=True, header=None)
    df_weight.columns = ["dV_kBT", "timestep", "dVkcalmol"]
    sum_total = df_Psi.shape[0]
    binsX = np.arange(float(Xdim[0]), (float(Xdim[1]) + binspace), binspace)
    hist, hist_edges = np.histogram(df_Psi[["Psi"]], bins=binsX, weights=None)
    pstarA = [i / sum_total for i in list(hist)]
    bins = create_bins(
        lower_bound=int(Xdim[0]), width=binspace, upper_bound=int(Xdim[1])
    )
    data = df_Psi["Psi"].values.tolist()
    binned_weights = []
    for value in data:
        bin_index = find_bin(value, bins)
        binned_weights.append(bin_index)
    df_index = pd.DataFrame(binned_weights)
    df_index.columns = ["index"]
    df = pd.concat([df_index, df_Psi, df_weight], axis=1)
    dV_c1 = []
    dV_c2 = []
    dV_c3 = []
    dV = []
    for i in range(len(bins)):
        df_i = df.loc[(df["index"] == i)]
        dV_list = df_i["dVkcalmol"].values.tolist()
        if len(dV_list) >= 10:
            dV_c1.append(statistics.mean(dV_list))
            dV_c2.append(
                statistics.mean([i ** 2 for i in dV_list])
                - (statistics.mean(dV_list)) ** 2
            )
            dV_c3.append(
                statistics.mean([i ** 3 for i in dV_list])
                - 3
                * (statistics.mean([i ** 2 for i in dV_list]))
                * (statistics.mean(dV_list))
                + 2 * (statistics.mean(dV_list)) ** 3
            )
        if len(dV_list) < 10:
            dV_c1.append(0)
            dV_c2.append(0)
            dV_c3.append(0)
        dV.append(dV_list)
    c1 = [i * beta for i in dV_c1]
    c2 = [i * ((beta ** 2) / 2) for i in dV_c2]
    c3 = [i * ((beta ** 3) / 6) for i in dV_c3]
    c1 = c1
    c12 = [a + b for a, b in zip(c1, c2)]
    c123 = [a + b for a, b in zip(c12, c3)]
    for i in range(len(c1)):
        if c1[i] >= 700:
            c1[i] = 700
    for i in range(len(c12)):
        if c12[i] >= 700:
            c12[i] = 700
    for i in range(len(c123)):
        if c123[i] >= 700:
            c123[i] = 700
    ensemble_average_c1 = [exp(i) for i in c1]
    ensemble_average_c12 = [exp(i) for i in c12]
    ensemble_average_c123 = [exp(i) for i in c123]
    numerator_c1 = [a * b for a, b in zip(pstarA, ensemble_average_c1)]
    numerator_c12 = [a * b for a, b in zip(pstarA, ensemble_average_c12)]
    numerator_c123 = [a * b for a, b in zip(pstarA, ensemble_average_c123)]
    #### c1
    denominatorc1 = []
    for i in range(len(bins)):
        product_c1 = pstarA[i] * ensemble_average_c1[i]
        denominatorc1.append(product_c1)
    denominator_c1 = sum(denominatorc1)
    pA_c1 = [i / denominator_c1 for i in numerator_c1]
    #### c12
    denominatorc12 = []
    for i in range(len(bins)):
        product_c12 = pstarA[i] * ensemble_average_c12[i]
        denominatorc12.append(product_c12)
    denominator_c12 = sum(denominatorc12)
    pA_c12 = [i / denominator_c12 for i in numerator_c12]
    #### c123
    denominatorc123 = []
    for i in range(len(bins)):
        product_c123 = pstarA[i] * ensemble_average_c123[i]
        denominatorc123.append(product_c123)
    denominator_c123 = sum(denominatorc123)
    pA_c123 = [i / denominator_c123 for i in numerator_c123]
    data_c1 = list(zip(bins, pA_c1))
    data_c12 = list(zip(bins, pA_c12))
    data_c123 = list(zip(bins, pA_c123))
    df_c1 = pd.DataFrame(data_c1, columns=["bins", "pA_c1"])
    df_c12 = pd.DataFrame(data_c12, columns=["bins", "pA_c12"])
    df_c123 = pd.DataFrame(data_c123, columns=["bins", "pA_c123"])
    ####c1
    df_c1.to_csv("c1_1d.txt", header=True, index=None, sep=" ", mode="w")
    with open("c1_1d.txt", "r") as f1, open("pA_c1_1d.txt", "w") as f2:
        for line in f1:
            f2.write(line.replace('"', "").replace("'", ""))
    os.system("rm -rf c1_1d.txt")
    ####c12
    df_c12.to_csv("c12_1d.txt", header=True, index=None, sep=" ", mode="w")
    with open("c12_1d.txt", "r") as f1, open("pA_c12_1d.txt", "w") as f2:
        for line in f1:
            f2.write(line.replace('"', "").replace("'", ""))
    os.system("rm -rf c12_1d.txt")
    ####c123
    df_c123.to_csv("c123_1d.txt", header=True, index=None, sep=" ", mode="w")
    with open("c123_1d.txt", "r") as f1, open("pA_c123_1d.txt", "w") as f2:
        for line in f1:
            f2.write(line.replace('"', "").replace("'", ""))
    os.system("rm -rf c123_1d.txt")
    ####c1_arranged
    df_c1_arranged = df_c1.sort_values(by="pA_c1", ascending=False)
    df_c1_arranged = df_c1_arranged[df_c1_arranged.pA_c1 > min_prob]
    df_c1_arranged.to_csv(
        "c1_arranged_1d.txt", header=True, index=None, sep=" ", mode="w"
    )
    with open("c1_arranged_1d.txt", "r") as f1, open(
        "pA_c1_arranged_1d.txt", "w"
    ) as f2:
        for line in f1:
            f2.write(line.replace('"', "").replace("'", ""))
    os.system("rm -rf c1_arranged_1d.txt")
    ####c12_arranged
    df_c12_arranged = df_c12.sort_values(by="pA_c12", ascending=False)
    df_c12_arranged = df_c12_arranged[df_c12_arranged.pA_c12 > min_prob]
    df_c12_arranged.to_csv(
        "c12_arranged_1d.txt", header=True, index=None, sep=" ", mode="w"
    )
    with open("c12_arranged_1d.txt", "r") as f1, open(
        "pA_c12_arranged_1d.txt", "w"
    ) as f2:
        for line in f1:
            f2.write(line.replace('"', "").replace("'", ""))
    os.system("rm -rf c12_arranged_1d.txt")
    ####c123_arranged
    df_c123_arranged = df_c123.sort_values(by="pA_c123", ascending=False)
    df_c123_arranged = df_c123_arranged[df_c123_arranged.pA_c123 > min_prob]
    df_c123_arranged.to_csv(
        "c123_arranged_1d.txt", header=True, index=None, sep=" ", mode="w"
    )
    with open("c123_arranged_1d.txt", "r") as f1, open(
        "pA_c123_arranged_1d.txt", "w"
    ) as f2:
        for line in f1:
            f2.write(line.replace('"', "").replace("'", ""))
    os.system("rm -rf c123_arranged_1d.txt")
    ####c1_arranged
    df_c1_arranged["index"] = df_c1_arranged.index
    index_list_c1 = df_c1_arranged["index"].tolist()
    df["frame_index"] = df.index
    df_frame_index = df[["frame_index", "index"]]
    frame_indices_c1 = []
    index_indces_c1 = []
    for i in index_list_c1:
        df_index_list_c1 = df_frame_index.loc[df_frame_index["index"] == i]
        frame_c1 = df_index_list_c1["frame_index"].tolist()
        frame_indices_c1.append(frame_c1)
        index_c1 = [i] * len(frame_c1)
        index_indces_c1.append(index_c1)
    frame_indices_c1 = [item for elem in frame_indices_c1 for item in elem]
    index_indces_c1 = [item for elem in index_indces_c1 for item in elem]
    df_c1_frame = pd.DataFrame(frame_indices_c1, columns=["frame_index"])
    df_c1_index = pd.DataFrame(index_indces_c1, columns=["index"])
    df_c1_frame_index = pd.concat([df_c1_frame, df_c1_index], axis=1)
    df_c1_frame_index = df_c1_frame_index.groupby("index").filter(
        lambda x: len(x) >= 10
    )
    df_c1_frame_index.to_csv(
        "c1_frame_index_1d.txt", header=True, index=None, sep=" ", mode="w"
    )
    with open("c1_frame_index_1d.txt", "r") as f1, open(
        "c1_frame_1d.txt", "w"
    ) as f2:
        for line in f1:
            f2.write(line.replace('"', "").replace("'", ""))
    os.system("rm -rf c1_frame_index_1d.txt")
    ####c12_arranged
    df_c12_arranged["index"] = df_c12_arranged.index
    index_list_c12 = df_c12_arranged["index"].tolist()
    df["frame_index"] = df.index
    df_frame_index = df[["frame_index", "index"]]
    frame_indices_c12 = []
    index_indces_c12 = []
    for i in index_list_c12:
        df_index_list_c12 = df_frame_index.loc[df_frame_index["index"] == i]
        frame_c12 = df_index_list_c12["frame_index"].tolist()
        frame_indices_c12.append(frame_c12)
        index_c12 = [i] * len(frame_c12)
        index_indces_c12.append(index_c12)
    frame_indices_c12 = [item for elem in frame_indices_c12 for item in elem]
    index_indces_c12 = [item for elem in index_indces_c12 for item in elem]
    df_c12_frame = pd.DataFrame(frame_indices_c12, columns=["frame_index"])
    df_c12_index = pd.DataFrame(index_indces_c12, columns=["index"])
    df_c12_frame_index = pd.concat([df_c12_frame, df_c12_index], axis=1)
    df_c12_frame_index = df_c12_frame_index.groupby("index").filter(
        lambda x: len(x) >= 10
    )
    df_c12_frame_index.to_csv(
        "c12_frame_index_1d.txt", header=True, index=None, sep=" ", mode="w"
    )
    with open("c12_frame_index_1d.txt", "r") as f1, open(
        "c12_frame_1d.txt", "w"
    ) as f2:
        for line in f1:
            f2.write(line.replace('"', "").replace("'", ""))
    os.system("rm -rf c12_frame_index_1d.txt")
    ####c123_arranged
    df_c123_arranged["index"] = df_c123_arranged.index
    index_list_c123 = df_c123_arranged["index"].tolist()
    df["frame_index"] = df.index
    df_frame_index = df[["frame_index", "index"]]
    frame_indices_c123 = []
    index_indces_c123 = []
    for i in index_list_c123:
        df_index_list_c123 = df_frame_index.loc[df_frame_index["index"] == i]
        frame_c123 = df_index_list_c123["frame_index"].tolist()
        frame_indices_c123.append(frame_c123)
        index_c123 = [i] * len(frame_c123)
        index_indces_c123.append(index_c123)
    frame_indices_c123 = [item for elem in frame_indices_c123 for item in elem]
    index_indces_c123 = [item for elem in index_indces_c123 for item in elem]
    df_c123_frame = pd.DataFrame(frame_indices_c123, columns=["frame_index"])
    df_c123_index = pd.DataFrame(index_indces_c123, columns=["index"])
    df_c123_frame_index = pd.concat([df_c123_frame, df_c123_index], axis=1)
    df_c123_frame_index = df_c123_frame_index.groupby("index").filter(
        lambda x: len(x) >= 10
    )
    df_c123_frame_index.to_csv(
        "c123_frame_index_1d.txt", header=True, index=None, sep=" ", mode="w"
    )
    with open("c123_frame_index_1d.txt", "r") as f1, open(
        "c123_frame_1d.txt", "w"
    ) as f2:
        for line in f1:
            f2.write(line.replace('"', "").replace("'", ""))
    os.system("rm -rf c123_frame_index_1d.txt")
    ####c1
    indices_c1_1d = df_c1_frame_index["index"].unique()
    frames_c1 = []
    for i in indices_c1_1d:
        x = df_c1_frame_index.loc[df_c1_frame_index["index"] == i]
        y = x["frame_index"].values.tolist()
        z = random.sample(y, n_structures)
        frames_c1.append(z)
    frames_c1_1d = [item for elem in frames_c1 for item in elem]
    with open("frames_c1_1d.pickle", "wb") as f:
        pk.dump(frames_c1_1d, f)
    with open("indices_c1_1d.pickle", "wb") as f:
        pk.dump(indices_c1_1d, f)
    ####c12
    indices_c12_1d = df_c12_frame_index["index"].unique()
    frames_c12 = []
    for i in indices_c12_1d:
        x = df_c12_frame_index.loc[df_c12_frame_index["index"] == i]
        y = x["frame_index"].values.tolist()
        z = random.sample(y, n_structures)
        frames_c12.append(z)
    frames_c12_1d = [item for elem in frames_c12 for item in elem]
    with open("frames_c12_1d.pickle", "wb") as f:
        pk.dump(frames_c12_1d, f)
    with open("indices_c12_1d.pickle", "wb") as f:
        pk.dump(indices_c12_1d, f)
    ####c123
    indices_c123_1d = df_c123_frame_index["index"].unique()
    frames_c123 = []
    for i in indices_c123_1d:
        x = df_c123_frame_index.loc[df_c123_frame_index["index"] == i]
        y = x["frame_index"].values.tolist()
        z = random.sample(y, n_structures)
        frames_c123.append(z)
    frames_c123_1d = [item for elem in frames_c123 for item in elem]
    with open("frames_c123_1d.pickle", "wb") as f:
        pk.dump(frames_c123_1d, f)
    with open("indices_c123_1d.pickle", "wb") as f:
        pk.dump(indices_c123_1d, f)
    # saving probabilities for  each selected frame
    ####c1
    prob_c1_1d_list = []
    for i in indices_c1_1d:
        prob_c1_1d_list.append(df_c1["pA_c1"][i])
    prob_c1_1d_list = list(
        itertools.chain.from_iterable(
            itertools.repeat(x, n_structures) for x in prob_c1_1d_list
        )
    )
    prob_c1_1d_list = [x / n_structures for x in prob_c1_1d_list]
    with open("prob_c1_1d_list.pickle", "wb") as f:
        pk.dump(prob_c1_1d_list, f)
    ####c12
    prob_c12_1d_list = []
    for i in indices_c12_1d:
        prob_c12_1d_list.append(df_c12["pA_c12"][i])
    prob_c12_1d_list = list(
        itertools.chain.from_iterable(
            itertools.repeat(x, n_structures) for x in prob_c12_1d_list
        )
    )
    prob_c12_1d_list = [x / n_structures for x in prob_c12_1d_list]
    with open("prob_c12_1d_list.pickle", "wb") as f:
        pk.dump(prob_c12_1d_list, f)
    ####c123
    prob_c123_1d_list = []
    for i in indices_c123_1d:
        prob_c123_1d_list.append(df_c123["pA_c123"][i])
    prob_c123_1d_list = list(
        itertools.chain.from_iterable(
            itertools.repeat(x, n_structures) for x in prob_c123_1d_list
        )
    )
    prob_c123_1d_list = [x / n_structures for x in prob_c123_1d_list]
    with open("prob_c123_1d_list.pickle", "wb") as f:
        pk.dump(prob_c123_1d_list, f)
    ref_df_1d = pd.DataFrame(bins, columns=["dim0", "dim1"])
    ref_df_1d["bins"] = ref_df_1d.agg(
        lambda x: f"[{x['dim0']} , {x['dim1']}]", axis=1
    )
    ref_df_1d = ref_df_1d[["bins"]]
    index_ref_1d = []
    for i in range(len(bins)):
        index_ref_1d.append(i)
    index_ref_df_1d = pd.DataFrame(index_ref_1d, columns=["index"])
    df_ref_1d = pd.concat([ref_df_1d, index_ref_df_1d], axis=1)
    df_ref_1d.to_csv("ref_1d.txt", header=True, index=None, sep=" ", mode="w")
    df.to_csv("df_1d.csv", index=False)
    os.system("rm -rf __pycache__")
    print("Successfully Completed Reweighing")


def reweight_2d(
    binspace=10,
    n_structures=4,
    Xdim=[-180, 180],
    Ydim=[-180, 180],
    T=300.0,
    min_prob=0.000001,
):

    """

    Reweights boosted potential energies in two-dimensions
    based on Maclaurin series expansion to one, two and
    three  degrees.

    Parameters
    ----------

    binspace: int
        Spacing between the bins

    n_structures: int
        Number of structures per bin chosen
        for Weighted Ensemble (WE) simulations

    Xdim: list
        Range of dihedral angles (1st dimension)

    Ydim: list
        Range of dihedral angles (2nd dimension)

    T: float
        MD simulation temperature

    min_prob: float
        minimum probability threshold

    """

    beta = 1.0 / (0.001987 * float(T))
    df_Phi_Psi = pd.read_csv("Phi_Psi.dat", delim_whitespace=True, header=None)
    df_Phi_Psi.columns = ["Phi", "Psi"]
    df_weight = pd.read_csv("weights.dat", delim_whitespace=True, header=None)
    df_weight.columns = ["dV_kBT", "timestep", "dVkcalmol"]
    sum_total = df_Phi_Psi.shape[0]
    binsX = np.arange(float(Xdim[0]), (float(Xdim[1]) + binspace), binspace)
    binsY = np.arange(float(Ydim[0]), (float(Ydim[1]) + binspace), binspace)
    hist2D, hist_edgesX, hist_edgesY = np.histogram2d(
        df_Phi_Psi["Phi"].values.tolist(),
        df_Phi_Psi["Psi"].values.tolist(),
        bins=(binsX, binsY),
        weights=None,
    )
    pstarA_2D = [i / sum_total for i in list(hist2D)]
    bins_tuple_X = create_bins(
        lower_bound=int(Xdim[0]), width=binspace, upper_bound=int(Xdim[1])
    )
    bins_tuple_Y = create_bins(
        lower_bound=int(Ydim[0]), width=binspace, upper_bound=int(Ydim[1])
    )
    bins = []
    for i in range(len(bins_tuple_X)):
        for j in range(len(bins_tuple_Y)):
            bins.append([bins_tuple_X[i], bins_tuple_Y[j]])
    pstarA = [item for elem in pstarA_2D for item in elem]
    hist = [item for elem in hist2D for item in elem]
    hist = [int(i) for i in hist]
    data_X = df_Phi_Psi["Phi"].values.tolist()
    binned_weights_X = []
    for value in data_X:
        bin_index_X = find_bin(value, bins_tuple_X)
        binned_weights_X.append(bin_index_X)
    data_Y = df_Phi_Psi["Psi"].values.tolist()
    binned_weights_Y = []
    for value in data_Y:
        bin_index_Y = find_bin(value, bins_tuple_Y)
        binned_weights_Y.append(bin_index_Y)
    binned_weights_2D = []
    for i in range(len(binned_weights_X)):
        binned_weights_2D.append([binned_weights_X[i], binned_weights_Y[i]])
    binned_weights = []
    for i in range(len(binned_weights_2D)):
        binned_weights.append(
            (binned_weights_2D[i][0] * len(bins_tuple_Y))
            + (binned_weights_2D[i][1] + 1)
        )
    df_index = pd.DataFrame(binned_weights)
    df_index.columns = ["index"]
    df_index["index"] = df_index["index"] - 1
    df = pd.concat([df_index, df_Phi_Psi, df_weight], axis=1)
    dV_c1 = []
    dV_c2 = []
    dV_c3 = []
    dV = []
    for i in range(len(bins)):
        df_i = df.loc[(df["index"] == i)]
        dV_list = df_i["dVkcalmol"].values.tolist()
        if len(dV_list) >= 10:
            dV_c1.append(statistics.mean(dV_list))
            dV_c2.append(
                statistics.mean([i ** 2 for i in dV_list])
                - (statistics.mean(dV_list)) ** 2
            )
            dV_c3.append(
                statistics.mean([i ** 3 for i in dV_list])
                - 3
                * (statistics.mean([i ** 2 for i in dV_list]))
                * (statistics.mean(dV_list))
                + 2 * (statistics.mean(dV_list)) ** 3
            )
        if len(dV_list) < 10:
            dV_c1.append(0)
            dV_c2.append(0)
            dV_c3.append(0)
        dV.append(dV_list)
    c1 = [i * beta for i in dV_c1]
    c2 = [i * ((beta ** 2) / 2) for i in dV_c2]
    c3 = [i * ((beta ** 3) / 6) for i in dV_c3]
    c1 = c1
    c12 = [a + b for a, b in zip(c1, c2)]
    c123 = [a + b for a, b in zip(c12, c3)]
    for i in range(len(c1)):
        if c1[i] >= 700:
            c1[i] = 700
    for i in range(len(c12)):
        if c12[i] >= 700:
            c12[i] = 700
    for i in range(len(c123)):
        if c123[i] >= 700:
            c123[i] = 700
    ensemble_average_c1 = [exp(i) for i in c1]
    ensemble_average_c12 = [exp(i) for i in c12]
    ensemble_average_c123 = [exp(i) for i in c123]
    numerator_c1 = [a * b for a, b in zip(pstarA, ensemble_average_c1)]
    numerator_c12 = [a * b for a, b in zip(pstarA, ensemble_average_c12)]
    numerator_c123 = [a * b for a, b in zip(pstarA, ensemble_average_c123)]
    #### c1
    denominatorc1 = []
    for i in range(len(bins)):
        product_c1 = pstarA[i] * ensemble_average_c1[i]
        denominatorc1.append(product_c1)
    denominator_c1 = sum(denominatorc1)
    pA_c1 = [i / denominator_c1 for i in numerator_c1]
    #### c12
    denominatorc12 = []
    for i in range(len(bins)):
        product_c12 = pstarA[i] * ensemble_average_c12[i]
        denominatorc12.append(product_c12)
    denominator_c12 = sum(denominatorc12)
    pA_c12 = [i / denominator_c12 for i in numerator_c12]
    #### c123
    denominatorc123 = []
    for i in range(len(bins)):
        product_c123 = pstarA[i] * ensemble_average_c123[i]
        denominatorc123.append(product_c123)
    denominator_c123 = sum(denominatorc123)
    pA_c123 = [i / denominator_c123 for i in numerator_c123]
    data_c1 = list(zip(bins, pA_c1))
    data_c12 = list(zip(bins, pA_c12))
    data_c123 = list(zip(bins, pA_c123))
    df_c1 = pd.DataFrame(data_c1, columns=["bins", "pA_c1"])
    df_c12 = pd.DataFrame(data_c12, columns=["bins", "pA_c12"])
    df_c123 = pd.DataFrame(data_c123, columns=["bins", "pA_c123"])
    df_c1.to_csv("c1_2d.txt", header=True, index=None, sep=" ", mode="w")
    with open("c1_2d.txt", "r") as f1, open("pA_c1_2d.txt", "w") as f2:
        for line in f1:
            f2.write(line.replace('"', "").replace("'", ""))
    os.system("rm -rf c1_2d.txt")
    ####c12
    df_c12.to_csv("c12_2d.txt", header=True, index=None, sep=" ", mode="w")
    with open("c12_2d.txt", "r") as f1, open("pA_c12_2d.txt", "w") as f2:
        for line in f1:
            f2.write(line.replace('"', "").replace("'", ""))
    os.system("rm -rf c12_2d.txt")
    ####c123
    df_c123.to_csv("c123_2d.txt", header=True, index=None, sep=" ", mode="w")
    with open("c123_2d.txt", "r") as f1, open("pA_c123_2d.txt", "w") as f2:
        for line in f1:
            f2.write(line.replace('"', "").replace("'", ""))
    os.system("rm -rf c123_2d.txt")
    ####c1_arranged
    df_c1_arranged = df_c1.sort_values(by="pA_c1", ascending=False)
    df_c1_arranged = df_c1_arranged[df_c1_arranged.pA_c1 > min_prob]
    df_c1_arranged.to_csv(
        "c1_arranged_2d.txt", header=True, index=None, sep=" ", mode="w"
    )
    with open("c1_arranged_2d.txt", "r") as f1, open(
        "pA_c1_arranged_2d.txt", "w"
    ) as f2:
        for line in f1:
            f2.write(line.replace('"', "").replace("'", ""))
    os.system("rm -rf c1_arranged_2d.txt")
    ####c12_arranged
    df_c12_arranged = df_c12.sort_values(by="pA_c12", ascending=False)
    df_c12_arranged = df_c12_arranged[df_c12_arranged.pA_c12 > min_prob]
    df_c12_arranged.to_csv(
        "c12_arranged_2d.txt", header=True, index=None, sep=" ", mode="w"
    )
    with open("c12_arranged_2d.txt", "r") as f1, open(
        "pA_c12_arranged_2d.txt", "w"
    ) as f2:
        for line in f1:
            f2.write(line.replace('"', "").replace("'", ""))
    os.system("rm -rf c12_arranged_2d.txt")
    ####c123_arranged
    df_c123_arranged = df_c123.sort_values(by="pA_c123", ascending=False)
    df_c123_arranged = df_c123_arranged[df_c123_arranged.pA_c123 > min_prob]
    df_c123_arranged.to_csv(
        "c123_arranged_2d.txt", header=True, index=None, sep=" ", mode="w"
    )
    with open("c123_arranged_2d.txt", "r") as f1, open(
        "pA_c123_arranged_2d.txt", "w"
    ) as f2:
        for line in f1:
            f2.write(line.replace('"', "").replace("'", ""))
    os.system("rm -rf c123_arranged_2d.txt")
    ####c1_arranged
    df_c1_arranged["index"] = df_c1_arranged.index
    index_list_c1 = df_c1_arranged["index"].tolist()
    df["frame_index"] = df.index
    df_frame_index = df[["frame_index", "index"]]
    frame_indices_c1 = []
    index_indces_c1 = []
    for i in index_list_c1:
        df_index_list_c1 = df_frame_index.loc[df_frame_index["index"] == i]
        frame_c1 = df_index_list_c1["frame_index"].tolist()
        frame_indices_c1.append(frame_c1)
        index_c1 = [i] * len(frame_c1)
        index_indces_c1.append(index_c1)
    frame_indices_c1 = [item for elem in frame_indices_c1 for item in elem]
    index_indces_c1 = [item for elem in index_indces_c1 for item in elem]
    df_c1_frame = pd.DataFrame(frame_indices_c1, columns=["frame_index"])
    df_c1_index = pd.DataFrame(index_indces_c1, columns=["index"])
    df_c1_frame_index = pd.concat([df_c1_frame, df_c1_index], axis=1)
    df_c1_frame_index = df_c1_frame_index.groupby("index").filter(
        lambda x: len(x) >= 10
    )
    df_c1_frame_index.to_csv(
        "c1_frame_index_2d.txt", header=True, index=None, sep=" ", mode="w"
    )
    with open("c1_frame_index_2d.txt", "r") as f1, open(
        "c1_frame_2d.txt", "w"
    ) as f2:
        for line in f1:
            f2.write(line.replace('"', "").replace("'", ""))
    os.system("rm -rf c1_frame_index_2d.txt")
    ####c12_arranged
    df_c12_arranged["index"] = df_c12_arranged.index
    index_list_c12 = df_c12_arranged["index"].tolist()
    df["frame_index"] = df.index
    df_frame_index = df[["frame_index", "index"]]
    frame_indices_c12 = []
    index_indces_c12 = []
    for i in index_list_c12:
        df_index_list_c12 = df_frame_index.loc[df_frame_index["index"] == i]
        frame_c12 = df_index_list_c12["frame_index"].tolist()
        frame_indices_c12.append(frame_c12)
        index_c12 = [i] * len(frame_c12)
        index_indces_c12.append(index_c12)
    frame_indices_c12 = [item for elem in frame_indices_c12 for item in elem]
    index_indces_c12 = [item for elem in index_indces_c12 for item in elem]
    df_c12_frame = pd.DataFrame(frame_indices_c12, columns=["frame_index"])
    df_c12_index = pd.DataFrame(index_indces_c12, columns=["index"])
    df_c12_frame_index = pd.concat([df_c12_frame, df_c12_index], axis=1)
    df_c12_frame_index = df_c12_frame_index.groupby("index").filter(
        lambda x: len(x) >= 10
    )
    df_c12_frame_index.to_csv(
        "c12_frame_index_2d.txt", header=True, index=None, sep=" ", mode="w"
    )
    with open("c12_frame_index_2d.txt", "r") as f1, open(
        "c12_frame_2d.txt", "w"
    ) as f2:
        for line in f1:
            f2.write(line.replace('"', "").replace("'", ""))
    os.system("rm -rf c12_frame_index_2d.txt")
    ####c123_arranged
    df_c123_arranged["index"] = df_c123_arranged.index
    df_c123_arranged["index"] = df_c123_arranged.index
    index_list_c123 = df_c123_arranged["index"].tolist()
    df["frame_index"] = df.index
    df_frame_index = df[["frame_index", "index"]]
    frame_indices_c123 = []
    index_indces_c123 = []
    for i in index_list_c123:
        df_index_list_c123 = df_frame_index.loc[df_frame_index["index"] == i]
        frame_c123 = df_index_list_c123["frame_index"].tolist()
        frame_indices_c123.append(frame_c123)
        index_c123 = [i] * len(frame_c123)
        index_indces_c123.append(index_c123)
    frame_indices_c123 = [item for elem in frame_indices_c123 for item in elem]
    index_indces_c123 = [item for elem in index_indces_c123 for item in elem]
    df_c123_frame = pd.DataFrame(frame_indices_c123, columns=["frame_index"])
    df_c123_index = pd.DataFrame(index_indces_c123, columns=["index"])
    df_c123_frame_index = pd.concat([df_c123_frame, df_c123_index], axis=1)
    df_c123_frame_index = df_c123_frame_index.groupby("index").filter(
        lambda x: len(x) >= 10
    )
    df_c123_frame_index.to_csv(
        "c123_frame_index_2d.txt", header=True, index=None, sep=" ", mode="w"
    )
    with open("c123_frame_index_2d.txt", "r") as f1, open(
        "c123_frame_2d.txt", "w"
    ) as f2:
        for line in f1:
            f2.write(line.replace('"', "").replace("'", ""))
    os.system("rm -rf c123_frame_index_2d.txt")
    ####c1
    indices_c1_2d = df_c1_frame_index["index"].unique()
    frames_c1 = []
    for i in indices_c1_2d:
        x = df_c1_frame_index.loc[df_c1_frame_index["index"] == i]
        y = x["frame_index"].values.tolist()
        z = random.sample(y, n_structures)
        frames_c1.append(z)
    frames_c1_2d = [item for elem in frames_c1 for item in elem]
    with open("frames_c1_2d.pickle", "wb") as f:
        pk.dump(frames_c1_2d, f)
    with open("indices_c1_2d.pickle", "wb") as f:
        pk.dump(indices_c1_2d, f)
    ####c12
    indices_c12_2d = df_c12_frame_index["index"].unique()
    frames_c12 = []
    for i in indices_c12_2d:
        x = df_c12_frame_index.loc[df_c12_frame_index["index"] == i]
        y = x["frame_index"].values.tolist()
        z = random.sample(y, n_structures)
        frames_c12.append(z)
    frames_c12_2d = [item for elem in frames_c12 for item in elem]
    with open("frames_c12_2d.pickle", "wb") as f:
        pk.dump(frames_c12_2d, f)
    with open("indices_c12_2d.pickle", "wb") as f:
        pk.dump(indices_c12_2d, f)
    ####c123
    indices_c123_2d = df_c123_frame_index["index"].unique()
    frames_c123 = []
    for i in indices_c123_2d:
        x = df_c123_frame_index.loc[df_c123_frame_index["index"] == i]
        y = x["frame_index"].values.tolist()
        z = random.sample(y, n_structures)
        frames_c123.append(z)
    frames_c123_2d = [item for elem in frames_c123 for item in elem]
    with open("frames_c123_2d.pickle", "wb") as f:
        pk.dump(frames_c123_2d, f)
    with open("indices_c123_2d.pickle", "wb") as f:
        pk.dump(indices_c123_2d, f)
    ##saving probabilities for  each selected frame
    ####c1
    prob_c1_2d_list = []
    for i in indices_c1_2d:
        prob_c1_2d_list.append(df_c1["pA_c1"][i])
    prob_c1_2d_list = list(
        itertools.chain.from_iterable(
            itertools.repeat(x, n_structures) for x in prob_c1_2d_list
        )
    )
    prob_c1_2d_list = [x / n_structures for x in prob_c1_2d_list]
    with open("prob_c1_2d_list.pickle", "wb") as f:
        pk.dump(prob_c1_2d_list, f)
    ####c12
    prob_c12_2d_list = []
    for i in indices_c12_2d:
        prob_c12_2d_list.append(df_c12["pA_c12"][i])
    prob_c12_2d_list = list(
        itertools.chain.from_iterable(
            itertools.repeat(x, n_structures) for x in prob_c12_2d_list
        )
    )
    prob_c12_2d_list = [x / n_structures for x in prob_c12_2d_list]
    with open("prob_c12_2d_list.pickle", "wb") as f:
        pk.dump(prob_c12_2d_list, f)
    ####c123
    prob_c123_2d_list = []
    for i in indices_c123_2d:
        prob_c123_2d_list.append(df_c123["pA_c123"][i])
    prob_c123_2d_list = list(
        itertools.chain.from_iterable(
            itertools.repeat(x, n_structures) for x in prob_c123_2d_list
        )
    )
    prob_c123_2d_list = [x / n_structures for x in prob_c123_2d_list]
    with open("prob_c123_2d_list.pickle", "wb") as f:
        pk.dump(prob_c123_2d_list, f)
    ref_df_2d = pd.DataFrame(bins, columns=["binsX", "binsY"])
    ref_df_2d["XY"] = ref_df_2d.agg(
        lambda x: f"{x['binsX']} , {x['binsX']}", axis=1
    )
    ref_df_2d = ref_df_2d[["XY"]]
    index_ref_2d = []
    for i in range(len(bins_tuple_X) * len(bins_tuple_Y)):
        index_ref_2d.append(i)
    index_ref_df_2d = pd.DataFrame(index_ref_2d, columns=["index"])
    df_ref_2d = pd.concat([ref_df_2d, index_ref_df_2d], axis=1)
    df_ref_2d.to_csv("ref_2d.txt", header=True, index=None, sep=" ", mode="w")
    df.to_csv("df_2d.csv", index=False)
    os.system("rm -rf __pycache__")
    print("Successfully Completed Reweighing")


def arrange_files():

    """

    Creates directories and move files to appropriate folders.

    """

    cwd = os.getcwd()
    os.system("rm -rf txt_csv_files")
    os.system("rm -rf dat_files")
    os.system("rm -rf pickle_files")
    os.system("rm -rf system_files")
    os.system("mkdir txt_csv_files")
    os.system("mkdir dat_files")
    os.system("mkdir pickle_files")
    os.system("mkdir system_files")
    shutil.move(
        cwd + "/" + "c1_frame_1d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "c1_frame_1d.txt",
    )
    shutil.move(
        cwd + "/" + "c12_frame_1d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "c12_frame_1d.txt",
    )
    shutil.move(
        cwd + "/" + "c123_frame_1d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "c123_frame_1d.txt",
    )
    shutil.move(
        cwd + "/" + "c1_frame_2d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "c1_frame_2d.txt",
    )
    shutil.move(
        cwd + "/" + "c12_frame_2d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "c12_frame_2d.txt",
    )
    shutil.move(
        cwd + "/" + "c123_frame_2d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "c123_frame_2d.txt",
    )
    shutil.move(
        cwd + "/" + "pA_c1_1d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "pA_c1_1d.txt",
    )
    shutil.move(
        cwd + "/" + "pA_c12_1d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "pA_c12_1d.txt",
    )
    shutil.move(
        cwd + "/" + "pA_c123_1d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "pA_c123_1d.txt",
    )
    shutil.move(
        cwd + "/" + "pA_c1_2d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "pA_c1_2d.txt",
    )
    shutil.move(
        cwd + "/" + "pA_c12_2d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "pA_c12_2d.txt",
    )
    shutil.move(
        cwd + "/" + "pA_c123_2d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "pA_c123_2d.txt",
    )
    shutil.move(
        cwd + "/" + "pA_c1_arranged_1d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "pA_c1_arranged_1d.txt",
    )
    shutil.move(
        cwd + "/" + "pA_c12_arranged_1d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "pA_c12_arranged_1d.txt",
    )
    shutil.move(
        cwd + "/" + "pA_c123_arranged_1d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "pA_c123_arranged_1d.txt",
    )
    shutil.move(
        cwd + "/" + "pA_c1_arranged_2d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "pA_c1_arranged_2d.txt",
    )
    shutil.move(
        cwd + "/" + "pA_c12_arranged_2d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "pA_c12_arranged_2d.txt",
    )
    shutil.move(
        cwd + "/" + "pA_c123_arranged_2d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "pA_c123_arranged_2d.txt",
    )
    shutil.move(
        cwd + "/" + "ref_1d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "ref_1d.txt",
    )
    shutil.move(
        cwd + "/" + "ref_2d.txt",
        cwd + "/" + "txt_csv_files" + "/" + "ref_2d.txt",
    )
    shutil.move(
        cwd + "/" + "df_1d.csv",
        cwd + "/" + "txt_csv_files" + "/" + "df_1d.csv",
    )
    shutil.move(
        cwd + "/" + "df_2d.csv",
        cwd + "/" + "txt_csv_files" + "/" + "df_2d.csv",
    )
    shutil.move(
        cwd + "/" + "weights.dat",
        cwd + "/" + "dat_files" + "/" + "weights.dat",
    )
    shutil.move(
        cwd + "/" + "Psi.dat", cwd + "/" + "dat_files" + "/" + "Psi.dat"
    )
    shutil.move(
        cwd + "/" + "Phi_Psi.dat",
        cwd + "/" + "dat_files" + "/" + "Phi_Psi.dat",
    )
    shutil.move(
        cwd + "/" + "prob_c1_1d_list.pickle",
        cwd + "/" + "pickle_files" + "/" + "prob_c1_1d_list.pickle",
    )
    shutil.move(
        cwd + "/" + "prob_c12_1d_list.pickle",
        cwd + "/" + "pickle_files" + "/" + "prob_c12_1d_list.pickle",
    )
    shutil.move(
        cwd + "/" + "prob_c123_1d_list.pickle",
        cwd + "/" + "pickle_files" + "/" + "prob_c123_1d_list.pickle",
    )
    shutil.move(
        cwd + "/" + "prob_c1_2d_list.pickle",
        cwd + "/" + "pickle_files" + "/" + "prob_c1_2d_list.pickle",
    )
    shutil.move(
        cwd + "/" + "prob_c12_2d_list.pickle",
        cwd + "/" + "pickle_files" + "/" + "prob_c12_2d_list.pickle",
    )
    shutil.move(
        cwd + "/" + "prob_c123_2d_list.pickle",
        cwd + "/" + "pickle_files" + "/" + "prob_c123_2d_list.pickle",
    )
    shutil.move(
        cwd + "/" + "frames_c1_1d.pickle",
        cwd + "/" + "pickle_files" + "/" + "frames_c1_1d.pickle",
    )
    shutil.move(
        cwd + "/" + "frames_c12_1d.pickle",
        cwd + "/" + "pickle_files" + "/" + "frames_c12_1d.pickle",
    )
    shutil.move(
        cwd + "/" + "frames_c123_1d.pickle",
        cwd + "/" + "pickle_files" + "/" + "frames_c123_1d.pickle",
    )
    shutil.move(
        cwd + "/" + "frames_c1_2d.pickle",
        cwd + "/" + "pickle_files" + "/" + "frames_c1_2d.pickle",
    )
    shutil.move(
        cwd + "/" + "frames_c12_2d.pickle",
        cwd + "/" + "pickle_files" + "/" + "frames_c12_2d.pickle",
    )
    shutil.move(
        cwd + "/" + "frames_c123_2d.pickle",
        cwd + "/" + "pickle_files" + "/" + "frames_c123_2d.pickle",
    )
    shutil.move(
        cwd + "/" + "indices_c1_1d.pickle",
        cwd + "/" + "pickle_files" + "/" + "indices_c1_1d.pickle",
    )
    shutil.move(
        cwd + "/" + "indices_c12_1d.pickle",
        cwd + "/" + "pickle_files" + "/" + "indices_c12_1d.pickle",
    )
    shutil.move(
        cwd + "/" + "indices_c123_1d.pickle",
        cwd + "/" + "pickle_files" + "/" + "indices_c123_1d.pickle",
    )
    shutil.move(
        cwd + "/" + "indices_c1_2d.pickle",
        cwd + "/" + "pickle_files" + "/" + "indices_c1_2d.pickle",
    )
    shutil.move(
        cwd + "/" + "indices_c12_2d.pickle",
        cwd + "/" + "pickle_files" + "/" + "indices_c12_2d.pickle",
    )
    shutil.move(
        cwd + "/" + "indices_c123_2d.pickle",
        cwd + "/" + "pickle_files" + "/" + "indices_c123_2d.pickle",
    )
    shutil.move(
        cwd + "/" + "system_final.inpcrd",
        cwd + "/" + "system_files" + "/" + "system_final.inpcrd",
    )
    shutil.move(
        cwd + "/" + "system_final.nc",
        cwd + "/" + "system_files" + "/" + "system_final.nc",
    )
    shutil.move(
        cwd + "/" + "system_final.out",
        cwd + "/" + "system_files" + "/" + "system_final.out",
    )
    shutil.move(
        cwd + "/" + "system_final.prmtop",
        cwd + "/" + "system_files" + "/" + "system_final.prmtop",
    )
    shutil.move(
        cwd + "/" + "system_final.rst",
        cwd + "/" + "system_files" + "/" + "system_final.rst",
    )
    shutil.move(
        cwd + "/" + "gamd.log", cwd + "/" + "system_files" + "/" + "gamd.log"
    )
    shutil.move(
        cwd + "/" + "md.in", cwd + "/" + "system_files" + "/" + "md.in"
    )
    shutil.move(
        cwd + "/" + "mdinfo", cwd + "/" + "system_files" + "/" + "mdinfo"
    )
    shutil.move(
        cwd + "/" + "gamd-restart.dat",
        cwd + "/" + "system_files" + "/" + "gamd-restart.dat",
    )


def run_reweigh():

    """
    Runs reweighing calculations systematically
    in the simulation folder.

    """

    dir_list = [
        "dihedral_threshold_lower",
        "dihedral_threshold_upper",
        "dual_threshold_lower",
        "dual_threshold_upper",
        "total_threshold_lower",
        "total_threshold_upper",
    ]
    cwd = os.getcwd()
    target_dir = cwd + "/" + "gamd_simulations" + "/"
    # run reweighting and analysis in each of the simulation folder
    for i in dir_list:
        os.chdir(target_dir + i)
        create_data_files()
        reweight_1d()
        reweight_2d()
        arrange_files()
    os.chdir(cwd)


def clean_for_analysis():

    """

    Rstructures the entire filetree to start reweighing
    analysis again. Used only when we want to run the analysis
    again.

    """

    dir_list = [
        "dihedral_threshold_lower",
        "dihedral_threshold_upper",
        "dual_threshold_lower",
        "dual_threshold_upper",
        "total_threshold_lower",
        "total_threshold_upper",
    ]
    cwd = os.getcwd()
    source_dir = cwd + "/"
    target_dir = cwd + "/" + "gamd_simulations" + "/"
    for i in dir_list:
        os.chdir(target_dir + i)
        os.system("rm -rf pickle_files dat_files txt_csv_files")
        os.chdir(cwd)
    for i in dir_list:
        shutil.move(
            cwd
            + "/"
            + "gamd_simulations"
            + "/"
            + i
            + "/"
            + "system_files"
            + "/"
            + "gamd.log",
            cwd + "/" + "gamd_simulations" + "/" + i + "/" + "gamd.log",
        )
        shutil.move(
            cwd
            + "/"
            + "gamd_simulations"
            + "/"
            + i
            + "/"
            + "system_files"
            + "/"
            + "gamd-restart.dat",
            cwd
            + "/"
            + "gamd_simulations"
            + "/"
            + i
            + "/"
            + "gamd-restart.dat",
        )
        shutil.move(
            cwd
            + "/"
            + "gamd_simulations"
            + "/"
            + i
            + "/"
            + "system_files"
            + "/"
            + "md.in",
            cwd + "/" + "gamd_simulations" + "/" + i + "/" + "md.in",
        )
        shutil.move(
            cwd
            + "/"
            + "gamd_simulations"
            + "/"
            + i
            + "/"
            + "system_files"
            + "/"
            + "mdinfo",
            cwd + "/" + "gamd_simulations" + "/" + i + "/" + "mdinfo",
        )
        shutil.move(
            cwd
            + "/"
            + "gamd_simulations"
            + "/"
            + i
            + "/"
            + "system_files"
            + "/"
            + "system_final.inpcrd",
            cwd
            + "/"
            + "gamd_simulations"
            + "/"
            + i
            + "/"
            + "system_final.inpcrd",
        )
        shutil.move(
            cwd
            + "/"
            + "gamd_simulations"
            + "/"
            + i
            + "/"
            + "system_files"
            + "/"
            + "system_final.nc",
            cwd + "/" + "gamd_simulations" + "/" + i + "/" + "system_final.nc",
        )
        shutil.move(
            cwd
            + "/"
            + "gamd_simulations"
            + "/"
            + i
            + "/"
            + "system_files"
            + "/"
            + "system_final.out",
            cwd
            + "/"
            + "gamd_simulations"
            + "/"
            + i
            + "/"
            + "system_final.out",
        )
        shutil.move(
            cwd
            + "/"
            + "gamd_simulations"
            + "/"
            + i
            + "/"
            + "system_files"
            + "/"
            + "system_final.prmtop",
            cwd
            + "/"
            + "gamd_simulations"
            + "/"
            + i
            + "/"
            + "system_final.prmtop",
        )
        shutil.move(
            cwd
            + "/"
            + "gamd_simulations"
            + "/"
            + i
            + "/"
            + "system_files"
            + "/"
            + "system_final.rst",
            cwd
            + "/"
            + "gamd_simulations"
            + "/"
            + i
            + "/"
            + "system_final.rst",
        )
    for i in dir_list:
        os.chdir(target_dir + i)
        os.system("rm -rf system_files")
        os.chdir(cwd)


"""
prepare_bpti()
run_equilibration()
create_starting_structures()
add_vec_inpcrd()
add_vec_prmtop()
create_filetree()
run_simulations()
run_reweigh()
clean_for_analysis()
"""
