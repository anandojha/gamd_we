from biopandas.pdb import PandasPdb
import mdtraj as md
import pytraj as pt
import pandas as pd
import numpy as np
import os
import re
##############################################################################################################
index_phi = [200, 202, 204, 207]
index_psi = [584, 586, 588, 591]
jump = 5
T = 300.0
factor = 0.001987*T
##############################################################################################################
# To make data consistent with gamd.log and .nc file
with open('md.in') as f:
    lines = f.readlines()
for i in lines:
    if "nstlim =" in i:
        nstlim_line = i
    if "ntcmd =" in i:
        ntcmd_line = i
    if "ntwx =" in i :
        ntwx_line = i
x = re.findall(r'\b\d+\b',ntcmd_line )
ntcmd = int(x[0])
x = re.findall(r'\b\d+\b',nstlim_line )
nstlim = int(x[0])
x = re.findall(r'\b\d+\b',ntwx_line )
ntwx = int(x[1])
# From the .nc trajectory files, we will not consider ntcmd trajectories
leave_frames = int(ntcmd/ntwx)
no_frames = int(nstlim / ntwx)
frame_indices = []
for i in range (leave_frames,no_frames,jump):
    frame_indices.append(i)
# Recheck conditions
file = open("gamd.log", "r")
number_of_lines = 0
for line in file:
  line = line.strip("\n")
  number_of_lines += 1
file.close()
f = open('gamd.log')
fourth_line = f.readlines()[3]
if str(ntcmd) in fourth_line:
    datapoints = number_of_lines - 4
if not str(ntcmd) in fourth_line:
    datapoints = number_of_lines - 3
print(datapoints == int((nstlim - ntcmd)/ntwx))
##############################################################################################################
# Creating Psi.dat and Phi_Psi.dat
traj = pt.load("system_final.nc", top = "system_final.prmtop", frame_indices = frame_indices)
index_phi_add = "@" + str(index_phi[0]) + " @" + str(index_phi[1]) + " @" + str(index_phi[2]) + " @" + str(index_phi[3])
index_psi_add = "@" + str(index_psi[0]) + " @" + str(index_psi[1]) + " @" + str(index_psi[2]) + " @" + str(index_psi[3])
phi = pt.dihedral(traj, index_phi_add)
psi = pt.dihedral(traj, index_psi_add)
df_psi = pd.DataFrame (phi,columns=['Psi'])
df_psi = df_psi.tail(int(datapoints))
df_psi.to_csv("Psi.dat", sep = "\t",index=False, header=False)
df_phi = pd.DataFrame (psi,columns=['Phi'])
df_phi = df_phi.tail(int(datapoints))
df_phi_psi = pd.concat([df_phi, df_psi], axis=1)
df_phi_psi.to_csv("Phi_Psi.dat", sep = "\t",index=False, header=False)
##############################################################################################################
# Creating weights.dat
with open('gamd.log') as f:
    lines = f.readlines()
column_names = lines[2]
column_names = column_names.replace('#', '')
column_names = column_names.replace('\n', '')
column_names = column_names.replace(' ', '')
column_names = column_names.split(',')
list_words = ["#"]
with open('gamd.log') as oldfile, open('data.log', 'w') as newfile:
    for line in oldfile:
        if not any(word in line for word in list_words):
            newfile.write(line)
df = pd.read_csv('data.log', delim_whitespace=True, header = None)
df.columns = column_names
df["dV(kcal/mol)"] = df["Boost-Energy-Potential"] + df["Boost-Energy-Dihedral"]
df["dV(kbT)"] = df["dV(kcal/mol)"] / factor
df_ = df[["dV(kbT)","total_nstep", "dV(kcal/mol)"]]
df_ = df_[::jump]
df_.to_csv("weights.dat", sep = "\t", index=False, header=False)
os.system("rm -rf data.log")
##############################################################################################################
print(df_phi_psi.shape)
print(df_phi.shape)
print(df_.shape)
##############################################################################################################
