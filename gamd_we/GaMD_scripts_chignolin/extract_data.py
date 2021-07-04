from biopandas.pdb import PandasPdb
import pytraj as pt
import numpy as np
import pandas as pd
import os
import re
##############################################################################################################
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
# Creating rmsd.dat and rmsd_rg.dat
# Creating rmsd.dat
traj = pt.load("system_final.nc", top = "system_final.prmtop", frame_indices = frame_indices, mask='@CA')
rmsd_array = pt.distance_rmsd(traj, ref = 0, mask='@CA')
print("maximum value of RMSD is :", max(rmsd_array))
print("minimum value of RMSD is :", min(rmsd_array))
df_rmsd = pd.DataFrame (data = rmsd_array,columns=['rmsd'])
df_rmsd = df_rmsd.tail(int(datapoints))
df_rmsd.to_csv("rmsd.dat", sep = "\t",index=False, header=False)
# Creating rg.dat
rg_array = pt.radgyr(traj, '@CA')
print("maximum value of Radius of gyration is :", max(rg_array))
print("minimum value of Radius of gyration is :", min(rg_array))
df_rg = pd.DataFrame (data = rg_array,columns=['rg'])
df_rg = df_rg.tail(int(datapoints))
df_rmsd_rg = pd.concat([df_rmsd, df_rg], axis=1)
df_rmsd_rg.to_csv("rmsd_rg.dat", sep = "\t",index=False, header=False)
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
print(df_rmsd_rg.shape)
print(df_rmsd.shape)
print(df_.shape)
print(df_rmsd[df_rmsd <0.6].count())
print(df_rmsd[df_rmsd >4.0].count())
print("Folded structures form", int(df_rmsd[df_rmsd <0.6].count())/ df_rmsd.shape[0], "fraction of the total structures")
print("Unfolded structures form", int(df_rmsd[df_rmsd >4.0].count())/ df_rmsd.shape[0], "fraction of the total structures")
##############################################################################################################
