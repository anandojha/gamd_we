import os
import shutil
import pandas as pd
import pickle as pk
import pytraj as pt
cwd = os.getcwd()
os.system("rm -rf we_structures")
os.system("mkdir we_structures")
os.chdir(cwd + "/" + "we_structures")
os.system("mkdir 1d_c1")
os.system("mkdir 1d_c12")
os.system("mkdir 1d_c123")
os.system("mkdir 2d_c1")
os.system("mkdir 2d_c12")
os.system("mkdir 2d_c123")
os.chdir(cwd)
df1 = pd.read_csv("df_1d.csv")
index = df1['index'].tolist()
frame = df1['frame_index'].tolist()
index_frame = dict(zip(frame, index))
df2 = pd.read_csv("ref_1d.txt", sep=' ',delimiter=None, header='infer')
index_ = df2['index'].tolist()
bins = df2['bins'].tolist()
index_bins = dict(zip(index_, bins)) 
#### 1d 
with open("frames_c1_1d.pickle", "rb") as input_file:
    frames_c1_1d = pk.load(input_file)
for i in frames_c1_1d:
    j = index_frame[i]
    frame_index = frames_c1_1d.index(i)
    k = index_bins[j]
    k = k.strip("[]")
    k = k.replace(" , ", "_")    
    traj = pt.load("system_final.nc", top = "system_final.prmtop", frame_indices = [i])
    frame_pdb = str(frame_index) + "_" + k + "_1d_c1_" + str(i) + ".pdb"
    pt.save(frame_pdb, traj, overwrite = True)
    target_dir  = cwd + "/" + "we_structures" + "/" + "1d_c1"
    shutil.move(cwd +  "/" + frame_pdb , target_dir + "/" + frame_pdb)    
with open("frames_c12_1d.pickle", "rb") as input_file:
    frames_c12_1d = pk.load(input_file)
for i in frames_c12_1d:
    j = index_frame[i]
    frame_index = frames_c12_1d.index(i)
    k = index_bins[j]
    k = k.strip("[]")
    k = k.replace(" , ", "_")
    traj = pt.load("system_final.nc", top = "system_final.prmtop", frame_indices = [i])
    frame_pdb = str(frame_index) + "_" + k + "_1d_c12_" + str(i) + ".pdb"
    pt.save(frame_pdb, traj, overwrite = True)
    target_dir  = cwd + "/" + "we_structures" + "/" + "1d_c12"
    shutil.move(cwd +  "/" + frame_pdb , target_dir + "/" + frame_pdb)   
with open("frames_c123_1d.pickle", "rb") as input_file:
    frames_c123_1d = pk.load(input_file)
for i in frames_c123_1d:
    j = index_frame[i]
    frame_index = frames_c123_1d.index(i)
    k = index_bins[j]
    k = k.strip("[]")
    k = k.replace(" , ", "_")
    traj = pt.load("system_final.nc", top = "system_final.prmtop", frame_indices = [i])
    frame_pdb = str(frame_index) + "_" + k + "_1d_c123_" + str(i) + ".pdb"
    pt.save(frame_pdb, traj, overwrite = True)
    target_dir  = cwd + "/" + "we_structures" + "/" + "1d_c123"
    shutil.move(cwd +  "/" + frame_pdb , target_dir + "/" + frame_pdb)   
df1 = pd.read_csv("df_2d.csv")
index = df1['index'].tolist()
frame = df1['frame_index'].tolist()
index_frame = dict(zip(frame, index))
df2 = pd.read_csv("ref_2d.txt", sep=' ',delimiter=None, header='infer')
index_ = df2['index'].tolist()
bins = df2['XY'].tolist()
index_bins = dict(zip(index_, bins)) 
#### 2d 
with open("frames_c1_2d.pickle", "rb") as input_file:
    frames_c1_2d = pk.load(input_file)
for i in frames_c1_2d:
    j = index_frame[i]
    frame_index = frames_c1_2d.index(i)
    k = index_bins[j]
    k = k.strip("[]")
    k = k.replace("] , [", "_")
    k = k.replace(", ", "_")
    traj = pt.load("system_final.nc", top = "system_final.prmtop", frame_indices = [i])
    frame_pdb = str(frame_index) + "_" + k + "_2d_c1_" + str(i) + ".pdb"
    pt.save(frame_pdb, traj, overwrite = True)
    target_dir  = cwd + "/" + "we_structures" + "/" + "2d_c1"
    shutil.move(cwd +  "/" + frame_pdb , target_dir + "/" + frame_pdb)    
with open("frames_c12_2d.pickle", "rb") as input_file:
    frames_c12_2d = pk.load(input_file)
for i in frames_c12_2d:
    j = index_frame[i]
    frame_index = frames_c12_2d.index(i)
    k = index_bins[j]
    k = k.strip("[]")
    k = k.replace("] , [", "_")
    k = k.replace(", ", "_")
    traj = pt.load("system_final.nc", top = "system_final.prmtop", frame_indices = [i])
    frame_pdb = str(frame_index) + "_" + k + "_2d_c12_" + str(i) + ".pdb"
    pt.save(frame_pdb, traj, overwrite = True)
    target_dir  = cwd + "/" + "we_structures" + "/" + "2d_c12"
    shutil.move(cwd +  "/" + frame_pdb , target_dir + "/" + frame_pdb)    
with open("frames_c123_2d.pickle", "rb") as input_file:
    frames_c123_2d = pk.load(input_file)
for i in frames_c123_2d:
    j = index_frame[i]
    frame_index = frames_c123_2d.index(i)
    k = index_bins[j]
    k = k.strip("[]")
    k = k.replace("] , [", "_")
    k = k.replace(", ", "_")
    traj = pt.load("system_final.nc", top = "system_final.prmtop", frame_indices = [i])
    frame_pdb = str(frame_index) + "_" + k + "_2d_c123_" + str(i) + ".pdb"
    pt.save(frame_pdb, traj, overwrite = True)
    target_dir  = cwd + "/" + "we_structures" + "/" + "2d_c123"
    shutil.move(cwd +  "/" + frame_pdb , target_dir + "/" + frame_pdb)
