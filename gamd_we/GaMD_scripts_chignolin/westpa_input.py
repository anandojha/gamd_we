import numpy as np
import os
import shutil
import fnmatch
cwd = os.getcwd()
list_dir = ["1d_c1", "1d_c12", "1d_c123", "2d_c1", "2d_c12", "2d_c123"]
for i in list_dir:
    os.chdir(cwd + "/" + "we_structures" + "/" + i)
    files = os.listdir('.')
    file_to_find = "*.pdb"
    pdb_list = []
    for x in files:
        if fnmatch.fnmatch(x, file_to_find):
            pdb_list.append(x)
    for j in pdb_list:
        inpcrd_file = j[:-4] + ".inpcrd"
        filename = "input_" + j[:-4]  + ".leap"
        file = open(filename,"w")
        file.write("source leaprc.protein.ff14SB"                                  + "\n")
        file.write("source leaprc.water.tip3p"                                     + "\n")
        file.write("set default FlexibleWater on"                                  + "\n")
        file.write("set default PBRadii mbondi2"                                   + "\n")
        file.write("pdb = loadpdb " + j                                            + "\n")
        file.write("saveamberparm pdb " + j[:-4] + ".prmtop " + j[:-4] + ".inpcrd" + "\n")
        file.write("quit"                                                          + "\n")
        file.close() 
    files = os.listdir('.')
    file_to_find = "*.leap"
    leap_list = []
    for y in files:
        if fnmatch.fnmatch(y, file_to_find):
               leap_list.append(y)
    for k in leap_list:
        command = "tleap -f {}".format(k)
        os.system(command)
    os.system("rm -rf leap.log")
    os.system("rm -rf *prmtop*")
    os.system("rm -rf *leap*") 
    os.system("rm -rf bstates")
    os.system("mkdir bstates")
    for j in pdb_list:
        shutil.move(cwd + "/" + "we_structures" + "/" + i + "/" + j[:-4] + ".inpcrd",cwd + "/" + "we_structures" + "/" + i + "/" + "bstates" + "/" + j[:-4] + ".inpcrd")
    os.chdir(cwd)
    
os.system("rm -rf westpa_inputs")
os.system("mkdir westpa_inputs")
for l in list_dir:
    os.chdir(cwd + "/" + "westpa_inputs")
    command = "rm -rf {}".format(l)
    os.system(command)
    command = "mkdir {}".format(l)
    os.system(command)
    shutil.move(cwd + "/" + "we_structures" + "/" + l + "/" + "bstates", cwd + "/" + "westpa_inputs" + "/" + l + "/" + "bstates")
    os.chdir(cwd)   
shutil.copy(cwd + "/" + "we_inputs" + "/" + "we_input_c1_1d.txt", cwd + "/" + "westpa_inputs" + "/" + list_dir[0]+ "/" + "we_input_c1_1d.txt")
shutil.copy(cwd + "/" + "we_inputs" + "/" + "we_input_c12_1d.txt", cwd + "/" + "westpa_inputs" + "/" + list_dir[1]+ "/" + "we_input_c12_1d.txt")
shutil.copy(cwd + "/" + "we_inputs" + "/" + "we_input_c123_1d.txt", cwd + "/" + "westpa_inputs" + "/" + list_dir[2]+ "/" + "we_input_c123_1d.txt")
shutil.copy(cwd + "/" + "we_inputs" + "/" + "we_input_c1_2d.txt", cwd + "/" + "westpa_inputs" + "/" + list_dir[3]+ "/" + "we_input_c1_2d.txt")
shutil.copy(cwd + "/" + "we_inputs" + "/" + "we_input_c12_2d.txt", cwd + "/" + "westpa_inputs" + "/" + list_dir[4]+ "/" + "we_input_c12_2d.txt")
shutil.copy(cwd + "/" + "we_inputs" + "/" + "we_input_c123_2d.txt", cwd + "/" + "westpa_inputs" + "/" + list_dir[5]+ "/" + "we_input_c123_2d.txt")

for i in list_dir:
    os.chdir(cwd + "/" + "westpa_inputs" + "/" + i)
    for file in os.listdir('.'):
        if fnmatch.fnmatch(file, "*.txt"):
            file_to_rename = file
            f = open(file_to_rename, "rt")
            data = f.read()
            data = data.replace('pdb', 'inpcrd')
            f.close()
            f = open(file_to_rename, "wt")
            f.write(data)
            f.close()
            os.rename(file_to_rename, "BASIS_STATES")
            os.chdir(cwd)
            
for i in list_dir:
    os.chdir(cwd + "/" + "westpa_inputs" + "/" + i)
    os.mkdir ("CONFIG")
    shutil.copy(cwd + "/" + "system_files" + "/" + "system_final.prmtop", cwd + "/" + "westpa_inputs" + "/" + i + "/" + "CONFIG" + "/" + "system_final.prmtop")
    os.chdir(cwd)  
