import pickle as pk
import fileinput
import shutil
import re
import os
cwd = os.getcwd()
target_dir = cwd + "/" + "starting_structures"
shutil.copy(cwd  + "/" + "equilibration" + "/" + "nvt_simulation_box_vectors.pkl" ,  target_dir + "/" + "nvt_simulation_box_vectors.pkl")
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
os.chdir(target_dir)
with open('nvt_simulation_box_vectors.pkl', 'rb') as f:
    nvt_simulation_box_vectors = pk.load(f)
nvt_simulation_box_vectors = create_vectors(nvt_simulation_box_vectors)
vectors = nvt_simulation_box_vectors[0][0], nvt_simulation_box_vectors[1][1], nvt_simulation_box_vectors[2][2]
vectors = round(vectors[0], 7), round(vectors[1], 7), round(vectors[2], 7)
oldbeta = "9.00000000E+01"
x = str(vectors[0]) + str(0) + "E+" + "01"
y = str(vectors[1]) + str(0) + "E+" + "01"
z = str(vectors[2]) + str(0) + "E+" + "01"
line1 = "%FLAG BOX_DIMENSIONS"
line2 = "%FORMAT(5E16.8)"
line3 = "  " + oldbeta + "  " + x + "  " + y + "  " + z
with open('system_final.prmtop') as i, open('system_intermediate_final.prmtop', 'w') as f:
    for line in i:
       if line.startswith('%FLAG RADIUS_SET'):
           line = line1 + "\n" + line2 + "\n" + line3 + "\n" + line
       f.write(line)
os.system("rm -rf system_final.prmtop")
os.system("mv system_intermediate_final.prmtop system_final.prmtop")
os.system("rm -rf nvt_simulation_box_vectors.pkl")
os.chdir(cwd)
