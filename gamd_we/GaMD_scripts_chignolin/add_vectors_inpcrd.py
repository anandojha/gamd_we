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
vectors = (nvt_simulation_box_vectors[0][0])*10, (nvt_simulation_box_vectors[1][1])*10, (nvt_simulation_box_vectors[2][2])*10
vectors = (round(vectors[0], 7), round(vectors[1], 7), round(vectors[2], 7))
last_line = "  " + str(vectors[0]) + "  " + str(vectors[1]) + "  " + str(vectors[2]) + "  90.0000000" + "  90.0000000"  + "  90.0000000" 
with open("system_final.inpcrd", "a+") as f:
    f.write(last_line)
os.system("rm -rf nvt_simulation_box_vectors.pkl")
os.chdir(cwd)  
