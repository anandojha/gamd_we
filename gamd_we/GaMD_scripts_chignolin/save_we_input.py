import os
import re
import shutil
import pickle as pk
import pandas as pd
cwd = os.getcwd()
target_dir = cwd + "/" + "we_structures" 
dir_list = ["1d_c1", "1d_c12", "1d_c123", "2d_c1", "2d_c12", "2d_c123"]
for i in dir_list:
    os.chdir(target_dir + "/" + i)
    pdbs = os.listdir('.') 
    pickle_file = "pdb_" + i + ".pickle"
    with open(pickle_file, "wb") as f:
        pk.dump(pdbs, f)
    shutil.move(target_dir + "/" + i + "/" + pickle_file, cwd + "/" + pickle_file)  
    os.chdir(cwd)
#c1_1d
with open("prob_c1_1d_list.pickle", "rb") as input_file:
    prob_c1_1d_list = pk.load(input_file)
prob_c1_1d_list = [i / min(prob_c1_1d_list) for i in prob_c1_1d_list]
prob_c1_1d_list = [i / sum(prob_c1_1d_list) for i in prob_c1_1d_list]  
with open("pdb_1d_c1.pickle", "rb") as input_file:
    pdb_1d_c1 = pk.load(input_file)
pdb_1d_c1_index = []
for i in range(len(pdb_1d_c1)):
    pdb_1d_c1_index.append(int(re.findall(r'\d+', pdb_1d_c1[i])[0]))
df = pd.DataFrame(list(zip(pdb_1d_c1, pdb_1d_c1_index)),columns=['pdb_name', 'pdb_index'])
df = df.sort_values(by=['pdb_index'])
index_list = df["pdb_index"].values.tolist()
pdb_list = df["pdb_name"].values.tolist()
df_merged = pd.DataFrame(list(zip(index_list, prob_c1_1d_list, pdb_list)),columns=['index','probability', 'pdb_name'])
df_merged.to_csv('we_input_c1_1d.txt', header=False, index=None, sep=" " , mode='w') 
#c12_1d
with open("prob_c12_1d_list.pickle", "rb") as input_file:
    prob_c12_1d_list = pk.load(input_file)
prob_c12_1d_list = [i / min(prob_c12_1d_list) for i in prob_c12_1d_list]
prob_c12_1d_list = [i / sum(prob_c12_1d_list) for i in prob_c12_1d_list]      
with open("pdb_1d_c12.pickle", "rb") as input_file:
    pdb_1d_c12 = pk.load(input_file)
pdb_1d_c12_index = []
for i in range(len(pdb_1d_c12)):
    pdb_1d_c12_index.append(int(re.findall(r'\d+', pdb_1d_c12[i])[0]))
df = pd.DataFrame(list(zip(pdb_1d_c12, pdb_1d_c12_index)),columns=['pdb_name', 'pdb_index'])
df = df.sort_values(by=['pdb_index'])
index_list = df["pdb_index"].values.tolist()
pdb_list = df["pdb_name"].values.tolist()
df_merged = pd.DataFrame(list(zip(index_list, prob_c12_1d_list, pdb_list)),columns=['index','probability', 'pdb_name'])
df_merged.to_csv('we_input_c12_1d.txt', header=False, index=None, sep=" " , mode='w') 
#c123_1d
with open("prob_c123_1d_list.pickle", "rb") as input_file:
    prob_c123_1d_list = pk.load(input_file)
prob_c123_1d_list = [i / min(prob_c123_1d_list) for i in prob_c123_1d_list]
prob_c123_1d_list = [i / sum(prob_c123_1d_list) for i in prob_c123_1d_list]      
with open("pdb_1d_c123.pickle", "rb") as input_file:
    pdb_1d_c123 = pk.load(input_file)
pdb_1d_c123_index = []
for i in range(len(pdb_1d_c123)):
    pdb_1d_c123_index.append(int(re.findall(r'\d+', pdb_1d_c123[i])[0]))   
df = pd.DataFrame(list(zip(pdb_1d_c123, pdb_1d_c123_index)),columns=['pdb_name', 'pdb_index'])
df = df.sort_values(by=['pdb_index'])
index_list = df["pdb_index"].values.tolist()
pdb_list = df["pdb_name"].values.tolist()
df_merged = pd.DataFrame(list(zip(index_list, prob_c123_1d_list, pdb_list)),columns=['index','probability', 'pdb_name'])
df_merged.to_csv('we_input_c123_1d.txt', header=False, index=None, sep=" " , mode='w')     
#c1_2d
with open("prob_c1_2d_list.pickle", "rb") as input_file:
    prob_c1_2d_list = pk.load(input_file)
prob_c1_2d_list = [i / min(prob_c1_2d_list) for i in prob_c1_2d_list]
prob_c1_2d_list = [i / sum(prob_c1_2d_list) for i in prob_c1_2d_list]  
with open("pdb_2d_c1.pickle", "rb") as input_file:
    pdb_2d_c1 = pk.load(input_file)
pdb_2d_c1_index = []
for i in range(len(pdb_2d_c1)):
    pdb_2d_c1_index.append(int(re.findall(r'\d+', pdb_2d_c1[i])[0]))
df = pd.DataFrame(list(zip(pdb_2d_c1, pdb_2d_c1_index)),columns=['pdb_name', 'pdb_index'])
df = df.sort_values(by=['pdb_index'])
index_list = df["pdb_index"].values.tolist()
pdb_list = df["pdb_name"].values.tolist()
df_merged = pd.DataFrame(list(zip(index_list, prob_c1_2d_list, pdb_list)),columns=['index','probability', 'pdb_name'])
df_merged.to_csv('we_input_c1_2d.txt', header=False, index=None, sep=" " , mode='w')     
#c12_2d
with open("prob_c12_2d_list.pickle", "rb") as input_file:
    prob_c12_2d_list = pk.load(input_file)
prob_c12_2d_list = [i / min(prob_c12_2d_list) for i in prob_c12_2d_list]
prob_c12_2d_list = [i / sum(prob_c12_2d_list) for i in prob_c12_2d_list]      
with open("pdb_2d_c12.pickle", "rb") as input_file:
    pdb_2d_c12 = pk.load(input_file)
pdb_2d_c12_index = []
for i in range(len(pdb_2d_c12)):
    pdb_2d_c12_index.append(int(re.findall(r'\d+', pdb_2d_c12[i])[0]))
df = pd.DataFrame(list(zip(pdb_2d_c12, pdb_2d_c12_index)),columns=['pdb_name', 'pdb_index'])
df = df.sort_values(by=['pdb_index'])
index_list = df["pdb_index"].values.tolist()
pdb_list = df["pdb_name"].values.tolist()
df_merged = pd.DataFrame(list(zip(index_list, prob_c12_2d_list, pdb_list)),columns=['index','probability', 'pdb_name'])
df_merged.to_csv('we_input_c12_2d.txt', header=False, index=None, sep=" " , mode='w')         
#c123_2d
with open("prob_c123_2d_list.pickle", "rb") as input_file:
    prob_c123_2d_list = pk.load(input_file)
prob_c123_2d_list = [i / min(prob_c123_2d_list) for i in prob_c123_2d_list]
prob_c123_2d_list = [i / sum(prob_c123_2d_list) for i in prob_c123_2d_list]      
with open("pdb_2d_c123.pickle", "rb") as input_file:
    pdb_2d_c123 = pk.load(input_file)
pdb_2d_c123_index = []
for i in range(len(pdb_2d_c123)):
    pdb_2d_c123_index.append(int(re.findall(r'\d+', pdb_2d_c123[i])[0]))
df = pd.DataFrame(list(zip(pdb_2d_c123, pdb_2d_c123_index)),columns=['pdb_name', 'pdb_index'])
df = df.sort_values(by=['pdb_index'])
index_list = df["pdb_index"].values.tolist()
pdb_list = df["pdb_name"].values.tolist()
df_merged = pd.DataFrame(list(zip(index_list, prob_c123_2d_list, pdb_list)),columns=['index','probability', 'pdb_name'])
df_merged.to_csv('we_input_c123_2d.txt', header=False, index=None, sep=" " , mode='w')   
