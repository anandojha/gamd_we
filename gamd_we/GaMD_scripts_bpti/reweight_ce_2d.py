import os
import random
import itertools
import statistics
import numpy as np
import pickle as pk
import pandas as pd
from math import exp
##############################################################################################################################################################################################
binspace = 10
n_structures = 4
Xdim = [-180,180]
Ydim = [-180,180]
T = 300.0            
beta = 1.0/(0.001987*float(T))
min_prob = 0.000001
##############################################################################################################################################################################################
def create_bins(lower_bound, width, upper_bound):
    bins = []
    for low in range(lower_bound, upper_bound, width):
        bins.append([low, low+width])
    return bins
def find_bin(value, bins):
    for i in range(0, len(bins)):
        if bins[i][0] <= value < bins[i][1]:
            return i
    return -1
##############################################################################################################################################################################################
df_Phi_Psi = pd.read_csv('Phi_Psi.dat', delim_whitespace=True, header = None)
df_Phi_Psi.columns = ["Phi","Psi"]
df_weight = pd.read_csv('weights.dat', delim_whitespace=True, header = None)
df_weight.columns = ["dV_kBT", "timestep", "dVkcalmol"]
##############################################################################################################################################################################################
sum_total = df_Phi_Psi.shape[0]
binsX = np.arange(float(Xdim[0]),(float(Xdim[1]) + binspace),binspace)
binsY = np.arange(float(Ydim[0]),(float(Ydim[1]) + binspace),binspace)
hist2D, hist_edgesX, hist_edgesY = np.histogram2d(df_Phi_Psi["Phi"].values.tolist(),df_Phi_Psi["Psi"].values.tolist(),bins = (binsX, binsY),weights=None)
pstarA_2D = [i / sum_total for i in list(hist2D)]
bins_tuple_X = create_bins(lower_bound=int(Xdim[0]),width=binspace,upper_bound=int(Xdim[1]))
bins_tuple_Y = create_bins(lower_bound=int(Ydim[0]),width=binspace,upper_bound=int(Ydim[1]))
bins = []
for i in range(len(bins_tuple_X)):
    for j in range(len(bins_tuple_Y)):
        bins.append([bins_tuple_X[i],bins_tuple_Y[j]])
pstarA = [item for elem in pstarA_2D for item in elem]
hist = [item for elem in hist2D for item in elem]
hist = [int(i) for i in hist]
#############################################################################################################################################################################################
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
    binned_weights.append((binned_weights_2D[i][0] * len(bins_tuple_Y)) + (binned_weights_2D[i][1] + 1))
df_index = pd.DataFrame(binned_weights)
df_index.columns = ["index"]
df_index["index"] = df_index["index"]  - 1
#############################################################################################################################################################################################
df = pd.concat([df_index, df_Phi_Psi, df_weight], axis=1)
dV_c1 = []
dV_c2 = []
dV_c3 = []
dV = []
for i in range(len(bins)):
    df_i = df.loc[(df['index'] == i)]
    dV_list = df_i["dVkcalmol"].values.tolist()
    if len(dV_list)>=10:
        dV_c1.append(statistics.mean(dV_list))
        dV_c2.append(statistics.mean([i**2 for i in dV_list]) - (statistics.mean(dV_list))**2)
        dV_c3.append(statistics.mean([i**3 for i in dV_list]) - 3 * (statistics.mean([i**2 for i in dV_list])) * (statistics.mean(dV_list)) + 2 * (statistics.mean(dV_list))**3)
    if len(dV_list)<10:
        dV_c1.append(0)   
        dV_c2.append(0)
        dV_c3.append(0)
    dV.append(dV_list)
#############################################################################################################################################################################################
c1 = [i * beta for i in dV_c1]
c2 = [i * ((beta**2)/2) for i in dV_c2]
c3 = [i * ((beta**3)/6) for i in dV_c3]
c1 = c1
c12 = [a + b for a, b in zip(c1, c2)]
c123 = [a + b for a, b in zip(c12, c3)]
for i in range(len(c1)):
    if c1[i]>=700:
        c1[i] = 700
for i in range(len(c12)):
    if c12[i]>=700:
        c12[i] = 700
for i in range(len(c123)):
    if c123[i]>=700:
        c123[i] = 700  
ensemble_average_c1 = [exp(i) for i in c1]
ensemble_average_c12 = [exp(i) for i in c12]
ensemble_average_c123 = [exp(i) for i in c123]
numerator_c1 = [a * b for a, b in zip(pstarA, ensemble_average_c1)]
numerator_c12 = [a * b for a, b in zip(pstarA, ensemble_average_c12)]
numerator_c123 = [a * b for a, b in zip(pstarA, ensemble_average_c123)]
#############################################################################################################################################################################################
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
#############################################################################################################################################################################################
data_c1 = list(zip(bins,pA_c1))
data_c12 = list(zip(bins,pA_c12))
data_c123 = list(zip(bins,pA_c123))
#############################################################################################################################################################################################
df_c1 = pd.DataFrame(data_c1, columns=['bins','pA_c1'])
df_c12 = pd.DataFrame(data_c12, columns=['bins','pA_c12'])
df_c123 = pd.DataFrame(data_c123, columns=['bins','pA_c123'])
#############################################################################################################################################################################################
df_c1.to_csv('c1_2d.txt', header=True, index=None, sep=" " , mode='w')
with open('c1_2d.txt', 'r') as f1, open('pA_c1_2d.txt', 'w') as f2:
    for line in f1:
        f2.write(line.replace('"', '').replace("'", ""))
os.system("rm -rf c1_2d.txt")
####c12
df_c12.to_csv('c12_2d.txt', header=True, index=None, sep=" " , mode='w')
with open('c12_2d.txt', 'r') as f1, open('pA_c12_2d.txt', 'w') as f2:
    for line in f1:
        f2.write(line.replace('"', '').replace("'", ""))
os.system("rm -rf c12_2d.txt")
####c123
df_c123.to_csv('c123_2d.txt', header=True, index=None, sep=" " , mode='w')
with open('c123_2d.txt', 'r') as f1, open('pA_c123_2d.txt', 'w') as f2:
    for line in f1:
        f2.write(line.replace('"', '').replace("'", ""))
os.system("rm -rf c123_2d.txt")
##############################################################################################################################################################################################
####c1_arranged
df_c1_arranged = df_c1.sort_values(by='pA_c1', ascending=False)
df_c1_arranged = df_c1_arranged[df_c1_arranged.pA_c1 > min_prob]
df_c1_arranged.to_csv('c1_arranged_2d.txt', header=True, index=None, sep=" " , mode='w')
with open('c1_arranged_2d.txt', 'r') as f1, open('pA_c1_arranged_2d.txt', 'w') as f2:
    for line in f1:
        f2.write(line.replace('"', '').replace("'", ""))
os.system("rm -rf c1_arranged_2d.txt")
####c12_arranged
df_c12_arranged = df_c12.sort_values(by='pA_c12', ascending=False)
df_c12_arranged = df_c12_arranged[df_c12_arranged.pA_c12 > min_prob]
df_c12_arranged.to_csv('c12_arranged_2d.txt', header=True, index=None, sep=" " , mode='w')
with open('c12_arranged_2d.txt', 'r') as f1, open('pA_c12_arranged_2d.txt', 'w') as f2:
    for line in f1:
        f2.write(line.replace('"', '').replace("'", ""))
os.system("rm -rf c12_arranged_2d.txt")
####c123_arranged
df_c123_arranged = df_c123.sort_values(by='pA_c123', ascending=False)
df_c123_arranged = df_c123_arranged[df_c123_arranged.pA_c123 > min_prob]
df_c123_arranged.to_csv('c123_arranged_2d.txt', header=True, index=None, sep=" " , mode='w')
with open('c123_arranged_2d.txt', 'r') as f1, open('pA_c123_arranged_2d.txt', 'w') as f2:
    for line in f1:
        f2.write(line.replace('"', '').replace("'", ""))
os.system("rm -rf c123_arranged_2d.txt")
##############################################################################################################################################################################################
####c1_arranged
df_c1_arranged['index'] = df_c1_arranged.index
index_list_c1 = df_c1_arranged['index'].tolist()
df['frame_index'] = df.index
df_frame_index = df[["frame_index", "index"]]
frame_indices_c1 = []
index_indces_c1 = []
for i in index_list_c1:
    df_index_list_c1 = df_frame_index.loc[df_frame_index['index'] == i]
    frame_c1 = df_index_list_c1["frame_index"].tolist()
    frame_indices_c1.append(frame_c1)
    index_c1 = [i] * len(frame_c1)
    index_indces_c1.append(index_c1)
frame_indices_c1 = [item for elem in frame_indices_c1 for item in elem]
index_indces_c1 = [item for elem in index_indces_c1 for item in elem]
df_c1_frame = pd.DataFrame (frame_indices_c1,columns=["frame_index"])
df_c1_index = pd.DataFrame (index_indces_c1,columns=["index"])
df_c1_frame_index = pd.concat([df_c1_frame, df_c1_index], axis=1)
df_c1_frame_index = df_c1_frame_index.groupby('index').filter(lambda x : len(x)>=10)
df_c1_frame_index.to_csv('c1_frame_index_2d.txt', header=True, index=None, sep=" " , mode='w')
with open('c1_frame_index_2d.txt', 'r') as f1, open('c1_frame_2d.txt', 'w') as f2:
    for line in f1:
        f2.write(line.replace('"', '').replace("'", ""))
os.system("rm -rf c1_frame_index_2d.txt")
####c12_arranged
df_c12_arranged['index'] = df_c12_arranged.index
index_list_c12 = df_c12_arranged['index'].tolist()
df['frame_index'] = df.index
df_frame_index = df[["frame_index", "index"]]
frame_indices_c12 = []
index_indces_c12 = []
for i in index_list_c12:
    df_index_list_c12 = df_frame_index.loc[df_frame_index['index'] == i]
    frame_c12 = df_index_list_c12["frame_index"].tolist()
    frame_indices_c12.append(frame_c12)
    index_c12 = [i] * len(frame_c12)
    index_indces_c12.append(index_c12)
frame_indices_c12 = [item for elem in frame_indices_c12 for item in elem]
index_indces_c12 = [item for elem in index_indces_c12 for item in elem]
df_c12_frame = pd.DataFrame (frame_indices_c12,columns=["frame_index"])
df_c12_index = pd.DataFrame (index_indces_c12,columns=["index"])
df_c12_frame_index = pd.concat([df_c12_frame, df_c12_index], axis=1)
df_c12_frame_index = df_c12_frame_index.groupby('index').filter(lambda x : len(x)>=10)
df_c12_frame_index.to_csv('c12_frame_index_2d.txt', header=True, index=None, sep=" " , mode='w')
with open('c12_frame_index_2d.txt', 'r') as f1, open('c12_frame_2d.txt', 'w') as f2:
    for line in f1:
        f2.write(line.replace('"', '').replace("'", ""))
os.system("rm -rf c12_frame_index_2d.txt")
####c123_arranged
df_c123_arranged['index'] = df_c123_arranged.index
df_c123_arranged['index'] = df_c123_arranged.index
index_list_c123 = df_c123_arranged['index'].tolist()
df['frame_index'] = df.index
df_frame_index = df[["frame_index", "index"]]
frame_indices_c123 = []
index_indces_c123 = []
for i in index_list_c123:
    df_index_list_c123 = df_frame_index.loc[df_frame_index['index'] == i]
    frame_c123 = df_index_list_c123["frame_index"].tolist()
    frame_indices_c123.append(frame_c123)
    index_c123 = [i] * len(frame_c123)
    index_indces_c123.append(index_c123)
frame_indices_c123 = [item for elem in frame_indices_c123 for item in elem]
index_indces_c123 = [item for elem in index_indces_c123 for item in elem]
df_c123_frame = pd.DataFrame (frame_indices_c123,columns=["frame_index"])
df_c123_index = pd.DataFrame (index_indces_c123,columns=["index"])
df_c123_frame_index = pd.concat([df_c123_frame, df_c123_index], axis=1)
df_c123_frame_index = df_c123_frame_index.groupby('index').filter(lambda x : len(x)>=10)
df_c123_frame_index.to_csv('c123_frame_index_2d.txt', header=True, index=None, sep=" " , mode='w')
with open('c123_frame_index_2d.txt', 'r') as f1, open('c123_frame_2d.txt', 'w') as f2:
    for line in f1:
        f2.write(line.replace('"', '').replace("'", ""))
os.system("rm -rf c123_frame_index_2d.txt")
##############################################################################################################################################################################################
####c1
indices_c1_2d = df_c1_frame_index["index"].unique()
frames_c1 = []
for i in indices_c1_2d:
    x = df_c1_frame_index.loc[df_c1_frame_index['index'] == i]
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
    x = df_c12_frame_index.loc[df_c12_frame_index['index'] == i]
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
    x = df_c123_frame_index.loc[df_c123_frame_index['index'] == i]
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
prob_c1_2d_list = list(itertools.chain.from_iterable(itertools.repeat(x, n_structures) for x in prob_c1_2d_list))
prob_c1_2d_list = [x / n_structures for x in prob_c1_2d_list]
with open("prob_c1_2d_list.pickle", "wb") as f:
    pk.dump(prob_c1_2d_list, f)
####c12
prob_c12_2d_list = []
for i in indices_c12_2d:
    prob_c12_2d_list.append(df_c12["pA_c12"][i])
prob_c12_2d_list = list(itertools.chain.from_iterable(itertools.repeat(x, n_structures) for x in prob_c12_2d_list))
prob_c12_2d_list = [x / n_structures for x in prob_c12_2d_list]
with open("prob_c12_2d_list.pickle", "wb") as f:
    pk.dump(prob_c12_2d_list, f)
####c123
prob_c123_2d_list = []
for i in indices_c123_2d:
    prob_c123_2d_list.append(df_c123["pA_c123"][i])
prob_c123_2d_list = list(itertools.chain.from_iterable(itertools.repeat(x, n_structures) for x in prob_c123_2d_list))
prob_c123_2d_list = [x / n_structures for x in prob_c123_2d_list]
with open("prob_c123_2d_list.pickle", "wb") as f:
    pk.dump(prob_c123_2d_list, f)
##############################################################################################################################################################################################
ref_df_2d = pd.DataFrame(bins, columns = ["binsX", "binsY"])
ref_df_2d['XY'] = ref_df_2d.agg(lambda x: f"{x['binsX']} , {x['binsX']}", axis=1)
ref_df_2d = ref_df_2d[['XY']]
index_ref_2d = []
for i in range(len(bins_tuple_X) * len(bins_tuple_Y)):
    index_ref_2d.append(i)
index_ref_df_2d = pd.DataFrame(index_ref_2d, columns = ["index"])
df_ref_2d = pd.concat([ref_df_2d, index_ref_df_2d], axis=1)
df_ref_2d.to_csv('ref_2d.txt', header=True, index=None, sep=" " , mode='w')    
##############################################################################################################################################################################################
df.to_csv('df_2d.csv', index=False)  
os.system("rm -rf __pycache__")
print("Successfully Completed Reweighing")
