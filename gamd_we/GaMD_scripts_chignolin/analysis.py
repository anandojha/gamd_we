import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.image as mpimg
import matplotlib.style as style 

cwd = os.getcwd()
os.chdir(cwd + "/" + "westpa_dir")
dir_list = ["dihedral_threshold_lower", "dihedral_threshold_upper", "dual_threshold_lower", "dual_threshold_upper", "total_threshold_lower", "total_threshold_upper"]
we_list = ["1d_c1", "1d_c12", "1d_c123", "2d_c1", "2d_c12", "2d_c123"]
confs = []
for i in dir_list:
    conf_within  = []
    for j in we_list:
        os.chdir(cwd + "/" + "westpa_dir" + "/" + str(i))
        os.chdir(cwd + "/" + "westpa_dir" + "/" + str(i) + "/" + str(j))
        if len(open("BASIS_STATES").readlines(  )) > 0 :
            count1 = len(open("BASIS_STATES").readlines(  ))
            count2 = len(open("BASIS_STATES_CORRECTED_RST").readlines(  ))
            conf = str(i),str(j),count1,count2
            conf_within.append(conf)
    confs.append(conf_within)
print(confs)
os.chdir(cwd)

corrected_list = []
for i in range(len(confs)):
    corrected_list_1 = []
    for j in range(len(confs[i])):
        corrected_list_1.append(confs[i][j][3])
    corrected_list.append(corrected_list_1)
print(corrected_list)

expanse_list = []
for i in range(len(confs)):
    expanse_list_1 = []
    for j in range(len(confs[i])):
        expanse_list_1.append(confs[i][j][1])
    expanse_list.append(expanse_list_1)
print(expanse_list)

x0 = expanse_list[0]
y0 = corrected_list[0]
x1 = expanse_list[1]
y1 = corrected_list[1]
x2 = expanse_list[2]
y2 = corrected_list[2]
x3 = expanse_list[3]
y3 = corrected_list[3]
x4 = expanse_list[4]
y4 = corrected_list[4]
x5 = expanse_list[5]
y5 = corrected_list[5]

y = y0
x = x0
title = "Configurations vs Different Expansions" + " for " + dir_list[0]
print(title)
sns.set(font_scale=1)
plt.rcParams['figure.figsize'] = (8,4)
plt.rcParams['font.family'] = "serif"
style.use('fivethirtyeight')
g = sns.barplot(y,x,palette=("binary"))
g.grid(False)
g.set_title(title)
g.set(xlabel='Configurations', ylabel='Expansion')
ax = g
for i, v in enumerate(y):
    ax.text(v + 1, i + .25, str(v), color='black', fontweight = "bold")
fig_name = dir_list[0]
plt.savefig(fig_name, bbox_inches='tight')
plt.show(block=False)
plt.pause(1)
plt.close()

y = y1
x = x1
title = "Configurations vs Different Expansions" + " for " + dir_list[1]
print(title)
sns.set(font_scale=1)
plt.rcParams['figure.figsize'] = (8,4)
plt.rcParams['font.family'] = "serif"
style.use('fivethirtyeight')
g = sns.barplot(y,x,palette=("binary"))
g.grid(False)
g.set_title(title)
g.set(xlabel='Configurations', ylabel='Expansion')
ax = g
for i, v in enumerate(y):
    ax.text(v + 1, i + .25, str(v), color='black', fontweight = "bold")
fig_name = dir_list[1]
plt.savefig(fig_name, bbox_inches='tight')
plt.show(block=False)
plt.pause(1)
plt.close()

y = y2
x = x2
title = "Configurations vs Different Expansions" + " for " + dir_list[2]
print(title)
sns.set(font_scale=1)
plt.rcParams['figure.figsize'] = (8,4)
plt.rcParams['font.family'] = "serif"
style.use('fivethirtyeight')
g = sns.barplot(y,x,palette=("binary"))
g.grid(False)
g.set_title(title)
g.set(xlabel='Configurations', ylabel='Expansion')
ax = g
for i, v in enumerate(y):
    ax.text(v + 1, i + .25, str(v), color='black', fontweight = "bold")
fig_name = dir_list[2]
plt.savefig(fig_name, bbox_inches='tight')
plt.show(block=False)
plt.pause(1)
plt.close()

y = y3
x = x3
title = "Configurations vs Different Expansions" + " for " + dir_list[3]
print(title)
sns.set(font_scale=1)
plt.rcParams['figure.figsize'] = (8,4)
plt.rcParams['font.family'] = "serif"
style.use('fivethirtyeight')
g = sns.barplot(y,x,palette=("binary"))
g.grid(False)
g.set_title(title)
g.set(xlabel='Configurations', ylabel='Expansion')
ax = g
for i, v in enumerate(y):
    ax.text(v + 1, i + .25, str(v), color='black', fontweight = "bold")
fig_name = dir_list[3]
plt.savefig(fig_name, bbox_inches='tight')
plt.show(block=False)
plt.pause(1)
plt.close()

y = y4
x = x4
title = "Configurations vs Different Expansions" + " for " + dir_list[4]
print(title)
sns.set(font_scale=1)
plt.rcParams['figure.figsize'] = (8,4)
plt.rcParams['font.family'] = "serif"
style.use('fivethirtyeight')
g = sns.barplot(y,x,palette=("binary"))
g.grid(False)
g.set_title(title)
g.set(xlabel='Configurations', ylabel='Expansion')
ax = g
for i, v in enumerate(y):
    ax.text(v + 1, i + .25, str(v), color='black', fontweight = "bold")
fig_name = dir_list[4]
plt.savefig(fig_name, bbox_inches='tight')
plt.show(block=False)
plt.pause(1)
plt.close()

y = y5
x = x5
title = "Configurations vs Different Expansions" + " for " + dir_list[5]
print(title)
sns.set(font_scale=1)
plt.rcParams['figure.figsize'] = (8,4)
plt.rcParams['font.family'] = "serif"
style.use('fivethirtyeight')
g = sns.barplot(y,x,palette=("binary"))
g.grid(False)
g.set_title(title)
g.set(xlabel='Configurations', ylabel='Expansion')
ax = g
for i, v in enumerate(y):
    ax.text(v + 1, i + .25, str(v), color='black', fontweight = "bold")
fig_name = dir_list[5]
plt.savefig(fig_name, bbox_inches='tight')
plt.show(block=False)
plt.pause(1)
plt.close()

rcParams['figure.figsize'] = 30,20
plt.rcParams['axes.grid'] = False
img_1 = mpimg.imread('dihedral_threshold_lower.png')
img_2 = mpimg.imread('dihedral_threshold_upper.png')
img_3 = mpimg.imread('dual_threshold_lower.png')
img_4 = mpimg.imread('dual_threshold_upper.png')
img_5 = mpimg.imread('total_threshold_lower.png')
img_6 = mpimg.imread('total_threshold_upper.png')
fig, ax = plt.subplots(3,2)
fig.suptitle('')
ax[0,1].imshow(img_1);
ax[1,1].imshow(img_2);
ax[0,0].imshow(img_3);
ax[1,0].imshow(img_4);
ax[2,0].imshow(img_5);
ax[2,1].imshow(img_6);
plt.savefig('analysis.jpeg')
plt.show(block=False)
plt.pause(3)
plt.close()

os.system("rm -rf *png*")
