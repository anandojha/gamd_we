import os
import fnmatch
cwd = os.getcwd()
source_dir = cwd
westpa_dir = cwd + "/" + "westpa_dir"
os.chdir(source_dir + "/" + "starting_structures")
with open('system_final.inpcrd') as f:
    for line in f:
        pass
    vector_information = line
print(vector_information)
os.chdir(source_dir)
dir_list = ["dihedral_threshold_lower", "dihedral_threshold_upper", "dual_threshold_lower",
            "dual_threshold_upper", "total_threshold_lower", "total_threshold_upper"]
we_list = ["1d_c1","1d_c12","1d_c123","2d_c1","2d_c12","2d_c123"]
for i in dir_list:
    os.chdir(westpa_dir + "/" + str(i))
    for j in we_list:
        os.chdir(westpa_dir + "/" + str(i) + "/" + str(j) + "/" + "bstates")
        files = os.listdir('.')
        file_to_find = "*.inpcrd"
        inpcrd_list = []
        for k in files:
            if fnmatch.fnmatch(k, file_to_find):
                inpcrd_list.append(k)  
        for l in inpcrd_list:
            with open(l, "a+") as f:
                f.write(vector_information)
