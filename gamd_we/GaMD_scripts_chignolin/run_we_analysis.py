import pandas as pd
import itertools
import shutil
import os
cwd = os.getcwd()
os.chdir(cwd + "/" + "westpa_dir")
dir_list = ["dihedral_threshold_lower", "dihedral_threshold_upper", "dual_threshold_lower", "dual_threshold_upper", "total_threshold_lower", "total_threshold_upper"]
we_list = ["1d_c1", "1d_c12", "1d_c123", "2d_c1", "2d_c12", "2d_c123"]
for i in dir_list:
    for j in we_list:
        os.chdir(cwd + "/" + "westpa_dir" + "/" + str(i))
        os.chdir(cwd + "/" + "westpa_dir" + "/" + str(i) + "/" + str(j))
        if len(open("BASIS_STATES").readlines(  )) > 0 :
            df = pd.read_csv("BASIS_STATES", delimiter = ' ', header = None)
            df.columns = [["descriptor", "probability", "file_name"]]
            df1 = df[["file_name"]]
            inpcrd_list = df1.values.tolist()
            inpcrd_list = list(itertools.chain(*inpcrd_list))
            os.system("rm -rf md_sims")
            os.system("mkdir md_sims")
            os.chdir(cwd + "/" + "westpa_dir" + "/" + str(i)+ "/" + str(j) + "/" + "md_sims")
            with open('md.in', 'w') as f:
                f.write("Run minimization followed by saving rst file"                + '\n')
                f.write("&cntrl"                                                      + '\n')
                f.write("  imin = 1, maxcyc = 10000, ntpr = 5, iwrap = 1, ntxo = 1"   + '\n')
                f.write("&end"                                                        + '\n')
            for k in inpcrd_list:
                source_dir = cwd + "/" + "westpa_dir" + "/" + str(i) + "/" + str(j) + "/" + "bstates"
                target_dir = cwd + "/" + "westpa_dir" + "/" + str(i) + "/" + str(j) + "/" + "md_sims"
                shutil.copy(source_dir + "/" + str(k) , target_dir + "/" + str(k))
            source_dir = cwd + "/" + "starting_structures"
            target_dir = cwd + "/" + "westpa_dir" + "/" + str(i) + "/" + str(j) + "/" + "md_sims"
            shutil.copy(source_dir + "/" + "system_final.prmtop", target_dir + "/" + "system_final.prmtop")
            for l in range(len(inpcrd_list)):
                command = "pmemd.cuda -O -i md.in -o " + inpcrd_list[l][:-6] + "out" + " -p system_final.prmtop -c " + inpcrd_list[l] + " -r " + inpcrd_list[l][:-6] + "rst"
                print(command)
                os.system(command)
