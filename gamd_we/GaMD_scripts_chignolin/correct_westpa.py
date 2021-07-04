import os
import re
import shutil
import fnmatch
import pandas as pd
cwd = os.getcwd()
os.chdir(cwd + "/" + "westpa_dir")
dir_list = ["dihedral_threshold_lower", "dihedral_threshold_upper", "dual_threshold_lower", "dual_threshold_upper", "total_threshold_lower", "total_threshold_upper"]
we_list = ["1d_c1", "1d_c12", "1d_c123", "2d_c1", "2d_c12", "2d_c123"]
for i in dir_list:
    for j in we_list:
        os.chdir(cwd + "/" + "westpa_dir" + "/" + str(i))
        os.chdir(cwd + "/" + "westpa_dir" + "/" + str(i) + "/" + str(j))
        if len(open("BASIS_STATES").readlines(  )) > 0 :    
            os.chdir (cwd + "/" + "westpa_dir" + "/" + str(i) + "/" + str(j) + "/" + "md_sims")
            files = os.listdir('.')
            file_to_find = "*.out"
            out_list = []
            for y in files:
                if fnmatch.fnmatch(y, file_to_find):
                    out_list.append(y)
            list_failed_jobs = []
            for out_file in out_list:
                with open(out_file, 'r') as f:
                    last_line = f.readlines()[-2]
                    if last_line.startswith('|') == False:
                        list_failed_jobs.append(out_file)
            for c in range(len(list_failed_jobs)):
                command = "rm -rf " + list_failed_jobs[c]
                os.system(command)
            for d in range(len(list_failed_jobs)):
                command = "rm -rf " + list_failed_jobs[d][:-3] + "rst"
                os.system(command)
            for e in range(len(list_failed_jobs)):
                command = "rm -rf " + list_failed_jobs[e][:-3] + "inpcrd"
                os.system(command)
            for f in range(len(list_failed_jobs)):
                command = "rm -rf " + list_failed_jobs[f][:-3] + "nc"
                os.system(command) 
                   
            files = os.listdir('.')
            file_to_find = "*.rst"
            rst_list = []
            for y in files:
                if fnmatch.fnmatch(y, file_to_find):
                    rst_list.append(y)
            rst_failed_jobs = []
            for rst_file in rst_list:
                with open(rst_file, 'r') as f:
                    req_line = f.readlines()[2]
                    if "NaN" in req_line:
                        rst_failed_jobs.append(rst_file)
            for g in range(len(rst_failed_jobs)):
                command = "rm -rf " + rst_failed_jobs[g]
                os.system(command)
            for h in range(len(rst_failed_jobs)):
                command = "rm -rf " + rst_failed_jobs[h][:-3] + "out"
                os.system(command)
            for u in range(len(rst_failed_jobs)):
                command = "rm -rf " + rst_failed_jobs[u][:-3] + "inpcrd"
                os.system(command)
            for v in range(len(rst_failed_jobs)):
                command = "rm -rf " + rst_failed_jobs[v][:-3] + "nc"
                os.system(command)
            
            files_2 = os.listdir('.')
            file_to_find_2 = "*.rst"
            rst_list_2 = []
            for y in files_2:
                if fnmatch.fnmatch(y, file_to_find_2):
                    rst_list_2.append(y)
            rst_failed_jobs_2 = []
            for rst_file_2 in rst_list_2:
                with open(rst_file_2, 'r') as f:
                    lines_file = f.readlines()
                    for req_line in lines_file:
                        if "*" in req_line:
                            rst_failed_jobs_2.append(rst_file_2)
            for g in range(len(rst_failed_jobs_2)):
                command = "rm -rf " + rst_failed_jobs_2[g]
                os.system(command)
            for h in range(len(rst_failed_jobs_2)):
                command = "rm -rf " + rst_failed_jobs_2[h][:-3] + "out"
                os.system(command)
            for u in range(len(rst_failed_jobs_2)):
                command = "rm -rf " + rst_failed_jobs_2[u][:-3] + "inpcrd"
                os.system(command)
            for v in range(len(rst_failed_jobs_2)):
                command = "rm -rf " + rst_failed_jobs_2[v][:-3] + "nc"
                os.system(command)
            
            os.system("rm -rf md.in")
            os.system("rm -rf system_final.prmtop")
            os.system("rm -rf mdinfo")  
            files = os.listdir('.')
            inpcrd_file_to_find = "*.inpcrd"
            rst_file_to_find = "*.rst"
            inpcrd_file_list = []
            for y in files:
                if fnmatch.fnmatch(y, inpcrd_file_to_find):
                    inpcrd_file_list.append(y)
            rst_file_list = []
            for z in files:
                if fnmatch.fnmatch(z, rst_file_to_find):
                    rst_file_list.append(z)
            os.chdir (cwd + "/" + "westpa_dir" + "/" + str(i) + "/" + str(j))
            os.system("rm -rf bstates_corrected_rst")
            os.system("mkdir bstates_corrected_rst")
            os.system("rm -rf bstates_corrected_inpcrd")
            os.system("mkdir bstates_corrected_inpcrd")
            for x in inpcrd_file_list:
                shutil.copy(cwd + "/" + "westpa_dir" + "/" + str(i) + "/" + str(j) + "/" + "md_sims" + "/" + str(x), cwd + "/" + "westpa_dir" + "/" + str(i) + "/" + str(j) + "/" + "bstates_corrected_inpcrd" + "/" + str(x))
            for y in rst_file_list:
                shutil.copy(cwd + "/" + "westpa_dir" + "/" + str(i) + "/" + str(j) + "/" + "md_sims" + "/" + str(y), cwd + "/" + "westpa_dir" + "/" + str(i) + "/" + str(j) + "/" + "bstates_corrected_rst" + "/" + str(y))
            df = pd.read_csv('BASIS_STATES', sep=" ", header = None)  
            df.columns = ["index_df", "probability", "inpcrd"]
            df = df[["probability", "inpcrd"]]
            df = df[df.inpcrd.str.contains('|'.join(inpcrd_file_list))]
            index_row_list = []
            for n in range(df.shape[0]):
                index_row_list.append(n)
            df = df.assign(index_ = index_row_list) 
            df = df[['index_', 'probability', 'inpcrd']]
            df.to_csv('BASIS_STATES_CORRECTED_INPCRD', header=False, index=None, sep=" " ,mode='w') 
            fin = open("BASIS_STATES_CORRECTED_INPCRD", "rt")
            fout = open("BASIS_STATES_CORRECTED_RST", "wt")
            for line in fin:
                fout.write(line.replace('inpcrd', 'rst'))
            fin.close()
            fout.close()
