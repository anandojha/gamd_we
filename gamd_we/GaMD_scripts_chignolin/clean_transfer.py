import os
import shutil
os.system("rm -rf westpa_dir")
os.system("mkdir westpa_dir")
dir_list = ["dihedral_threshold_lower", "dihedral_threshold_upper", "dual_threshold_lower", "dual_threshold_upper", "total_threshold_lower", "total_threshold_upper"]
cwd = os.getcwd()
source_dir = cwd + "/"
target_dir = cwd + "/" + "gamd_simulations" + "/"
for i in dir_list:
    os.chdir(target_dir + i)
    os.system("rm -rf  extract_data.py")
    os.system("rm -rf reweight_ce_1d.py")
    os.system("rm -rf reweight_ce_2d.py")
    os.system("rm -rf save_frames.py")
    os.system("rm -rf save_we_input.py")
    os.system("rm -rf arrange_files.py")
    os.system("rm -rf westpa_input.py")
for i in dir_list:
    os.chdir(source_dir + "westpa_dir")
    command = "mkdir {}".format(i)
    os.system(command)
    os.chdir(cwd)
for i in dir_list:
    shutil.copytree(target_dir + i + "/" + "westpa_inputs", source_dir + "westpa_dir" + "/" + i  + "/"  "westpa_inputs")
    we_list = ["1d_c1","1d_c12","1d_c123","2d_c1","2d_c12","2d_c123"]
    for j in we_list:
        shutil.copytree(source_dir + "westpa_dir" + "/"  + i + "/" + "westpa_inputs" + "/" + j, source_dir + "westpa_dir" + "/" + i + "/" + j)
    dest_dir = source_dir + "westpa_dir" + "/"  + i 
    os.chdir(dest_dir)
    os.system("rm -rf westpa_inputs")
    os.chdir(cwd)
