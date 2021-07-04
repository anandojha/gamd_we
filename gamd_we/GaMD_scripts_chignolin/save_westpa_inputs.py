import os
import shutil
dir_list = ["dihedral_threshold_lower", "dihedral_threshold_upper", "dual_threshold_lower", "dual_threshold_upper", "total_threshold_lower", "total_threshold_upper"]
cwd = os.getcwd()
source_dir = cwd + "/" 
target_dir = cwd + "/" + "gamd_simulations" + "/" 
#copy relevant python files to simulation folders
for i in dir_list:
    shutil.copy(source_dir + "westpa_input.py" , target_dir + i + "/" + "westpa_input.py")
for i in dir_list:
    os.chdir(target_dir + i)
    os.system("python westpa_input.py")
