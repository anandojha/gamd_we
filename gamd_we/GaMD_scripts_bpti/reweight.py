import os
import shutil
dir_list = ["dihedral_threshold_lower", "dihedral_threshold_upper", "dual_threshold_lower", "dual_threshold_upper", "total_threshold_lower", "total_threshold_upper"]
cwd = os.getcwd()
source_dir = cwd + "/" 
target_dir = cwd + "/" + "gamd_simulations" + "/" 
#copy relevant python files to simulation folders
for i in dir_list:
    shutil.copy(source_dir + "extract_data.py" , target_dir + i + "/" + "extract_data.py")
    shutil.copy(source_dir + "reweight_ce_1d.py" , target_dir + i + "/" + "reweight_ce_1d.py")
    shutil.copy(source_dir + "reweight_ce_2d.py" , target_dir + i + "/" + "reweight_ce_2d.py")
    shutil.copy(source_dir + "arrange_files.py" , target_dir + i + "/" + "arrange_files.py")
#run reweighting and analysis in each of the simulation folder    
for i in dir_list:
    os.chdir(target_dir + i)
    os.system("python extract_data.py")
    os.system("python reweight_ce_1d.py")
    os.system("python reweight_ce_2d.py")
    os.system("python arrange_files.py")
