import os
import shutil
dir_list = ["dihedral_threshold_lower", "dihedral_threshold_upper", "dual_threshold_lower", "dual_threshold_upper", "total_threshold_lower", "total_threshold_upper"]
cwd = os.getcwd()
source_dir = cwd + "/"
target_dir = cwd + "/" + "gamd_simulations" + "/"
for i in dir_list:
    os.chdir(target_dir + i)
    os.system("rm -rf pickle_files dat_files txt_csv_files")
    os.chdir(cwd)
for i in dir_list:
    shutil.move(cwd + "/" + "gamd_simulations" + "/" + i  + "/" + "system_files" + "/" + "gamd.log", cwd + "/" + "gamd_simulations" + "/" + i  + "/" + "gamd.log")
    shutil.move(cwd + "/" + "gamd_simulations" + "/" + i  + "/" + "system_files" + "/" + "gamd-restart.dat", cwd + "/" + "gamd_simulations" + "/" + i  + "/" + "gamd-restart.dat")
    shutil.move(cwd + "/" + "gamd_simulations" + "/" + i  + "/" + "system_files" + "/" + "md.in", cwd + "/" + "gamd_simulations" + "/" + i  + "/" + "md.in")
    shutil.move(cwd + "/" + "gamd_simulations" + "/" + i  + "/" + "system_files" + "/" + "mdinfo", cwd + "/" + "gamd_simulations" + "/" + i  + "/" + "mdinfo")
    shutil.move(cwd + "/" + "gamd_simulations" + "/" + i  + "/" + "system_files" + "/" + "system_final.inpcrd", cwd + "/" + "gamd_simulations" + "/" + i  + "/" + "system_final.inpcrd")
    shutil.move(cwd + "/" + "gamd_simulations" + "/" + i  + "/" + "system_files" + "/" + "system_final.nc", cwd + "/" + "gamd_simulations" + "/" + i  + "/" + "system_final.nc")
    shutil.move(cwd + "/" + "gamd_simulations" + "/" + i  + "/" + "system_files" + "/" + "system_final.out", cwd + "/" + "gamd_simulations" + "/" + i  + "/" + "system_final.out")
    shutil.move(cwd + "/" + "gamd_simulations" + "/" + i  + "/" + "system_files" + "/" + "system_final.prmtop", cwd + "/" + "gamd_simulations" + "/" + i  + "/" + "system_final.prmtop")
    shutil.move(cwd + "/" + "gamd_simulations" + "/" + i  + "/" + "system_files" + "/" + "system_final.rst", cwd + "/" + "gamd_simulations" + "/" + i  + "/" + "system_final.rst")
for i in dir_list:
    os.chdir(target_dir + i)
    os.system("rm -rf system_files")
    os.chdir(cwd)

