import os
import shutil
cwd = os.getcwd()
target_dir = cwd + "/" + "equilibration"
os.system("rm -rf equilibration")
os.system("mkdir equilibration")
shutil.copy(cwd + "/" + "system_inputs" + "/" + "system_TIP3P.inpcrd" , target_dir + "/" + "system_TIP3P.inpcrd")
shutil.copy(cwd + "/" + "system_inputs" + "/" + "system_TIP3P.parm7" , target_dir + "/" + "system_TIP3P.parm7")
shutil.copy(cwd + "/" + "system_inputs" + "/" + "system_TIP3P.pdb" , target_dir + "/" + "system_TIP3P.pdb")
shutil.copy(cwd + "/" + "system_inputs" + "/" + "system_TIP3P.prmtop" , target_dir + "/" + "system_TIP3P.prmtop")
shutil.copy(cwd + "/" + "system_inputs" + "/" + "system_TIP3P.rst7" , target_dir + "/" + "system_TIP3P.rst7")
shutil.copy(cwd + "/" + "system_inputs" + "/" + "system.pdb", target_dir + "/" + "system.pdb")
shutil.copy(cwd + "/" + "system_inputs" + "/" +  "chignolin.pdb" , target_dir + "/" + "chignolin.pdb")
shutil.copy(cwd + "/" + "system_inputs" + "/" + "input_TIP3P.leap" , target_dir + "/" + "input_TIP3P.leap")
shutil.copy(cwd + "/" + "input_parameters.py" , target_dir + "/" + "input_parameters.py")
shutil.copy(cwd + "/" + "equilibration_I.py" , target_dir + "/" + "equilibration_I.py")
shutil.copy(cwd + "/" + "equilibration_II.py" , target_dir + "/" + "equilibration_II.py")
shutil.copy(cwd + "/" + "equilibration_III.py" , target_dir + "/" + "equilibration_III.py")
os.chdir(target_dir)
os.system("python equilibration_I.py > equilibration_I.out")
os.system("python equilibration_II.py > equilibration_II.out")
os.system("python equilibration_III.py > equilibration_III.out")
os.system("rm -rf system_TIP3P.inpcrd")
os.system("rm -rf system_TIP3P.parm7")
os.system("rm -rf system_TIP3P.pdb")
os.system("rm -rf system_TIP3P.inpcrd")
os.system("rm -rf system_TIP3P.rst7")
os.system("rm -rf system_TIP3P.prmtop")
os.system("rm -rf system.pdb")
os.system("rm -rf chignolin.pdb")
os.system("rm -rf input_TIP3P.leap")
os.system("rm -rf equilibration_I.py")
os.system("rm -rf equilibration_II.py")
os.system("rm -rf equilibration_III.py")
os.system("rm -rf input_parameters.py")
os.chdir(cwd)
