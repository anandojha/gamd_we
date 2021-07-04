import os
import fileinput
import shutil
cwd = os.getcwd()
target_dir = cwd + "/" + "starting_structures"
os.system("rm -rf starting_structures")
os.system("mkdir starting_structures")
shutil.copy(cwd  + "/" + "equilibration" + "/" + "system_nvt_output_last_frame.pdb" ,  target_dir + "/" + "system_nvt_output_last_frame.pdb")
os.chdir(target_dir)
os.system("pdb4amber -i system_nvt_output_last_frame.pdb -o intermediate_temp.pdb")
os.system("rm -rf intermediate_temp_renum.txt")
os.system("rm -rf intermediate_temp_sslink")
os.system("rm -rf intermediate_temp_nonprot.pdb")
remove_words = ['H   ARG A   1']
with open('intermediate_temp.pdb') as oldfile, open('intermediate.pdb', 'w') as newfile:
    for line in oldfile:
        if not any(word in line for word in remove_words):
            newfile.write(line)
# Save the tleap script to file
with open('final_input_TIP4P.leap', 'w') as f:
    f.write('''
source leaprc.protein.ff14SB
source leaprc.water.tip4pew
pdb = loadpdb intermediate.pdb
charge pdb
saveamberparm pdb system_final.prmtop system_final.inpcrd
saveamberparm pdb system_final.parm7 system_final.rst7
savepdb pdb system_final.pdb
quit
''')
os.system("tleap -f final_input_TIP4P.leap")
os.system("rm -rf leap.log")
os.system("rm -rf leap.log")
os.system("rm -rf intermediate.pdb")
os.system("rm -rf intermediate_temp.pdb")
os.system("rm -rf system_nvt_output_last_frame.pdb")
os.chdir(cwd)
