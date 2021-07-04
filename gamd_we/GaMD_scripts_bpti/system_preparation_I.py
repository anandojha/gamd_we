import os
import shutil
#  Go to http://ambermd.org/tutorials/advanced/tutorial22/files/5PTI-DtoH-dry.pdb and download : 5PTI-DtoH-dry.pdb
os.system("curl -O http://ambermd.org/tutorials/advanced/tutorial22/files/5PTI-DtoH-dry.pdb")
os.system("rm -rf system_inputs")
os.system("mkdir system_inputs")
cwd = os.getcwd()
target_dir = cwd + "/" + "system_inputs"
# save the tleap script to file
with open('input_TIP4P.leap', 'w') as f:
    f.write('''
source leaprc.protein.ff14SB
loadOff solvents.lib
loadOff tip4pbox.off
loadOff tip4pewbox.off
source leaprc.water.tip4pew
HOH = TP4
pdb = loadpdb 5PTI-DtoH-dry.pdb
bond pdb.55.SG pdb.5.SG
bond pdb.30.SG pdb.51.SG
bond pdb.14.SG pdb.38.SG
charge pdb
addions2 pdb Cl- 6
charge pdb
solvatebox pdb TIP4PEWBOX 12.0
saveamberparm pdb system_TIP4P.prmtop system_TIP4P.inpcrd
saveamberparm pdb system_TIP4P.parm7 system_TIP4P.rst7
savepdb pdb system_TIP4P.pdb
quit
''')
os.system("tleap -f input_TIP4P.leap")
os.system("rm -rf leap.log")
shutil.copy(cwd + "/" + "system_TIP4P.inpcrd" , target_dir + "/" + "system_TIP4P.inpcrd")
shutil.copy(cwd + "/" + "system_TIP4P.parm7" , target_dir + "/" + "system_TIP4P.parm7")
shutil.copy(cwd + "/" + "system_TIP4P.pdb" , target_dir + "/" + "system_TIP4P.pdb")
shutil.copy(cwd + "/" + "system_TIP4P.prmtop" , target_dir + "/" + "system_TIP4P.prmtop")
shutil.copy(cwd + "/" + "system_TIP4P.rst7" , target_dir + "/" + "system_TIP4P.rst7")
shutil.copy(cwd + "/" + "input_TIP4P.leap" , target_dir + "/" + "input_TIP4P.leap")
shutil.copy(cwd + "/" + "5PTI-DtoH-dry.pdb" , target_dir + "/" + "5PTI-DtoH-dry.pdb")
os.system("rm -rf system_TIP4P.inpcrd")
os.system("rm -rf system_TIP4P.parm7")
os.system("rm -rf system_TIP4P.pdb")
os.system("rm -rf system_TIP4P.rst7")
os.system("rm -rf system_TIP4P.prmtop")
os.system("rm -rf input_TIP4P.leap")
os.system("rm -rf 5PTI-DtoH-dry.pdb")
