import os
import shutil
#  Go to https://files.rcsb.org/download/1UAO.pdb1.gz and download : 1UAO.pdb1.gz
os.system("curl -O https://files.rcsb.org/download/1UAO.pdb1.gz")
os.system("gunzip 1UAO.pdb1.gz")
os.system("mv 1UAO.pdb1 chignolin.pdb")
os.system("rm -rf system_inputs")
os.system("mkdir system_inputs")
cwd = os.getcwd()
target_dir = cwd + "/" + "system_inputs"
os.system("pdb4amber -i chignolin.pdb -o system.pdb")
# save the tleap script to file
with open('input_TIP3P.leap', 'w') as f:
    f.write('''
source leaprc.protein.ff14SB
source leaprc.water.tip3p
set default FlexibleWater on
set default PBRadii mbondi2
pdb = loadpdb system.pdb
solvateBox pdb TIP3PBOX 15
charge pdb
addions2 pdb Na+ 2
charge pdb
saveamberparm pdb system_TIP3P.prmtop system_TIP3P.inpcrd
saveamberparm pdb system_TIP3P.parm7 system_TIP3P.rst7
savepdb pdb system_TIP3P.pdb
quit
''') 
os.system("tleap -f input_TIP3P.leap")
os.system("rm -rf leap.log")
shutil.copy(cwd + "/" + "system_TIP3P.inpcrd" , target_dir + "/" + "system_TIP3P.inpcrd")
shutil.copy(cwd + "/" + "system_TIP3P.parm7" , target_dir + "/" + "system_TIP3P.parm7")
shutil.copy(cwd + "/" + "system_TIP3P.pdb" , target_dir + "/" + "system_TIP3P.pdb")
shutil.copy(cwd + "/" + "system_TIP3P.prmtop" , target_dir + "/" + "system_TIP3P.prmtop")
shutil.copy(cwd + "/" + "system_TIP3P.rst7" , target_dir + "/" + "system_TIP3P.rst7")
shutil.copy(cwd + "/" + "system.pdb", target_dir + "/" + "system.pdb")
shutil.copy(cwd + "/" + "input_TIP3P.leap" , target_dir + "/" + "input_TIP3P.leap")
shutil.copy(cwd + "/" + "chignolin.pdb" , target_dir + "/" + "chignolin.pdb")
os.system("rm -rf system_sslink")
os.system("rm -rf system_nonprot.pdb")
os.system("rm -rf system.pdb")
os.system("rm -rf system_renum.txt")
os.system("rm -rf system_TIP3P.inpcrd")
os.system("rm -rf system_TIP3P.parm7")
os.system("rm -rf system_TIP3P.pdb")
os.system("rm -rf system_TIP3P.inpcrd")
os.system("rm -rf system_TIP3P.rst7")
os.system("rm -rf system_TIP3P.prmtop")
os.system("rm -rf system.pdb")
os.system("rm -rf input_TIP3P.leap")
os.system("rm -rf chignolin.pdb")
