import os
import shutil
from input_parameters import *
cwd = os.getcwd()
os.system("rm -rf gamd_simulations")
os.system("mkdir gamd_simulations")
os.chdir (cwd + "/" + "gamd_simulations")
source_dir = cwd + "/" + "starting_structures" 
target_dir = cwd + "/" + "gamd_simulations" 
dir_list = ["dihedral_threshold_lower", "dihedral_threshold_upper", "dual_threshold_lower", "dual_threshold_upper", "total_threshold_lower", "total_threshold_upper"] 
for i in range(len(dir_list)):
    os.mkdir (dir_list[i])
    os.chdir(target_dir + "/" + dir_list[i])
    shutil.copy(source_dir + "/" + "system_final.inpcrd", target_dir + "/" + dir_list[i] + "/" + "system_final.inpcrd")
    shutil.copy(source_dir + "/" + "system_final.prmtop", target_dir + "/" + dir_list[i] + "/" + "system_final.prmtop")
    #shutil.copy(cwd + "/"  + "input_parameters.py", target_dir + "input_parameters.py")
    if "lower" in dir_list[i]:
        i_E = 1
    if "upper" in dir_list[i]:
        i_E = 2   
    if "total" in dir_list[i]:
        i_gamd = 1
    if "dihedral" in dir_list[i]:
        i_gamd = 2
    if "dual" in dir_list[i]:
        i_gamd = 3
    with open('md.in', 'w') as f:
        f.write("&cntrl"                                                                                                 + '\n')
        f.write("  imin = 0, irest = 0, ntx = 1,"                                                                        + '\n')
        f.write("  nstlim = " + str(nst_lim) + ", dt = 0.002,"                                                           + '\n')
        f.write("  ntc = 2, ntf = 2, tol = 0.000001,"                                                                    + '\n')
        f.write("  iwrap = 1, ntb = 1, cut = 8.0,"                                                                       + '\n')
        f.write("  ntt = 3, temp0 = 300.0, gamma_ln = 1.0, "                                                             + '\n')
        f.write("  ntpr = 500, ntwx = " + str(ntw_x) + ", ntwr = 500,"                                                   + '\n')    
        f.write("  ntxo = 2, ioutfm = 1, ig = -1, ntwprt = 0,"                                                           + '\n')    
        f.write("  igamd = " + str(i_gamd) + ", iE = " + str(i_E) + ", irest_gamd = 0,"                                  + '\n')    
        f.write("  ntcmd = " + str(nt_cmd) + ", nteb = " + str(n_teb) + ", ntave = " + str(n_tave) + ","                 + '\n')   
        f.write("  ntcmdprep = " + str(ntcmd_prep) + ", ntebprep = " + str(nteb_prep) + ","                              + '\n')    
        f.write("  sigma0D = 6.0, sigma0P = 6.0"                                                                         +' \n')    
        f.write("&end"                                                                                                   + '\n')  
    os.chdir(target_dir)
