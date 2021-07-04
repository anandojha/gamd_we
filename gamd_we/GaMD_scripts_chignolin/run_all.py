import os
# module load amber
# module load cuda-10.1
os.system("python system_preparation_I.py")
os.system("python equilibration.py")
os.system("python system_preparation_II.py")
os.system("python add_vectors_inpcrd.py")
os.system("python add_vectors_prmtop.py")
os.system("python filetree.py")
os.system("python run_simulations.py")
# module unload amber
# module unload cuda-10.1
os.system("python reweight.py")
# module load amber
# module load cuda-10.1
os.system("python save_westpa_inputs.py")
os.system("python clean_transfer.py")
os.system("python add_vectors_westpa.py")
os.system("python run_we_analysis.py")
os.system("python correct_westpa.py")
os.system("python analysis.py")
# To run the entire analysis again,run the following file
os.system("python clean_reanalysis.py")
