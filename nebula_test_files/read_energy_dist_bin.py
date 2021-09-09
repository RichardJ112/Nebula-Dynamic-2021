import numpy as np
import os
import sys
import pickle
from statistics import median 
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import rc
from mpl_toolkits import mplot3d
import time
import random
from matplotlib import cm
from matplotlib import animation
from struct import *
import math
import subprocess
import ffmpeg

#Parameters for Energy Distributions

axis_range = 5 # change +- x/y limits in nm
desktop_toggle = True

if desktop_toggle:
    path_to_test_files = "C:/Users/Richard/source/repos/Nebula/nebula_test_files/"



# IO parameters
date = "/September/1_9_2021/"
deposit_type = 'point'#types supported are line and point
material = 'tungsten' #other important material are : tungsten,silicon,graphite-phonon
cross_section_type = '_SM_' #other cross-sections include: Alman-Winters(_AM_), Winters(_W_), Alman (_A_) and Smith (_SM_) 
parameter_summary = "1keV_1_1kpp_pitch_0_401_401_1001_sb_1000_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_"
#extra = "no_deposition/"
#extra = "deposition/"
extra = "_2s"
input_path = path_to_test_files + "output"+date+parameter_summary+material+cross_section_type+"energy_dist"+extra+".bin"

#Parameter Retrieval
path_list = (input_path.split('/')[-1]).split('.')
title =path_list[0]
title_list = title.split('_')
energy_str = title_list[0]
pillar_num = int(title_list[1])
pitch = int(title_list[4])
electron_num = pillar_num * int(title_list[2].replace('kpp',''))*1000


output_path_folder = path_to_test_files+"figures"+date+title+"/"
general_path = output_path_folder+title
fignum = 1

if not os.path.exists(output_path_folder):
    os.makedirs(output_path_folder)

# General path strings
general_energy_dist_str = general_path+"_general_energy_dist"

#Text Loading
tic =  time.perf_counter() 
print("Loading file...")

with open(input_path, mode='rb') as file: # b is important -> binary
   fileContent = file.read()

#Binary Unpacking ----------------------------------------------
lengrid = unpack("l", fileContent[:4])[0]
counter = 4

energies = unpack("f"*lengrid,fileContent[counter:counter+lengrid*4])
counter+=lengrid*4
x = unpack("f"*lengrid,fileContent[counter:counter+lengrid*4])
counter+=lengrid*4
y = unpack("f"*lengrid,fileContent[counter:counter+lengrid*4])
counter+=lengrid*4
z = unpack("f"*lengrid,fileContent[counter:counter+lengrid*4])
counter+=lengrid*4

toc = time.perf_counter()
print(f"Finished Loading in {toc - tic:0.4f} seconds\n")

#The stored bytes in order c++ ---------------------------------------
""" 
//Length 
std::ofstream output_bin_file(file_name, std::ios::binary);
int64_t len = boundary_energy_dist.boundary_energies.size();
output_bin_file.write( (char*)&len, sizeof(len) );

//Vectors
output_bin_file.write( (char*)&boundary_energy_dist.boundary_energies[0], len * sizeof(float_t) );
output_bin_file.write( (char*)&boundary_energy_dist.x_pos[0], len * sizeof(float_t) );
output_bin_file.write( (char*)&boundary_energy_dist.y_pos[0], len * sizeof(float_t) );
output_bin_file.write( (char*)&boundary_energy_dist.z_pos[0], len * sizeof(float_t) );

"""

#Data Manipulation --------------------------------------------------

z_depth = int(title_list[7])
zero_point = int(title_list[6])*0.3/2
energies_array = np.asarray(energies)

#Electron Efficiency
substrate_se_electrons = np.count_nonzero(energies_array < 50)
substrate_bse_electrons = np.count_nonzero(energies_array > 50)
se_efficiency = round(substrate_se_electrons/electron_num,2)
bse_efficiency = round(substrate_bse_electrons/electron_num,2)

print(str(substrate_se_electrons) + " electrons under 50eV detected")
print(str(substrate_bse_electrons) + " electrons above 50eV detected")
print("The simulated SE efficiency is : " + str(se_efficiency))
print(str(se_efficiency))
print("The simulated BSE efficiency is : " + str(bse_efficiency))
print(str(bse_efficiency))


#Timing
tic =  time.perf_counter() 
print("Sorting into Primary Cascade Groups...")

x_array = np.asarray(x)
y_array = np.asarray(y)
z_array = np.asarray(z)

#Index Finding
min_index = int(float(title_list[6])/2-axis_range/0.3)
max_index = int(float(title_list[6])/2+axis_range/0.3)

# Axis Shifting
x_array = zero_point-x_array
y_array = zero_point-y_array
z_array = (z_depth*0.3)-z_array
    
print("Creating Plots...")
tic = time.perf_counter()

#General Energy Distribution ------------------------------
fig = plt.figure(fignum)
fignum+=1
low_energies = [x for x in energies_array if x <= 100]
#density = stats.gaussian_kde(energies_array)
density = stats.gaussian_kde(low_energies)
n, x, _ = plt.hist(low_energies, bins=500,histtype=u'step', density=False)    
plt.xlabel('Energy(eV)')
plt.xlim(left = 0)
plt.ylabel('(# of electrons)')

fig.savefig(general_energy_dist_str, dpi = 800, bbox_inches="tight")

toc = time.perf_counter()
print(f"Finished  Plots {toc - tic:0.4f} seconds")
