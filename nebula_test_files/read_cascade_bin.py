import numpy as np
import os
import sys
import pickle
from statistics import median 
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

# Version Looking at Overall Deep Scattering 2D
# External Function

#Parameters Cascade Diagram

primaries_to_plot = 10
f_range = 60 # change +- x/y limits in nm (this is where mirroring occurs)
z_roi_limits = [-10,400] #limit of interest z-direction in nm
overlay_toggle = True # toggle inclusion of simple geometric overlay
deposition_toggle = True # Only Show Primary Cacsades that include Deposition ELectrons
HPC_toggle = False
show_secondary_toggle = False
desktop_toggle = True

# IO parameters
date = "/July/7_7_2021/"
deposit_type = 'line'#types supported are line and point
material = 'tungsten'
cross_section_type = '_SM_'
parameter_summary = "20keV_1_1kpp_pitch_0_401_401_1001_sb_1000_sd_34000_sh_996_detect_dome_"

#Input Paths
if HPC_toggle: #Adjustments for HPC environment
    #HPC path
    path_to_test_files = "/home/richarddejong/nebula_test_files/"
elif desktop_toggle:
    path_to_test_files = "C:/Users/Richard/source/repos/Nebula/nebula_test_files/"
else: 
    # Laptop path
    path_to_test_files= "C:/Users/richa/Documents/repos/nebula_test_files/"

input_path = path_to_test_files+"output"+date+parameter_summary+material+cross_section_type+"cascade.bin"
data_path = path_to_test_files+"output"+date+parameter_summary+material+cross_section_type+"geometry_data.p"#Primary Tags 
geometry_data = pickle.load(open( data_path, "rb" ) )

#Geometry_data definition from read_geometry_bin #Debug tags are rule based 
""" 
geometry_data = [tag_set,heights,debug_tags] 
pickle.dump( geometry_data, open( filename_geometry_data, "wb" ) )
"""

deposition_tags = geometry_data[0]
height_map = geometry_data[1]
debug_tags = geometry_data[2]

#Cascade Choices
primary_select = 'random'
data_select = 'normal'
data_tags = [] # fill this with relevant primary tags if data_select = True
cascade_choices = [] 
voxel_size = 0.3

#Parameter Retrieval
path_list = (input_path.split('/')[-1]).split('.')
title =path_list[0]
title_list = title.split('_')
pillar_num = int(title_list[1])
pitch = int(title_list[4])
electron_num = pillar_num * int(title_list[2].replace('kpp',''))*1000
energy_str = title_list[0]
simulation_height = int(title_list[13])*voxel_size
z_height = (int(title_list[11])+int(title_list[7]))*voxel_size

output_path_folder = path_to_test_files+"figures"+date+title+"/"
general_path = output_path_folder+title
fignum = 1

if not os.path.exists(output_path_folder):
    os.makedirs(output_path_folder)

# General
line_scatter_str = general_path+"_np_"+str(primaries_to_plot)+"line_scatter"
line_scatter_side_str = general_path+"_np_"+str(primaries_to_plot)+"_line_scatter_zoom"

#2d
line_scatter_2d_str = general_path+"_np_"+str(primaries_to_plot)+"line_scatter_2_prim"
if show_secondary_toggle:
    line_scatter_2d_str = general_path+"_np_"+str(primaries_to_plot)+"line_scatter_2d_sec"


#Text Loading
tic =  time.perf_counter() 
print("Loading file...")

with open(input_path, mode='rb') as file: # b is important -> binary
   fileContent = file.read()

#Binary Unpacking ----------------------------------------------
lengrid = unpack("l", fileContent[:4])[0]
counter = 4

primaries= unpack("I"*lengrid,fileContent[counter:counter+lengrid*4])
counter+=lengrid*4
secondaries = unpack("I"*lengrid,fileContent[counter:counter+lengrid*4])
counter+=lengrid*4
x = unpack("f"*lengrid,fileContent[counter:counter+lengrid*4])
counter+=lengrid*4
y = unpack("f"*lengrid,fileContent[counter:counter+lengrid*4])
counter+=lengrid*4
z = unpack("f"*lengrid,fileContent[counter:counter+lengrid*4])
counter+=lengrid*4

toc = time.perf_counter()
print(f"Finished Loading in {toc - tic:0.4f} seconds")

#The stored bytes in order c++ ---------------------------------------
""" 
//Length 
std::ofstream output_bin_file(file_name, std::ios::binary);
int64_t len = cascade_diagram.primaries.size();
output_bin_file.write( (char*)&len, sizeof(len) );

//Vectors
output_bin_file.write( (char*)&cascade_diagram.primaries[0], len * sizeof(uint32_t) );
output_bin_file.write( (char*)&cascade_diagram.secondaries[0], len * sizeof(uint32_t) );
output_bin_file.write( (char*)&cascade_diagram.x_pos[0], len * sizeof(float_t) );
output_bin_file.write( (char*)&cascade_diagram.y_pos[0], len * sizeof(float_t) );
output_bin_file.write( (char*)&cascade_diagram.z_pos[0], len * sizeof(float_t) );
"""

#Data Manipulation --------------------------------------------------

#Timing
tic =  time.perf_counter() 
print("Sorting into Primary Cascade Groups...")

plt.rcParams['font.size'] = 16
zero_point_voxel = int(float(title_list[6])/2)
zero_point = zero_point_voxel*voxel_size
primaries_array = np.asarray(primaries)
secondaries_array =  np.asarray(secondaries)
secondaries_array[0] = 0
x_array = np.asarray(x)
y_array = np.asarray(y)
z_array = np.asarray(z)

#Index Finding
min_index = int(zero_point-f_range/voxel_size)
max_index = int(zero_point+f_range/voxel_size)

# Axis Shifting
x_array = zero_point-x_array
y_array = zero_point-y_array

#Further Manipulation for Deposition Primaries -------------------------

active_cascade_array = np.transpose(np.vstack((primaries_array,secondaries_array,x_array, y_array,z_array)))

if deposition_toggle: #Replaces Full Combination Array with Array Consisting of Only Cascades Resulting in Depositions(Also supports further specification)
    if data_select=='debug':
        dep_bool = (np.in1d(primaries_array, debug_tags)).astype(int)
    else:
        dep_bool = (np.in1d(primaries_array, deposition_tags)).astype(int)
    dep_indices = np.where(dep_bool==0)

    active_cascade_array = np.delete(active_cascade_array, dep_indices, axis=0) # all deposition information for deposition primaries

unique_primaries, event_counts = np.unique(active_cascade_array[:,0], return_counts=True)

# Retrieving Individual Cascades----------------------

#Rough Cascade Splitter -Simple Primary Families
total_events =  len(active_cascade_array)
primary_cascades = []
total_counter = 0
for i in range(len(unique_primaries)):
    cascade_data = []
    for j in range(event_counts[i]):
        cascade_data.append(active_cascade_array[total_counter+j,:])
    primary_cascades.append(cascade_data)
    total_counter+=event_counts[i]
    if total_events%100 == 0:
        print(str(total_counter/total_events*100)+"% ")
    cascade_data = []

toc = time.perf_counter()
print(f"Finished Sorting Primary Cascade Groups in {toc - tic:0.4f} seconds")

#Secondary Splitter
def secondary_splitter(primary_cascades):
    """
    This Method splits a primary cascade into seperate secondary cascades
    """
    cascade_collection = []
    cascade_array = np.vstack(primary_cascades)
    n_paths = int(np.max(cascade_array[:,1]))+1

    # Simple Split
    for i in range(n_paths):
        cascade_collection.append(cascade_array[np.where(cascade_array[:,1] == i)])
    return cascade_collection
    
#Line Plot 3D ------------------------------
fig = plt.figure(fignum)
fignum+=1
plt.xlabel('X(nm)')
plt.ylabel('Simulation Depth Z(nm)')

cascade_choices = [] # fill this with chosen primaries if needed

#Primary Selection

if primary_select == 'random':
    cascade_choices = []
    for i in range(primaries_to_plot):
        cascade_choices.append(random.choice(primary_cascades))

#Manual selection
if primary_select == 'manual':
    cascade_choices = [] 
    for primary_tag in cascade_choices:
        cascade_choices.append(primary_cascades[0])# fill in primary ids here

#Need to change this to allow for unpickling of previous primaries or ect.
if primary_select == 'load_previous':
    cascade_choices = []
    for i in range(primaries_to_plot):
        cascade_choices.append(random.choice(primary_cascades))

# Improved Path Diagram

print("Creating Plots...")
tic = time.perf_counter()

plt.rcParams['font.size'] = 16
if overlay_toggle:
    print("Generating Surface...")
    tic_surf = time.perf_counter()

    # Deposition Surface Overlay
    x_surf = np.linspace(-zero_point,zero_point,num =int(title_list[6]))
    y_surf = np.linspace(zero_point,-zero_point,num =int(title_list[6]))
    #z_surf = simulation_height-height_map[:,zero_point_voxel]
    height_map[height_map<0] = 0 
    flat_map = np.amax(height_map,axis=1)
    z_surf = simulation_height-flat_map

    plt.plot(x_surf,z_surf,color = 'k',alpha = 0.5)
    toc_surf = time.perf_counter()
    print(f"Finished Generating Surface {toc_surf - tic_surf:0.4f} seconds")

for primary_cascade in cascade_choices:
    secondary_cascades = secondary_splitter(primary_cascade)
    #Random Color Generation per Primary
    r = random.uniform(0, 0.5)
    b = random.uniform(0, 0.5)
    g = random.uniform(0, 0.5)
    primary_color = (r, g, b)
    secondary_color = (r+0.1,g+0.1,b+0.1)
    num_secondaries = len(secondary_cascades)
    color_it = 0.4/num_secondaries
    iteration_color = primary_color
    if show_secondary_toggle:
        for s_cascade in secondary_cascades:
            plt.plot(s_cascade[:,2], s_cascade[:,4],ms = 1,color = iteration_color,linewidth = 0.3)
            plt.scatter(s_cascade[-1,2], s_cascade[-1,4],s = 1,color = iteration_color,linewidth = 0.1) #terminations
            iteration_color = (iteration_color[0]+color_it,iteration_color[1]+color_it,iteration_color[2]+color_it)
    else:
        s_cascade = secondary_cascades[0]
        plt.plot(s_cascade[:,2], s_cascade[:,4],ms = 1,color = iteration_color,linewidth = 0.3)
        plt.scatter(s_cascade[-1,2], s_cascade[-1,4],s = 1,color = iteration_color,linewidth = 0.1) #terminations
        iteration_color = (iteration_color[0]+color_it,iteration_color[1]+color_it,iteration_color[2]+color_it)

#Horizontal Surface Line        
xs = np.linspace(-f_range,f_range,20)
horiz_line_data = np.array([simulation_height for i in xs])
plt.plot(xs, horiz_line_data, 'r',linewidth = 0.5) 

plt.xlim(-f_range,+f_range)
plt.ylim(round(simulation_height+z_roi_limits[0],-1),round(simulation_height+z_roi_limits[1],-1))   
#plt.ylim(-z_roi_limit,+z_roi_limit)   
plt.gca().invert_yaxis()
if show_secondary_toggle:
    title_scat = str(primaries_to_plot)+" Deposition Relevant Primary and Secondary Scattering Cascades for " + energy_str
else:
    title_scat = str(primaries_to_plot)+" Deposition Relevant Primary Scattering Cascades for " + energy_str
#plt.title(title_scat)
fig.savefig(line_scatter_2d_str, dpi = 400, bbox_inches="tight")



toc = time.perf_counter()
print(f"Finished  Plots in {toc - tic:0.4f} seconds")
