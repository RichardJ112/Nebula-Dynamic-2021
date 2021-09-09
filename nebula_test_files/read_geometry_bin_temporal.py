import numpy as np
import os
import pickle
import sys
from statistics import median 
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import rc
from matplotlib import cm
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import time
from struct import *

#This file is for reading line and pillar binary deposition geometries but only creates plots for temporal distributions for electron classification
# and dimensional analysis, many things are stripped out to make execution quicker 
# See read_geometry_bin for general analysis

# IO parameters
#date = "/July/30_7_2021/"
date = "/August/30_8_2021/"
deposit_type = 'point'#types supported are line and point
material = 'tungsten'
#material = 'graphite-phonon'
#cross_section_type = '_SMC_'
cross_section_type = '_SM_'
HPC_toggle = False #Set to true when running on HPC
direct_pathing = False #allow for direct input of file through command line
legacy_class_toggle = True # turn off old classifcation plots
desktop_toggle = True
point_contour = True
enable_e_grid = False
HELIOS_data = False
nm_check = 2
memory_saver = True
max_toggle = True #if true takes max height whole substrate
electron_bin = 10
helios_s_run = 10e3 #amount of electrons until restriction 90e3 --GH, 380e3 --NGH


gas_handling_str = "_e20_em0_m"
#gas_handling_str = "_e20_em9"
#gas_handling_str = ""

#output_extra_str = ""
output_extra_str = ""

#parameter_summary = "1keV_1_10kpp_pitch_0_161_161_1001_sb_1000_vs_357_sd_34000_sh_996_detect_dome_mirror_vs_357_" #point
#parameter_summary = "1keV_1_600kpp_pitch_0_251_251_1001_sb_3300_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_" #point
#parameter_summary = "1keV_1_400kpp_pitch_0_161_161_1001_sb_1000_vs_357_sd_34000_sh_996_detect_dome_mirror_vs_357_"
parameter_summary = "1keV_1_10kpp_pitch_0_401_401_1001_sb_1000_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_"
#parameter_summary = "1keV_1_500kpp_pitch_0_161_161_1001_sb_1000_sd_34000_sh_996_detect_dome_mirror_vs_357_"
#parameter_summary = "1keV_1_600kpp_pitch_0_251_251_1001_sb_3300_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_"
#parameter_summary = "20keV_1_50kpp_pitch_0_401_401_1001_sb_1000_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_" #point

#Input Paths
if HPC_toggle: #Adjustments for HPC environment
    #HPC path
    path_to_test_files = "/home/richarddejong/nebula_test_files/"
elif desktop_toggle:
    path_to_test_files = "C:/Users/Richard/source/repos/Nebula/nebula_test_files/"
else: 
    # Laptop path
    path_to_test_files= "C:/Users/richa/Documents/repos/nebula_test_files/"

#Direct Pathing via system argument (overwrites input_path above(does not change IO parameters))

input_path = path_to_test_files+"output"+date+parameter_summary+material+cross_section_type+"output"+gas_handling_str+".bin"

#Arg Tests
print(sys.argv[0])
if direct_pathing:
    if(len(sys.argv) > 0):
	    input_path = str(sys.argv[-1])
    #Direct Pathing Translator (Just changes the \\ to / for interpretation)
    #input_path = input_path.replace("\\","/")

#Path Interpreter
path_list = (input_path.split('/')[-1]).split('.')
title =path_list[0]
title_list = title.split('_')
pillar_num = int(title_list[1])
pitch = float(int(title_list[4])/1000)
electrons_per_pillar = int(title_list[2].replace('kpp',''))*1000
electron_num = pillar_num * electrons_per_pillar
energy_str = title_list[0]

#Plotting Slices
#electron_bin = int(50e3)

#max_t = 10000

max_t = electron_num


output_path_folder = path_to_test_files+"figures"+date+title
#output_path_folder = "/mnt/c/Users/richa/Documents/repos/nebula_test_files/figures"+date+title #WSL
  
general_path = output_path_folder+"/"+title
fignum = 1

if not os.path.exists(output_path_folder):
    os.makedirs(output_path_folder)

# General path strings
T_WH_str = general_path+"_WH_Temporal" # temporal width/height figure
T_EC_str = general_path+"_EC_Temporal" #temporal electron classification
T_WH_ratio_str = general_path +"_WH_Ratio_Temporal"

#Binary Checks
unpack_flag = sys.byteorder
print('matplotlib: {}'.format(matplotlib.__version__))
#Binary Unpacking
print("System Byte Order: "+ unpack_flag)

#Text Loading
print("Loading file...")

tic =  time.perf_counter() 
###Change this to allow for bin outputs

with open(input_path, mode='rb') as file: # b is important -> binary
   fileContent = file.read()

#Binary Unpacking -HPC Changes made
#lengrid = unpack("l", fileContent[:4])[0]
lengrid = unpack("q", fileContent[:8])[0]
voxel_size = unpack("f", fileContent[8:12])[0]
dim = unpack("iii", fileContent[12:24])
counter = 24

mat_grid_in = unpack("h"*lengrid,fileContent[counter:counter+lengrid*2])
counter+=lengrid*2
tag_grid_in = unpack("i"*lengrid,fileContent[counter:counter+lengrid*4])
counter+=lengrid*4
if enable_e_grid:
    e_grid_in = unpack("f"*lengrid,fileContent[counter:counter+lengrid*4])
    counter+=lengrid*4
electron_type_in = unpack("h"*lengrid,fileContent[counter:counter+lengrid*2])
counter+=lengrid*2
new_species_in = unpack("h"*lengrid,fileContent[counter:])


#The stored bytes in order c++
""" 
std::ofstream output_bin_file(file_name, std::ios::binary);
	int64_t len = _mat_grid.size();
	output_bin_file.write( (char*)&len, sizeof(len) );
	output_bin_file.write( (char*)&_voxel_size, sizeof(real) );
	output_bin_file.write( (char*)&_size_x, sizeof(int32_t) );
	output_bin_file.write( (char*)&_size_y, sizeof(int32_t) );
	output_bin_file.write( (char*)&_size_z, sizeof(int32_t) );

	//Vectors
	output_bin_file.write( (char*)&_mat_grid[0], len * sizeof(int16_t) );
	output_bin_file.write( (char*)&_tag_grid[0], len * sizeof(int32_t) );
	output_bin_file.write( (char*)&_e_grid[0], len * sizeof(real) );
	output_bin_file.write( (char*)&_species_grid[0], len * sizeof(int16_t) );
    output_bin_file.write( (char*)&_new_species_grid_slice_int16[0], len * sizeof(int16_t) );

    
    //output_bin_file.write((char*)&_new_species_grid_flat_slice[0], len*8*sizeof(char);) //new bin files for classification
"""

toc = time.perf_counter()
print(f"Finished Loading in {toc - tic:0.4f} seconds")


#Start Timing Data Manipulation
print("Data Manipulation...")

tic = time.perf_counter()

#Grid Work
tag_grid = np.reshape(tag_grid_in, dim, order='F')   
new_species_grid = np.reshape(new_species_in, dim, order='F')

#Fix Plots --switching x and y axes (Comes from z axis flip)

tag_grid_t = tag_grid.transpose(1, 0, 2)
new_species_grid_t = new_species_grid.transpose(1,0,2)

toc = time.perf_counter()
print(f"Finished Data Manipulation in {toc - tic:0.4f} seconds")

def height_map(tag_grid,axis,invalid_val = dim[2]):
    
    """ 
    from tag_grid create 2d numpy map with maximum height values on each column

    Parameters
    ----------
    tag_grid_f : primary tag for each deposited electron

    Returns
    -------
    height_map : 2d numpy array of maximum values

    """
    mask = tag_grid!=0
    height_indices = np.where(mask.any(axis=axis), mask.argmax(axis=axis), invalid_val)
    height_map = (dim[2]-height_indices)*voxel_size
    return height_map

def temporal_grids(species_grid,tag_grid,tag_limits):
    """ 
    Splits species_grid into seperate temporal grids for dimensional analysis and classificiation

    Parameters
    ----------
    species_grid : np.ndarray (3d) - consisting of electron classifications
    tag_limits: np.array (1d) - consists of tag_limits per temporal slice
    tag_grid : primary tag for each deposited electron

    Returns
    -------
    temp_grids : list of np.ndarray (3d) -returns grids

    """
    temp_grids =[]
    for tag_limit in tag_limits:
        grid_i = np.where(tag_grid < tag_limit, species_grid, 0)
        temp_grids.append(grid_i)
    return temp_grids

def type_conc(species_grid,conc_type):
    """ 
    Smith et Al Amount Concentrations (Using Smith definitions of PE,FSE,BSE,SE1,SE2)

    Parameters
    ----------
    species_grid : np.ndarray - 3 dimensional array consisting of electron classifications
    conc_type : choose electron type model

    Returns
    -------
    type_counts : 5 integers indicating the concentration of each electron type caused deposition

    """
    if conc_type == "smith":
        S_T = np.sum(species_grid != 0)
        S_PE = np.sum(species_grid == 1)/S_T
        S_FSE = np.sum((species_grid == 2) | (species_grid == 5))/S_T
        S_BSE = np.sum((species_grid == 3) | (species_grid == 4))/S_T
        S_SE_1 = np.sum((species_grid == 8) | (species_grid == 9))/S_T #secondary electrons from intial traversal
        S_SE_2 = np.sum((species_grid == 10) | (species_grid == 11))/S_T # secondary electrons due to BSE
        type_conc = [S_PE,S_FSE,S_BSE,S_SE_1,S_SE_2]
        return type_conc

def find_lower_bound_index(pillar_vals,zero_index):
    """
    Function that finds first non-sequential index in list of indices and pillar values
    """
    zero_indices = zero_index[0]
    for i in range(len(zero_indices)-1):
        if zero_indices[i+1] != zero_indices[i]+1:
            if pillar_vals[zero_indices[i+1]+1] != 0:
                lsi = int(zero_indices[i+1])
                return lsi
    return 0

tag_limits = np.arange(electron_bin,max_t+1,electron_bin)

e_counts = []
pillar_heights = []
pillar_widths_nm_check =[]
pillar_HDWs = []
midpoint = int(dim[1]/2)

if memory_saver:
    tic = time.perf_counter()
    print("Processing Slices...")
    for tag_limit in tag_limits:
        species_grid = np.where(tag_grid < tag_limit, new_species_grid, 0)
        #Dimensions
        e_counts.append(type_conc(species_grid,"smith"))
        vert_slice_x =  np.transpose(species_grid[ : , midpoint, ::-1]) #takes vertical_slice at pillar position along x = 0
        i = 1
        j = 0
        pillar_vals = vert_slice_x[:,midpoint]# takes single column at x = pillar_pos , y = 0 
        zero_index = np.nonzero(pillar_vals)
        if zero_index[0].size == 0:
            pillar_heights.append(0)
            pillar_widths_nm_check.append(0)
            pillar_HDWs.append(0)
            continue
        height_index = np.max(zero_index) 
        pillar_height_i = height_index*voxel_size
        if max_toggle:
            h_map = height_map(species_grid,axis=2)
            pillar_heights.append(np.max(h_map))
        else:
            pillar_heights.append(pillar_height_i)
        #lower_bound_index = find_lower_bound_index(pillar_vals,zero_index) <-- we use 0 as this is pillar
        lower_bound_index = 0
        half_diameter_index = int((height_index-lower_bound_index)/2 + lower_bound_index) #index of half diameter
        row_at_nm_check = vert_slice_x[int(nm_check/voxel_size)]
        row_half = vert_slice_x[half_diameter_index]
        zero_index_row_nm_check = np.nonzero(row_at_nm_check)
        zero_index_row_half =  np.nonzero(row_half)
        if pillar_height_i/2 > 1:
            pillar_HDWs.append((np.max(zero_index_row_half)-np.min(zero_index_row_half))*voxel_size)
        else:
            pillar_HDWs.append(0)
        if zero_index_row_nm_check[0].size == 0:
            pillar_widths_nm_check.append(0)
            continue
        pillar_widths_nm_check.append((np.max(zero_index_row_nm_check)-np.min(zero_index_row_nm_check))*voxel_size)

    toc = time.perf_counter()
    print(f"Finished Processing in {toc - tic:0.4f} seconds")
else:
    tic = time.perf_counter()
    print("Creating Temporal Slices...")

    t_grids = temporal_grids(new_species_grid_t,tag_grid_t,tag_limits)

    toc = time.perf_counter()
    print(f"Finished Temporal Slices in {toc - tic:0.4f} seconds")

    tic = time.perf_counter()
    print("Procesing Slices...")

    #Temporal Study Pillars
    for species_grid in t_grids:
        #Dimensions
        e_counts.append(type_conc(species_grid,"smith"))
        vert_slice_x =  np.transpose(species_grid[ : , midpoint, ::-1]) #takes vertical_slice at pillar position along x = 0
        i = 1
        j = 0
        pillar_vals = vert_slice_x[:,midpoint]# takes single column at x = pillar_pos , y = 0 
        zero_index = np.nonzero(pillar_vals)
        if zero_index[0].size == 0:
            pillar_heights.append(0)
            pillar_widths_nm_check.append(0)
            pillar_HDWs.append(0)
            continue
        height_index = np.max(zero_index) 
        pillar_height_i = height_index*voxel_size
        
        if max_toggle:
            h_map = height_map(species_grid,axis=2)
            pillar_heights.append(np.max(h_map))
        else:
            pillar_heights.append(pillar_height_i)
        #lower_bound_index = find_lower_bound_index(pillar_vals,zero_index) <-- we use 0 as this is pillar
        lower_bound_index = 0
        half_diameter_index = int((height_index-lower_bound_index)/2 + lower_bound_index) #index of half diameter
        row_at_nm_check = vert_slice_x[int(nm_check/voxel_size)]
        row_half = vert_slice_x[half_diameter_index]
        zero_index_row_nm_check = np.nonzero(row_at_nm_check)
        zero_index_row_half =  np.nonzero(row_half)
        if pillar_height_i/2 > 1:
            pillar_HDWs.append((np.max(zero_index_row_half)-np.min(zero_index_row_half))*voxel_size)
        else:
            pillar_HDWs.append(0)
        if zero_index_row_nm_check[0].size == 0:
            pillar_widths_nm_check.append(0)
            continue
        pillar_widths_nm_check.append((np.max(zero_index_row_nm_check)-np.min(zero_index_row_nm_check))*voxel_size)

    toc = time.perf_counter()
    print(f"Finished  Processing in {toc - tic:0.4f} seconds")

# Parameter Study

print(title)
tic = time.perf_counter()
print("Creating Plots...")

# Dimensions Over Time
T_WH_plot = plt.figure(fignum)
fignum+=1

labels_distro = ["Pillar Height","Pillar Width at "+str(nm_check)+"nm","Pillar FWHM"]

plt.rcParams['font.size'] = 16

plt.plot(tag_limits,pillar_heights)
plt.plot(tag_limits,pillar_widths_nm_check)
plt.plot(tag_limits,pillar_HDWs)

#Dimensions over time extra -HELIOS data
if HELIOS_data:
    helios_tags = helios_s_run*np.linspace(1,6,6)
    helios_tags[2] = 2.5*helios_s_run
    helios_heights = np.array([39.8,52.2,61.1,73.6,89.3,100.9])
    helios_widths = np.array([40,53.4,55.4,62.6,65.3,64.9])
    plt.plot(helios_tags,helios_heights)
    plt.plot(helios_tags,helios_widths)
    plt.xticks(fontsize = 12)
    labels_distro.append("HELIOS Heights")
    labels_distro.append("HELIOS Widths")


plt.legend(labels_distro,loc='lower right')
plt.xlabel(" Dwell Time (# of electrons)")
plt.ylabel(" Dimension (nm)")
#plt.title("Temporal Height and Width of "+energy_str +" Beam")
ymax = max(max(np.max(pillar_heights),np.max(pillar_widths_nm_check)),np.max(pillar_HDWs))
plt.ylim(0,ymax)
plt.xlim(0,max_t)
T_WH_plot.savefig(T_WH_str, dpi = 200, bbox_inches="tight")

# Distributions Over Time
T_EC_plot = plt.figure(fignum)
fignum+=1

labels_distro = ["PE","FSE","BSE","SE I","SE II"]

plt.rcParams['font.size'] = 16
e_class_array = np.asarray(e_counts)
plt.plot(tag_limits,e_class_array[:,0], color = 'red' ) #PE
plt.plot(tag_limits,e_class_array[:,1], color = 'green' ) #FSE
plt.plot(tag_limits,e_class_array[:,2], color = 'blue' ) #BSE
plt.plot(tag_limits,e_class_array[:,3], color = '#ffa020' ) #SE_1
plt.plot(tag_limits,e_class_array[:,4], color = 'magenta' ) #SE_2

#plt.legend(labels_distro,loc='upper left')
plt.legend(labels_distro,loc='lower right')
plt.xlabel(" Dwell Time (# of electrons)")
plt.ylabel(" Deposition Type Concentration (-)")
#plt.title("Temporal Deposition Electron Concentration of "+energy_str +" Beam")
#max_t_alt= 4000
max_t_alt = max_t
plt.xlim(0,max_t_alt)
plt.ylim(0,1)
T_EC_plot.savefig(T_EC_str, dpi = 200, bbox_inches="tight")

if HELIOS_data:
    T_WH_ratio_plot = plt.figure(fignum)
    fignum+=1
    plt.xticks(fontsize = 12)
    labels_distro = ["NBD21 2nm","NBD21 FWHM","HELIOS"]
    pillar_heights = np.asarray(pillar_heights)
    pillar_widths_nm_check = np.asarray(pillar_widths_nm_check)
    pillar_HDWs = np.asarray(pillar_HDWs)
    plt.plot(tag_limits,np.divide(pillar_heights,pillar_widths_nm_check))
    plt.plot(tag_limits,np.divide(pillar_heights,pillar_HDWs))
    plt.plot(helios_tags,np.divide(helios_heights,helios_widths))
    plt.legend(labels_distro,loc='upper left')
    plt.xlabel(" Dwell Time (# of electrons)")
    plt.ylabel(" Height to Width Ratio [-]")
    plt.xlim(0,max_t)
    
    #plt.xticks(np.arange(np.min(tag_limits), np.max(tag_limits)+1, 1.0))
    T_WH_ratio_plot.savefig(T_WH_ratio_str, dpi = 200, bbox_inches="tight")

toc = time.perf_counter()
print(f"Finished  Plots in {toc - tic:0.4f} seconds")
