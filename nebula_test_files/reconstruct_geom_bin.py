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
import matplotlib
import time
from struct import *

#This file is for reconstructing deposits into format for reinsertion into simulator(for simulated SEM images for instance)


# IO parameters
date = "/August/29_8_2021/"
material = 'tungsten'
deposit_type = "point"
#material = 'graphite-phonon'
cross_section_type = '_SM_'
direct_pathing = False #allow for direct input of file through command line
HPC_toggle = False
desktop_toggle = True
SEM_toggle = True #creates geometry suited to SEM (has code for simple elevation atm) with mirrors
#SEM_routine = "_elevation" # "_elevation" or "_mirrors" are only available options (_mirrors is just id for SEM_Toggle = true)
#SEM_routine = "_mirrors"
SEM_routine = "_mirrors"


parameter_summary = "1keV_1_5kpp_pitch_0_161_161_1001_sb_1000_vs_250_sd_34000_sh_996_detect_dome_mirror_vs_250_" #point
#parameter_summary = "1keV_21_1kpp_pitch_2400_seq_p_401_401_1001_sb_1000_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_" #line
#parameter_summary = "1keV_21_1kpp_pitch_3000_seq_l_401_401_1001_sb_1000_sd_34000_sh_996_detect_dome_mirror_vs_300_" #line
#parameter_summary = "1keV_1_5kpp_pitch_0_81_81_101_sb_1000_sd_34000_sh_96_detect_dome_vs_250_" # pillar
#parameter_summary = "1keV_11_21_1kpp_pitchx_6000_pitchy_1700_seq_d_401_401_1001_sb_1000_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_" #wall
#parameter_summary = "1keV_21_1kpp_pitchy_1700_pitchx_5500_seq_d_401_401_1001_sb_1000_sd_1000_sh_996_detect_dome_" #wall
#parameter_summary = "1keV_rpit_1700_ppit_4000_1kpp_425p_pass_1_radii_30_4_cone_dep_random_401_401_1001_sb_1000_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_"

index = parameter_summary.find('kpp')
index_l = 3
if deposit_type == 'line':
    index = parameter_summary.find('seq')
    index_l = 5
if deposit_type == 'wall':
    index = parameter_summary.find('seq')
    index_l = 5
if deposit_type == 'cone':
    index = parameter_summary.find('401')
    index_l = -1
SEM_routine_extra_str = "_RC_"+parameter_summary[:index+index_l]
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

input_path = path_to_test_files+"output"+date+parameter_summary+material+cross_section_type+"output.bin"

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
electrons_per_pillar = int(title_list[2].replace('kpp',''))*1000

output_path_folder = path_to_test_files+"figures"+date+title
#output_path_folder = "/mnt/c/Users/richa/Documents/repos/nebula_test_files/figures"+date+title #WSL
  
general_path = output_path_folder+"/"+title
fignum = 1

#Text Loading
print("Loading file...")

tic =  time.perf_counter() 
###Change this to allow for bin outputs

with open(input_path, mode='rb') as file: # b is important -> binary
   fileContent = file.read()

#Binary Checks
unpack_flag = sys.byteorder
print('matplotlib: {}'.format(matplotlib.__version__))
#Binary Unpacking
print("System Byte Order: "+ unpack_flag)

lengrid = unpack("q", fileContent[:8])[0]
voxel_size = unpack("f", fileContent[8:12])[0]
dim = unpack("iii", fileContent[12:24])
counter = 24

mat_grid_in = unpack("h"*lengrid,fileContent[counter:counter+lengrid*2])
counter+=lengrid*2
tag_grid_in = unpack("i"*lengrid,fileContent[counter:counter+lengrid*4])
counter+=lengrid*4
electron_type_in = unpack("h"*lengrid,fileContent[counter:counter+lengrid*2])
counter+=lengrid*2
new_species_in = unpack("h"*lengrid,fileContent[counter:])

toc = time.perf_counter()
print(f"Finished Loading in {toc - tic:0.4f} seconds")


#Start Timing Data Manipulation
tic = time.perf_counter()
print("Data Manipulation...")

#Grid Work (Only need mat_grid)
mat_grid = np.reshape(mat_grid_in, dim, order='F')/123 + 1
mat_grid_og = np.reshape(mat_grid_in, dim, order='F')

mat_grid_t = mat_grid.transpose(1, 0, 2)

toc = time.perf_counter()
print(f"Finished Data Manipulation in {toc - tic:0.4f} seconds")

tic = time.perf_counter()
print("Reconstructing Geometry around Slice...")


# Target Geometry
voxel_size = 0.25; # voxel size in nanometers (0.27 nm is appr. one atom of Si)
size_x = 161 # horizontal size in the x direction in voxels (now for +/- x)
size_y = 161 # horizontal size in the y direction in voxels (now for +/- y)
size_z = 1001# vertical size in voxels
volume = size_x*size_y*size_z # total voxel volume
sim_depth = 34000 # simulation depth under the voxels at z < 0 for SEM bulk samples in voxels
sample_height = size_z-5 # 5 voxel ~ 1.5nm distance

#params to write
param_nums = [voxel_size,size_x,size_y,size_z,sample_height,sim_depth]
param_strs = [str(param) for param in param_nums]

# Geometry Construction

#Base Routine
ini_geom_3d = np.ones([size_x,size_y,size_z],dtype=np.int16)*np.int16(-123)
ini_geom_3d[:,:,sample_height:] = np.int16(0)

ini_geom_3d[:,:,0] = np.int16(-126) #creates top detector layer

#Reconstruction
slice_thickness = dim[2]
ini_geom_3d[:,:,sample_height-slice_thickness:sample_height] = mat_grid_og

#SEM Routines
#Mirrors
if SEM_toggle:
    ini_geom_3d[0,:,:] = np.int16(-122)
    ini_geom_3d[size_x-1,:,:] = np.int16(-122)
    ini_geom_3d[:,0,:] = np.int16(-122)
    ini_geom_3d[:,size_y-1,:] = np.int16(-122)

#Detect Dome
if SEM_routine == "_detect_dome":
    ini_geom_3d[0,:,:sample_height] = np.int16(-126)
    ini_geom_3d[size_x-1,:,:sample_height] = np.int16(-126)
    ini_geom_3d[:,0,:sample_height] = np.int16(-126)
    ini_geom_3d[:,size_y-1,:sample_height] = np.int16(-126)

#Elevation
if SEM_routine == "_elevation":
    center_vox = int(size_x/2)
    spread = 80
    elevation_height = 100
    vox_x_min = center_vox-spread
    vox_x_max = center_vox+spread
    mat_choice = 1
    ini_geom_3d[vox_x_min:vox_x_max,:,sample_height-elevation_height:sample_height] = np.int16(mat_choice)
    SEM_routine_extra_str = "_"+str(spread)+"_"+str(elevation_height)+"_"+str(mat_choice)

#Top Detector
ini_geom_3d[:,:,0] = np.int16(-126) #creates top detector layer

#Ready for binary write
flat_3d = ini_geom_3d.flatten('F')
len_ini_vec = np.int64(len(flat_3d)) #

#Timing Finish
toc =  time.perf_counter() #timing
print(f"Reconstructed Geometry in {toc - tic:0.4f} seconds")


tic = time.perf_counter()
print("Saving Geometry...")

# This is a numpy datatype that corresponds to io_structure for binary input files
"""
dt = np.dtype([
	("len_ini_vec",int),('voxel_size', float), ('size_x', int), ('size_y', float), ('size_z',int),    # Size
	('sim_depth', float), ('sample_height', int)]) #Depth   
"""
param_nums = [len_ini_vec,voxel_size,size_x,size_y,size_z,sample_height,sim_depth]
dt = np.dtype([
	("len_ini_vec",np.int64),('voxel_size', np.float32), ('size_x', np.int32), ('size_y', np.int32), ('size_z',np.int32),    # Size
	('sim_depth', np.float32), ('sample_height', np.int32)]) #Depth   

#Input Paths
if HPC_toggle: #Adjustments for HPC environment
    #HPC path
    output_path = "/home/richarddejong/nebula_test_files/vox_tri_pri/"
elif desktop_toggle:
    #Desktop path
    output_path = "C:/Users/Richard/source/repos/Nebula/nebula_test_files/vox_tri_pri/"
else: 
    # Laptop path
    output_path = "C:/Users/richa/Documents/repos/nebula_test_files/vox_tri_pri/"

# Naming for geometries
voxel_size_pm = int(voxel_size*1000)
name_str = str(size_x)+'_'+str(size_y)+'_'+str(size_z)+'_sd_'+str(sim_depth)+"_sh_"+str(sample_height)+SEM_routine+"_vs_"+str(voxel_size_pm)+SEM_routine_extra_str+'.bin'

with open(output_path+name_str, 'wb') as file:
    array = np.empty(1,dtype=dt)

    # Write buffer to file
    array['len_ini_vec'] = len_ini_vec
    array['voxel_size'] = voxel_size
    array['size_x'] = size_x
    array['size_y'] = size_y
    array['size_z'] = size_z
    array['sim_depth'] = sim_depth
    array['sample_height'] = sample_height

    array.tofile(file)
    flat_3d.tofile(file)

#Timing Finish
toc =  time.perf_counter() #timing
print(f"Saved Reconstructed Geometry in {toc - tic:0.4f} seconds")

