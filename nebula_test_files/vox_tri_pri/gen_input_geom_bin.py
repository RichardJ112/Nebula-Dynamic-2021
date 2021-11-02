import numpy as np
import sys
import pandas as pd
import struct
import time

"""
IO file generation

file structure:
1st line: voxel_size size_x size_y size_z sample_height sim_depth
	 * other lines: material 

saves materials in c order (i+j*y+k*x*y)
"""
print(sys.version)

HPC_toggle = False
desktop_toggle = True
SEM_toggle = False #creates geometry suited to SEM (has code for simple elevation atm)
#SEM_routine = "_elevation" # "_elevation" or "_mirrors" are only available options (mirrors is just id for SEM_Toggle = true)
SEM_routine = "_detect_dome"

SEM_routine_extra_str = ""

tic =  time.perf_counter() #timing

voxel_size = 0.3; # voxel size in nanometers (0.27 nm is appr. one atom of Si)
voxel_size_pm = int(voxel_size*1000)

size_x = 401 # horizontal size in the x direction in voxels (now for +/- x)
size_y = 401 # horizontal size in the y direction in voxels (now for +/- y)
size_z = 1001# vertical size in voxels
volume = size_x*size_y*size_z # total voxel volume
sim_depth = 34000 # simulation depth under the voxels at z < 0 for SEM bulk samples in voxels
sample_height = size_z-5 # 5 voxel ~ 1.5nm distance

#params to write
param_nums = [voxel_size,size_x,size_y,size_z,sample_height,sim_depth]
param_strs = [str(param) for param in param_nums]

# Geometry Construction

#Base Routine
ini_geom_3d = np.ones([size_x,size_y,size_z],dtype=np.int16)*np.int16(-123) # Everything starts as vacuum
ini_geom_3d[:,:,sample_height:] = np.int16(0) #sample height and below is material

#SEM Routines

#Mirrors
if SEM_toggle:
    ini_geom_3d[0,:,:] = np.int16(-122)
    ini_geom_3d[size_x-1,:,:] = np.int16(-122)
    ini_geom_3d[:,0,:] = np.int16(-122)
    ini_geom_3d[:,size_y-1,:] = np.int16(-122)

#Detect Dome
if SEM_routine == "_detect_dome":
    ini_geom_3d[0,:,:sample_height-1] = np.int16(-126)
    ini_geom_3d[size_x-1,:,:sample_height-1] = np.int16(-126)
    ini_geom_3d[:,0,:sample_height-1] = np.int16(-126)
    ini_geom_3d[:,size_y-1,:sample_height-1] = np.int16(-126)

#Detect Dome with Mirror Along seam
if SEM_routine == "_detect_dome_mirror":
    #Detection over area
    ini_geom_3d[0,:,:sample_height] = np.int16(-126)
    ini_geom_3d[size_x-1,:,:sample_height] = np.int16(-126)
    ini_geom_3d[:,0,:sample_height] = np.int16(-126)
    ini_geom_3d[:,size_y-1,:sample_height] = np.int16(-126)

    #Mirrors along seam
    ini_geom_3d[0,:,sample_height:] = np.int16(-122)
    ini_geom_3d[size_x-1,:,sample_height:] = np.int16(-122)
    ini_geom_3d[:,0,sample_height:] = np.int16(-122)
    ini_geom_3d[:,size_y-1,sample_height:] = np.int16(-122)

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
    output_path = "C:/Users/Richard/source/repos/Nebula/nebula_dynamic_2021/nebula_test_files/vox_tri_pri/"
else: 
    # Laptop path
    output_path = "C:/Users/richa/Documents/repos/nebula_test_files/vox_tri_pri/"

# Naming for geometries
name_str = str(size_x)+'_'+str(size_y)+'_'+str(size_z)+'_sd_'+str(sim_depth)+"_sh_"+str(sample_height)+SEM_routine+SEM_routine_extra_str+"_vs_"+str(voxel_size_pm)+'.bin'

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
toc = time.perf_counter()
print(f"Finished Generating in {toc - tic:0.4f} seconds")


