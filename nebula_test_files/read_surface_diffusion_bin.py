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
from collections import OrderedDict
import time
import math
from struct import *
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

#This file is for reading surface voxels

#IDs
"""
0 - no ID
1 - surface voxels
2 - open site
3 - occupied site
4 - sources

"""


# IO parameters
date = "/August/30_8_2021/"
deposit_type = 'point'#types supported are line and point
material = 'tungsten'
cross_section_type = '_SM_'
HPC_toggle = False #Set to true when running on HPC
direct_pathing = False #allow for direct input of file through command line
legacy_class_toggle = False # turn off old classifcation plots
desktop_toggle = True
distance_toggle = False
tracking = False
#source = 'walls' #walls or point
source = 'walls'
extra_plots = True
f_range = 20#fixed range in nm

#output_extra_str = ""
output_extra_str = "_e17_em10_mb"


parameter_summary = "1keV_1_400kpp_pitch_0_161_161_1001_sb_1000_vs_357_sd_34000_sh_996_detect_dome_mirror_vs_357_"
#parameter_summary = "5keV_21_10kpp_pitch_1000_401_401_601_sb_1000_sd_10000_"
#parameter_summary = "1keV_1_0-31kpp_pitch_0_57_57_101_sb_1000_vs_357_sd_34000_sh_96_detect_dome_vs_357_" # 2D diffusion wall
#parameter_summary = "1keV_1_0-193kpp_pitch_0_113_113_101_sb_1000_vs_357_sd_34000_sh_96_detect_dome_vs_357_" # 2D diffusion point
#parameter_summary = "1keV_1_5kpp_pitch_0_161_161_1001_sb_1000_vs_250_sd_34000_sh_996_detect_dome_mirror_vs_250_" #small pillar
#parameter_summary = "1keV_1_2kpp_pitch_0_161_161_1001_sb_1000_vs_250_sd_34000_sh_996_mirrors_vs_250_RC_1keV_1_5kpp_" #small pillar reconstructed
#parameter_summary = "1keV_1_0-100kpp_pitch_0_85_85_101_sb_1000_sd_34000_sh_96_detect_dome_vs_960_"
#parameter_summary = "1keV_1_10kpp_pitch_0_81_81_101_sb_1000_sd_34000_sh_96_detect_dome_vs_250_RC_1keV_1_5kpp_"
#parameter_summary = "1keV_1_5kpp_pitch_0_81_81_101_sb_1000_sd_34000_sh_96_detect_dome_vs_250_"
#parameter_summary = "1keV_1_250kpp_pitch_0_161_161_1001_sb_1000_sd_34000_sh_996_detect_dome_mirror_vs_357_"
#parameter_summary = "1keV_1_0-1kpp_pitch_0_401_401_101_sb_1000_sd_34000_sh_96_detect_dome_vs_250_"

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

input_path = path_to_test_files+"output"+date+parameter_summary+material+cross_section_type+"surface"+output_extra_str+".bin"

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
if title_list[2][0] == '0':
	electrons_per_pillar = int(title_list[2].split('-')[1].replace('kpp','')) #<1000 electrons
else:
	electrons_per_pillar = int(title_list[2].replace('kpp',''))*1000 #1000s of electrons
electron_num = pillar_num * electrons_per_pillar

output_path_folder = path_to_test_files+"figures"+date+title
#output_path_folder = "/mnt/c/Users/richa/Documents/repos/nebula_test_files/figures"+date+title #WSL
  
general_path = output_path_folder+"/"+title
fignum = 1

if not os.path.exists(output_path_folder):
    os.makedirs(output_path_folder)

# General path strings
cs_mat_vert_x0_str = general_path +"_cs_vert_mat_X0" # vertical cross section of materials at X=0
cs_surf_vert_x0_str = general_path +"_cs_vert_surf_X0" # vertical cross section of surface identifiers at X=0

cs_mat_vert_mx_str = general_path +"_cs_vert_mat_Xm" # vertical cross section at a specific -X
cs_surf_vert_mx_str = general_path +"_cs_vert_surf_Xm" # vertical cross section at a specific -X
cs_mat_vert_px_str = general_path +"_cs_vert_mat_Xp" # vertical cross section at a specific +X
cs_surf_vert_px_str = general_path +"_cs_vert_surf_Xp" # vertical cross section at a specific +X

cs_surf_hor_rms = general_path + "_cs_hor_surf_rms" #hotizontal cross-section with rms
cs_surf_hist = general_path +"_cs_hor_surf_hist" # vertical cross section at a specific -X
cs_surf_hor_adsorb = general_path +"_cs_hor_surf_adsorb" # vertical cross section at a specific -X

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

#Binary Unpacking -HPC Changes made
lengrid = unpack("q", fileContent[:8])[0]
voxel_size = unpack("f", fileContent[8:12])[0]
dim = unpack("iii", fileContent[12:24])
counter = 24


mat_grid_in = unpack("h"*lengrid,fileContent[counter:counter+lengrid*2])
counter+=lengrid*2
surface_grid_in = unpack("h"*lengrid,fileContent[counter:counter+lengrid*2])
counter+=lengrid*2
lendist = unpack("q", fileContent[counter:counter+8])[0] #i[4 bits] normally q[8 bits] after 30/8/2021
counter+=8

if tracking:
	ids = unpack("i"*lendist,fileContent[counter:counter+lendist*4])
	counter+=lendist*4
	start_index_list = unpack("i"*lendist,fileContent[counter:counter+lendist*4])
	counter+=lendist*4
	stop_index_list = unpack("i"*lendist,fileContent[counter:counter+lendist*4])
	counter+=lendist*4

#The stored bytes in order c++
""" 
	//Open and determine length of integer stream 
	std::ofstream output_bin_file(file_name, std::ios::binary);
	int64_t len = _mat_grid_slice.size();
	output_bin_file.write( (char*)&len, sizeof(len) );

	
	output_bin_file.write( (char*)&_voxel_size, sizeof(real) );
	output_bin_file.write( (char*)&_size_x, sizeof(int32_t) );
	output_bin_file.write( (char*)&_size_y, sizeof(int32_t) );
	output_bin_file.write( (char*)&_save_height, sizeof(int32_t) );

	//Vectors
	output_bin_file.write( (char*)&_mat_grid_slice_int16[0], len * sizeof(int16_t) );
	output_bin_file.write( (char*)&_surface_grid_slice_int16[0], len * sizeof(int16_t) );

	//Tracking
	output_bin_file.write( (char*)&track_sz, sizeof(int32_t) );
	output_bin_file.write( (char*)&track_ids[0], track_sz * sizeof(int) );
	output_bin_file.write( (char*)&start_pos[0], track_sz * sizeof(int) );
	output_bin_file.write( (char*)&end_pos[0], track_sz * sizeof(int) );

    output_bin_file.close();
"""

toc = time.perf_counter()
print(f"Finished Loading in {toc - tic:0.4f} seconds")


#Start Timing Data Manipulation
print("Data Manipulation...")

tic = time.perf_counter()

#Grid Work
mat_grid = np.reshape(mat_grid_in, dim, order='F')
surface_grid = np.reshape(surface_grid_in, dim, order='F')

#Fix Plots --switching x and y axes

mat_grid_t = mat_grid.transpose(1, 0, 2)
surface_grid_t = surface_grid.transpose(1, 0, 2)

# Finding Appropriate Plotting Restrictions-----------------------------------------

#Range Conversions to around 0 (keep pixels same size) (scales)
s_dims = [dim[0]/2,dim[1]/2]
s_dims_nm = [s_dims[0]*voxel_size,s_dims[1]*voxel_size]

#Fixed Width Ranges --------------------------------------------------------

f_x = np.arange(start =-f_range,stop = f_range,step = voxel_size)
f_y = np.arange(start = -f_range,stop =f_range,step = voxel_size)
f_steps = len(f_x)/2
f_x_min = int(s_dims[0]-f_steps+1)
f_x_max = int(s_dims[0]+f_steps)
f_y_min = int(s_dims[1]-f_steps)
f_y_max = int(s_dims[1]+f_steps)
mat_grid_f = mat_grid_t[f_x_min:f_x_max,f_y_min:f_y_max,:]
surface_grid_f = surface_grid_t[f_x_min:f_x_max,f_y_min:f_y_max,:]

f_x = np.arange(start =-f_range,stop = f_range,step = voxel_size)
f_y = np.arange(start = -f_range,stop =f_range,step = voxel_size)

toc = time.perf_counter()
print(f"Finished Data Manipulation in {toc - tic:0.4f} seconds")



#Extra Plots (Disabled for Now)

#Custom Cmap
viridis = cm.get_cmap('viridis', 256)

#Tracking
_size_x = dim[0]
_size_y = dim[1]
_size_z = dim[2]
def index_to_pos(index):
	count = index
	z = math.floor(count/(_size_x*_size_y))
	count -= z*(_size_x*_size_y)
	y = math.floor(count/(_size_y))
	count -= y*(_size_y)
	x = math.floor(count)
	position = np.array([x,y,z])
	return position

def distance_to_source(start_index,stop_index):
	start_pos = index_to_pos(start_index)
	stop_pos = index_to_pos(stop_index)
	diff = stop_pos-start_pos
	dist = math.sqrt(diff[0]**2 + diff[1]**2)
	#dist = math.sqrt(diff[0]**2 + diff[1]**2+diff[2]**2)
	return dist

if extra_plots:
	# Vertical cross section at X = 0 -----------------------------------------

	plt.rcParams['font.size'] = 16
	cs_plot = plt.figure(fignum)
	fignum+=1
	vert_slice_x_f_mat =  np.transpose(mat_grid_f[ : , int((f_x_max-f_x_min)/2), ::1])
	#Process
	vert_slice_x_f_mat[vert_slice_x_f_mat==-123] = -1
	plt.imshow(vert_slice_x_f_mat, extent=(-f_range, f_range, 0, dim[2] * voxel_size))
	plt.xlabel("y (nm)")
	plt.colorbar()
	plt.ylabel("z (nm)")
	#plt.title("Vertical cross section through center at X = 0")
	cs_plot.savefig(cs_mat_vert_x0_str, dpi = 200,  bbox_inches="tight")

	plt.rcParams['font.size'] = 16
	cs_plot = plt.figure(fignum)
	fignum+=1
	vert_slice_x_f_surf =  np.transpose(surface_grid_f[ : , int((f_x_max-f_x_min)/2), ::1])
	plt.imshow(vert_slice_x_f_surf, extent=(-f_range, f_range, 0, dim[2] * voxel_size))
	#plt.colorbar(ticks=[0, 1, 2])
	plt.xlabel("y (nm)")
	plt.ylabel("z (nm)")
	#plt.title("Vertical cross section through center at X = 0")
	cs_plot.savefig(cs_surf_vert_x0_str, dpi = 200,  bbox_inches="tight")

	# Vertical cross section at specific X -----------------------------------------
	nm_displacement = voxel_size

	# px 
	cs_plot = plt.figure(fignum)
	fignum+=1
	plt.rcParams['font.size'] = 16
	vert_slice_x_f =  np.transpose(surface_grid_f[ : , int((f_x_max-f_x_min)/2+nm_displacement/voxel_size),  ::1])
	plt.imshow(vert_slice_x_f, extent=(-f_range, f_range, 0, dim[2] * voxel_size))
	plt.xlabel("y (nm)")
	plt.ylabel("z (nm)")
	#plt.title("Vertical cross section at X = "+str(nm_displacement)+"nm")
	cs_plot.savefig(cs_surf_vert_px_str, dpi = 200,extent=(-f_range, f_range, 0, dim[2] * voxel_size), bbox_inches="tight")

	# mx 
	cs_plot = plt.figure(fignum)
	fignum+=1
	plt.rcParams['font.size'] = 16
	vert_slice_x_f =  np.transpose(surface_grid_f[ : , int((f_x_max-f_x_min)/2-nm_displacement/voxel_size),  ::1])
	plt.imshow(vert_slice_x_f, extent=(-f_range, f_range, 0, dim[2] * voxel_size))
	plt.xlabel("y (nm)")
	plt.ylabel("z (nm)")
	#plt.title("Vertical cross section at X = -"+ str(nm_displacement)+"nm")
	cs_plot.savefig(cs_surf_vert_mx_str, dpi = 200,extent=(-f_range, f_range, 0, dim[2] * voxel_size), bbox_inches="tight")


print("Calculating Distances...")
tic = time.perf_counter()

#Distance Distribution Plot
hom_surface_grid_t = surface_grid_t.copy()
hom_surface_grid_t[hom_surface_grid_t == 5] = 3
adsorbates_index_all = np.where(hom_surface_grid_t == 3)
adsorbates_index_tracked = np.where(surface_grid_t == 5)

if tracking:
	#Implementation to Individual Sources
	min_distance_t_list = []
	for i in range(len(ids)):
		dist_i = distance_to_source(start_index_list[i],stop_index_list[i])
		min_distance_t_list.append(dist_i)
	min_distance_t = np.asarray(min_distance_t_list)

	if source == 'point':
		x_distance = np.abs(s_dims[0] - adsorbates_index_all[0])
		y_distance = np.abs(s_dims[1] - adsorbates_index_all[1])
		min_distance_t = np.sqrt(np.square(x_distance)+np.square(y_distance))

	toc = time.perf_counter()
	print(f"Finished  Calculating Distances in {toc - tic:0.4f} seconds")


	print("Creating Plots...")
	tic = time.perf_counter()

	min_distance_t_nm = min_distance_t*voxel_size
	rms_wall_t_nm = np.sqrt(np.mean(min_distance_t_nm)**2)


	ul_plot = plt.figure(fignum)
	fignum+=1
	bin_num = 20
	plt.rcParams['font.size'] = 16
	if source == 'walls':
		#plt.title(" Minimum Distance of Tracked Adsorbates From Source Distribution", pad = 8)
		plt.hist(min_distance_t_nm, bins = 'auto',edgecolor='black', linewidth=1.2)
		plt.axvline(x=rms_wall_t_nm, color = 'r',label = 'RMS Tracked: '+str(round(rms_wall_t_nm,2))+"nm",linestyle='-.')
		plt.xlim(0,10)

	else:
		#plt.title(" Distance From Source Distribution", pad = 8)
		plt.hist(min_distance_t_nm, bins = 'auto',edgecolor='black', linewidth=1.2)
		plt.axvline(x=rms_wall_t_nm, color = 'r',label = 'RMS Tracked: '+str(round(rms_wall_t_nm,2))+"nm",linestyle='-.')
		plt.xlim(0,20)

	plt.xlabel("Distance from Source (nm)")
	plt.ylabel(" # of Adsorbates")
	plt.legend(loc='upper right')
	ul_plot.savefig(cs_surf_hist, dpi = 200, bbox_inches="tight",extent=(-f_range, f_range, 0, dim[2] * voxel_size))

	#Surface Plot
	ul_plot = plt.figure(fignum)
	ax = ul_plot.gca()
	plt.rcParams['font.size'] = 16
	fignum+=1
	if source =="walls":
		plt.imshow(surface_grid_t[:, :, -1], extent=(-s_dims_nm[0], s_dims_nm[0], -s_dims_nm[1], s_dims_nm[1]))
	else:
		plt.imshow(hom_surface_grid_t[:, :, -1], extent=(-s_dims_nm[0], s_dims_nm[0], -s_dims_nm[1], s_dims_nm[1]))
	plt.xlabel("x (nm)")
	plt.ylabel("y (nm)")
	if source == 'walls':
		plt.axvline(x=-s_dims_nm[0]+rms_wall_t_nm, color = 'r',linestyle='-.')
		plt.axvline(x=s_dims_nm[0]-rms_wall_t_nm, color = 'r',linestyle='-.')
		plt.axhline(y=-s_dims_nm[1]+rms_wall_t_nm, color='r',linestyle='-.')
		plt.axhline(y=s_dims_nm[1]-rms_wall_t_nm, color='r',label = 'RMS Tracked: '+str(round(rms_wall_t_nm,2))+"nm",linestyle='-.')
	elif source == 'point':
		circ = plt.Circle((0,0),radius=rms_wall_t_nm,color ='r',fill=False,label = 'RMS: '+str(round(rms_wall_t_nm,2))+"nm",linestyle = '--')
		ax.add_patch(circ)
	plt.legend(loc='upper right')
	#plt.title(" Non Deposition 2D Surface Diffusion", pad = 8)

	ul_plot.savefig(cs_surf_hor_rms, dpi = 200, bbox_inches="tight",extent=(-f_range, f_range, 0, dim[2] * voxel_size))

	#Timing and Plotting

toc = time.perf_counter()
print(f"Finished  Plots in {toc - tic:0.4f} seconds")







