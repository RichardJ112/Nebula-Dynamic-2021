import numpy as np
import os
import pickle
import sys
from statistics import median 
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import rc
from matplotlib import cm
from mpl_toolkits import mplot3d
import matplotlib.patches as mpatches
import matplotlib
import time
import matplotlib.ticker as ticker
from struct import *

# This file is for reading line,pillar and wall output binary deposition geometries - cones handled seperately
# This version removes any use of the old classifications / streamlines some aspects

# IO parameters
date = "/August/30_8_2021/"
#date = "/July/10_7_2021/"
#deposit_type = 'point'#types supported are line,diagonal,point,wall
deposit_type = 'point'
material = 'tungsten'
#material = 'graphite-phonon'
#cross_section_type = '_AWM_'
#cross_section_type = '_SMC_'
cross_section_type = '_SMC_'
HPC_toggle = False #Set to true when running on HPC
direct_pathing = False #allow for direct input of file through command line
desktop_toggle = True
enable_e_grid = False
gas_handling = True
nm_check = 2
f_range = 20#fixed range in nm for plotting section of full domain
f_width_auto = 0
max_toggle = True

#Plotting Toggles
surface_map_3d_toggle = True
height_map_3d_toggle = False
class_legend_toggle = False

#gas_handling_str = "_e20_em9"
#gas_handling_str = "_e18_em8_b"
#gas_handling_str = ""
gas_handling_str = "_e17_em9_mb"

#output_extra_str = ""
output_extra_str = ""

#parameter_summary ="1keV_1_10kpp_pitch_0_401_401_1001_sb_3000_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_" #point small pillar
parameter_summary = "1keV_1_500kpp_pitch_0_161_161_1001_sb_1273_vs_250_sd_34000_sh_996_detect_dome_mirror_vs_250_" #point
#parameter_summary = "1keV_1_400kpp_pitch_0_161_161_1001_sb_1000_vs_357_sd_34000_sh_996_detect_dome_mirror_vs_357_" #point
#parameter_summary ="1keV_1_100kpp_pitch_0_401_401_1001_sb_3000_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_" #point
#parameter_summary = "10keV_1_100kpp_pitch_0_401_401_1001_sb_1000_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_" #point
#parameter_summary = "1keV_1_10kpp_pitch_0_401_401_1001_sb_1000_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_" #point
#parameter_summary = "5keV_21_5kpp_pitch_3100_seq_p_401_401_1001_sb_1000_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_" #line
#parameter_summary = "1keV_21_1kpp_pitch_3000_seq_l_401_401_1001_sb_1000_sd_34000_sh_996_detect_dome_mirror_vs_300_" #line
#parameter_summary = "1keV_21_1kpp_pitch_2400_seq_p_401_401_1001_sb_1000_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_" #line
#parameter_summary = "1keV_21_1kpp_pitch_800_seq_p_401_401_1001_sb_1000_sd_34000_sh_996_detect_dome_mirror_vs_300_" #line
#parameter_summary = "1keV_11_21_1kpp_pitchx_6000_pitchy_1700_seq_d_401_401_1001_sb_1000_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_" #wall
#parameter_summary = "1keV_11_21_1kpp_pitchx_4000_pitchy_1700_seq_ll_401_401_1001_sb_1000_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_" #wall


#Input Paths -- User PSA: Change to ones relevant for your environments
if HPC_toggle: 
    # HPC Path
    path_to_test_files = "/home/richarddejong/nebula_test_files/"
elif desktop_toggle:
    # Desktop Path
    path_to_test_files = "C:/Users/Richard/source/repos/Nebula/nebula_test_files/"
else: 
    # Laptop Path
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
energy_str = title_list[0]

# Shape Specific Interpreter
if deposit_type == 'point':
    electrons_per_pillar = int(title_list[2].replace('kpp',''))*1000
    #electrons_per_pillar = 3500 #override
    pillar_num = int(title_list[1])
    beam_std_nm = int(title_list[9])/1000
    pitch = float(int(title_list[4])/1000)
    pitch_nm = float(int(title_list[4])/1000)

    electron_num = pillar_num * electrons_per_pillar
    f_width_auto = np.around(beam_std_nm*2/5, decimals=0)*5
    f_range = 20 #override
    
if deposit_type == 'line':

    electrons_per_pillar = int(title_list[2].replace('kpp',''))*1000
    pillar_num = int(title_list[1])
    pitch = float(int(title_list[4])/1000)
    pitch_nm = float(int(title_list[4])/1000)
    
    electron_num = pillar_num * electrons_per_pillar
    dep_length = (pitch_nm*pillar_num/2)+5
    f_width_auto = np.around(dep_length/5, decimals=0)*5
    f_range = f_width_auto

if deposit_type == 'wall':
    pitch_x_nm = float(int(title_list[5])/1000)
    pitch_y_nm = float(int(title_list[7])/1000)

    pillar_num_x = int(title_list[1])
    pillar_num_y = int(title_list[2])
    pillar_num = pillar_num_x*pillar_num_y

    electrons_per_pillar = int(title_list[3].replace('kpp',''))*1000

    electron_num = pillar_num * electrons_per_pillar 
    dep_length = (pitch_x_nm*pillar_num_x/2)+15
    f_width_auto = np.around(np.ceil(dep_length/5), decimals=0)*5
    f_range = f_width_auto

    #f_range = 45 #override

output_path_folder = path_to_test_files+"figures"+date+title
#output_path_folder = "/mnt/c/Users/richa/Documents/repos/nebula_test_files/figures"+date+title #WSL
  
general_path = output_path_folder+"/"+title
fignum = 1

if not os.path.exists(output_path_folder):
    os.makedirs(output_path_folder)

# General path strings
td_profile_str = general_path+"_TD_profile" #top down profile of growth
t_evo_x_section_countour_str =  general_path+"_time_evolution_X_contours" #countour string of time evolution in vertical cross section(contours)
t_evo_y_section_countour_str =  general_path+"_time_evolution_Y_contours" #countour string of time evolution in vertical cross section(contours)
t_evo_x_str = general_path+"_time_evolution_X" #time evolution in vertical cross section
t_evo_y_str = general_path+"_time_evolution_Y" #time evolution in vertical cross section
cs_hor_surface_str_new = general_path +"_cs_hor_surface_new"# horizontal cross section near surface
cs_hor_surface0_str_new = general_path +"_cs_hor_halo_new"# horizontal cross section near surface
t_evo_surface = general_path +"_surface_time_evolution"# horizontal cross section near surface
t_evo_halo = general_path +"_surface_halo_time_evolution"# horizontal cross section at z = 0nm
cs_vert_x0_str_new = general_path +"_cs_vert_X0_new" # vertical cross section at X=0(New Classification)
cs_vert_y0_str_new = general_path +"_cs_vert_Y0_new" # vertical cross section at X=0(New Classification)
cs_vert_x_str_new = general_path +"_cs_vert_X_new" # vertical cross section at specific X (New Classification)
rad_dist_str = general_path + "_rad_dist" # radial distribution for pillar
lat_dist_str = general_path + "_lat_dist" # lateral distribution for lines

#Experimental path strings
height_map_3d_str = general_path + "_3d_height_map"
surf_map_3d_str = general_path + "_3d_surf_map"


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
mat_grid = np.reshape(mat_grid_in, dim, order='F')/123 + 1
tag_grid = np.reshape(tag_grid_in, dim, order='F')   
new_species_grid = np.reshape(new_species_in, dim, order='F')

# Debug new species
max_new_species = np.max(new_species_grid)
min_new_species = np.min(new_species_grid)

#Fix Plots --switching x and y axes

mat_grid_t = mat_grid.transpose(1, 0, 2)
tag_grid_t = tag_grid.transpose(1, 0, 2)
new_species_grid_t = new_species_grid.transpose(1,0,2)
profile = np.sum(mat_grid_t, axis=2)

#Fixed Width Ranges --------------------------------------------------------
if f_width_auto>(dim[1]*voxel_size/2): # small correction in case of overreach
    f_range -= 5
f_x = np.arange(start =-f_range,stop = f_range,step = voxel_size)
f_y = np.arange(start = -f_range,stop =f_range,step = voxel_size)
f_steps = len(f_x)/2
s_dims = [dim[0]/2,dim[1]/2]
f_x_min = int(s_dims[0]-f_steps)
f_x_max = int(s_dims[0]+f_steps)
f_y_min = int(s_dims[1]-f_steps)
f_y_max = int(s_dims[1]+f_steps)

mat_grid_f = mat_grid_t[f_x_min:f_x_max,f_y_min:f_y_max,:]
tag_grid_f =  tag_grid_t[f_x_min:f_x_max,f_y_min:f_y_max,:]
new_species_grid_f = new_species_grid_t[f_x_min:f_x_max,f_y_min:f_y_max,:]

f_x = np.arange(start =-f_range,stop = f_range,step = voxel_size)
f_y = np.arange(start = -f_range,stop =f_range,step = voxel_size)

if enable_e_grid:
    e_grid = np.reshape(e_grid_in, dim, order='F')  
    e_grid_t = e_grid.transpose(1, 0, 2)
    e_grid_f =  e_grid_t[f_x_min:f_x_max,f_y_min:f_y_max,:]

toc = time.perf_counter()
print(f"Finished Data Manipulation in {toc - tic:0.4f} seconds")


#Functions

# Height Map Generation
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

def surface_index(tag_grid,dim,remove_surface = 0.4,remove_surface_toggle=True,range_limit_toggle = True, remove_range = f_range):
    """ 
    from tag_grid create 3d numpy map with surface indices of each column

    Parameters
    ----------
    tag_grid_f : primary tag for each deposited electron

    Returns
    -------
    plottable indices : 3d numpy array of all surfaces

    """

    mask = tag_grid!=0

    mask_m = np.rot90(np.rot90(mask))
    # Positive Direction
    x_indices = dim[0] - np.where(mask.any(axis=0), mask.argmax(axis=0), dim[0])
    y_indices = dim[1] - np.where(mask.any(axis=1), mask.argmax(axis=1), dim[1])
    z_indices = dim[2] - np.where(mask.any(axis=2), mask.argmax(axis=2), dim[2])

    xm_indices = dim[0] - np.where(mask_m.any(axis=0), mask_m.argmax(axis=0), dim[0])
    ym_indices = dim[1] - np.where(mask_m.any(axis=1), mask_m.argmax(axis=1), dim[1])

    #Zeros
    x_p_z = np.nonzero(x_indices)
    y_p_z = np.nonzero(y_indices)
    z_p_z = np.nonzero(z_indices)

    x_m_z = np.nonzero(xm_indices)
    y_m_z = np.nonzero(ym_indices)

    #Surface Indices (Dont forget that x and y are inverted due to z definition)
    x_p_i = np.vstack((x_p_z[0]-dim[1]/2,dim[0]/2-x_indices[x_p_z],dim[2]-x_p_z[1]))
    y_p_i = np.vstack((y_indices[y_p_z]-dim[1]/2,y_p_z[0]-dim[0]/2,dim[2]-y_p_z[1],))
    z_p_i = np.vstack((z_p_z[1]-dim[1]/2,z_p_z[0]-dim[0]/2,z_indices[z_p_z]))

    x_m_i = np.vstack((x_m_z[0]-dim[1]/2,xm_indices[x_m_z]-dim[1]/2,dim[2]-x_m_z[1]))
    y_m_i = np.vstack((ym_indices[y_m_z]-dim[0]/2,dim[1]/2-y_p_z[0],dim[2]-y_m_z[1]))

    #Large Stack
    #surface_i = np.hstack((x_p_i,y_p_i,z_p_i,x_m_i,y_m_i))*voxel_size # all scans
    surface_i = np.hstack((x_p_i,y_p_i,z_p_i,x_m_i))*voxel_size #no y_m
    #surface_i = np.hstack((x_p_i,y_p_i,z_p_i))*voxel_size #only positive scans
    #surface_i = z_p_i*voxel_size
    #surface_i = np.hstack((y_p_i,y_m_i))*voxel_size #only positive scans
    #surface_i = y_m_i*voxel_size

    #Removes excess points lower than remove surface
    if remove_surface_toggle:
        removal_i = np.where(surface_i[2]<remove_surface)
        surface_i = np.delete(surface_i,removal_i,axis=1)
    
    #Removes excess points past axis limits in xy
    if range_limit_toggle:
        removal_x = np.where(np.abs(surface_i[0])>remove_range)
        surface_i = np.delete(surface_i,removal_x,axis=1)
        removal_y = np.where(np.abs(surface_i[1])>remove_range)
        surface_i = np.delete(surface_i,removal_y,axis=1)
    return surface_i

print("Creating Plots...")
tic = time.perf_counter()
#ColorMap comparision

"""
    Old Dingemanse:(colours for cmap_comp)
    1 == empty white
    dist_PE = R[electron_type == 2] red
    dist_FSE = R[electron_type == 3] orange
    dist_BSE = R[electron_type == 4] blue
    dist_VE = R[electron_type == 5] purple
    dist_SE_PE = R[electron_type == 6] #ffa020
    dist_SE_FSE = R[electron_type == 7]  #80ff80
    dist_SE_BSE = R[electron_type == 8]  cyan
    dist_SE_VE = R[electron_type == 9] magenta

    New:(colours for cmap_comp)
    species_labels = 
    ['0 Empty', white
    ' 1 PE', red
    '2 PE_FSE_D', green
    '3  PE_BSE_D', blue
    '4  PE_BSE_S', blue
    '5  PE_FSE_DVD', #dc5858
    '6  PE_BSE_DVD', blue
    '7  PE_BSE_SVD',blue
    '8  SE_PE', #ffa020
    '9  SE_FSE', '#80ff80'
    '10  SE_BSE_D', magenta
    '11  SE_BSE_S', purple
    '12  Uncategorized'] gray
"""
cmap_comp = LinearSegmentedColormap.from_list('mycmap', ['white', 'red', 'green','blue','blue','#dc5858','blue','blue','#ffa020','#80ff80','magenta','purple','gray']) #13 colors for Smith
cmap = LinearSegmentedColormap.from_list('mycmap', ['white', 'red', 'blue', 'cyan', 'green', 'xkcd:olive brown'])
cmap8 = LinearSegmentedColormap.from_list('mycmap', ['white', 'gray', 'red', 'green', 'blue', 'purple', '#FFA020', '#80FF80', 'cyan', 'magenta']) #10 colours for 8 types
cmap12 = LinearSegmentedColormap.from_list('mycmap', ['white', 'gray', 'red', 'green', 'blue', 'purple', '#FFA020', '#80FF80', 'cyan', 'magenta','orange','black','#301713']) #13 colors for 13 types

# Height Map
plt.rcParams['font.size'] = 16
height_map_str = general_path + "_height_map"
h_map_t = height_map(tag_grid_t,axis=2)
h_map_f = height_map(tag_grid_f,axis=2)
plt.figure(fignum)
fignum += 1
max_h = np.max(h_map_t)
plt.imshow(h_map_f, extent=(-f_range, f_range, -f_range, f_range), origin={'upper', 'upper'})
plt.colorbar(shrink = 1,label = "Height(nm)")
plt.xlabel("x (nm)")
plt.ylabel("y (nm)")
#plt.title(r"Height Map", pad= 15)
plt.savefig(height_map_str, aspect = 'auto', bbox_inches="tight")

#Material profile -----------------------------------------------------------------
plt.figure(fignum)
fignum += 1
plt.rcParams['font.size'] = 16
profile = np.sum(mat_grid_t, axis=2)
profile_f = np.sum(mat_grid_f, axis =2)

plt.imshow(profile_f*voxel_size, extent=(-f_range, f_range, -f_range, f_range), origin={'upper', 'upper'})
plt.xlabel("x (nm)")
plt.ylabel("y (nm)")
#plt.title(r"Material Sum Map", pad= 15)
plt.colorbar(shrink = 1,label = "Pixel Material Sum(nm)")
plt.savefig(td_profile_str, aspect = 'auto', bbox_inches="tight")

#3D Height Map --------------------------------------------------------------------------------

if height_map_3d_toggle:
    plt.figure(fignum)
    fignum += 1
    x = np.linspace(-voxel_size*dim[0]/2,voxel_size*dim[0]/2, dim[0])
    y = np.linspace(-voxel_size*dim[1]/2,voxel_size*dim[1]/2, dim[1])
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='ij')
    xv = np.ravel(xv)
    yv = np.ravel(yv)
    z = np.ravel(h_map_t)
    ax = plt.axes(projection='3d')
    ax.scatter(xv, yv, z,c=z, cmap='viridis', linewidth=0.5)
    plt.savefig(height_map_3d_str, aspect = 'auto', bbox_inches="tight")

#3D Surface Map ---------------------------------------------------------------------------------

if surface_map_3d_toggle:
    plt.figure(fignum)
    fignum += 1
    surface_list = surface_index(tag_grid_t,dim,remove_surface=1,remove_range=f_range)
    ax = plt.axes(projection='3d')
    #plt.title("Deposit 3D Surface Scan", pad= 15)
    ax.view_init(elev=30, azim=60)
    ax.scatter(surface_list[0], surface_list[1], surface_list[2],c=surface_list[2], cmap='viridis', linewidth=0.5,alpha=0.1)
    #ax.plot_trisurf(surface_list[0], surface_list[1], surface_list[2], cmap='viridis', linewidth=0.3)
    #plt.show()
    ax.set_xlim(-f_range,f_range)
    if deposit_type == 'line':
        #ax.set_xlim(-f_range/3,f_range/3)
        ax.set_xlim(-f_range,f_range)
    ax.set_ylim(-f_range,f_range)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_zlabel("z (nm)")
    plt.xticks(np.arange(-f_range, f_range+1, f_range/2))
    plt.yticks(np.arange(-f_range, f_range+1, f_range/2))
    plt.tight_layout()
    plt.savefig(surf_map_3d_str, aspect = 'auto',bbox_inches="tight")

# Horizontal cross section at surface (new classification)-----------------------------------------------

ul_plot = plt.figure(fignum)
fignum+=1
plt.rcParams['font.size'] = 16
plt.imshow(np.flipud(new_species_grid_f[:, :, -1]), cmap=cmap_comp, extent=(-f_range, f_range, -f_range, f_range),interpolation = "nearest",vmin= 0,vmax = 12)
plt.xlabel("x (nm)")
plt.ylabel("y (nm)")
#plt.title(" Horizontal cross section at surface", pad = 8)
ul_plot.savefig(cs_hor_surface0_str_new, dpi = 200, bbox_inches="tight")

# Time evolution graph at a specific X Contours*-------------------------------------------------------
tag_cs_vert_x = np.transpose(tag_grid_f[ : , int((f_x_max-f_x_min)/2), :])
time_evo_plot_X_C = plt.figure(fignum)
fignum+=1
plt.rcParams['font.size'] = 16

if pillar_num == 1:
    electron_increment = int(electron_num/10) 
    contour_levels = np.arange(start=0, stop=electron_num+1, step=electron_increment)
if deposit_type == 'line':
    contour_levels  = np.arange(start=0, stop=electron_num+1, step=electrons_per_pillar)
if deposit_type == 'wall':
    electron_increment = int(electron_num/pillar_num_y) 
    contour_levels = np.arange(start=0, stop=electron_num+1, step=electron_increment)

f_z = np.arange(start =0,stop = (dim[2])*voxel_size,step =voxel_size)
contour_Z,contour_Y = np.meshgrid(f_y,f_z)

#plt.title("Time evolution at X = 0 Section Countour Plot")
plt.contour(contour_Z,contour_Y,np.flipud(tag_cs_vert_x),levels = contour_levels,colors = 'black',linewidths = 0.3)
plt.contourf(contour_Z,contour_Y,np.flipud(tag_cs_vert_x),levels = contour_levels)
plt.colorbar(shrink = 1,label = "Primary Electron ID")
plt.xlabel("y (nm)")
plt.ylabel("z (nm)")
time_evo_plot_X_C.savefig(t_evo_x_section_countour_str, dpi = 200,bbox_inches="tight")

# Time evolution graph at a specific Y Contours*-------------------------------------------------------
tag_cs_vert_y = np.transpose(tag_grid_f[ int((f_x_max-f_x_min)/2) , :, :])
time_evo_plot_Y_C = plt.figure(fignum)
fignum+=1 
plt.rcParams['font.size'] = 16

plt.contour(contour_Z,contour_Y,np.flipud(tag_cs_vert_y),levels = contour_levels,colors = 'black',linewidths = 0.3)
plt.contourf(contour_Z,contour_Y,np.flipud(tag_cs_vert_y),levels = contour_levels)
plt.colorbar(shrink = 1,label = "Primary Electron ID")
plt.xlabel("x (nm)")
plt.ylabel("z (nm)")
#plt.title("Time evolution at Y = 0 Section Countour Plot")
time_evo_plot_Y_C.savefig(t_evo_y_section_countour_str, dpi = 200,bbox_inches="tight")


# Time evolution graph at a specific Z -------------------------------------------------------

ul_plot = plt.figure(fignum)
fignum+=1
plt.rcParams['font.size'] = 16
plt.imshow(np.flipud(tag_grid_f[:, :, -1]), extent=(-f_range, f_range, -f_range, f_range))
plt.xlabel("x (nm)")
plt.ylabel("y (nm)")
plt.colorbar(shrink = 1,label = " Primary Electron ID")
#plt.title(" Time evolution at Z = 1nm above surface", pad = 8)
ul_plot.savefig(t_evo_surface, dpi = 200,bbox_inches="tight")

# Cross Sections Electron Classification -------------------------------------------------------------------------------

# Vertical cross section at X=0 (New Classification)-----------------------------------------
cs_plot = plt.figure(fignum)
fignum+=1
plt.rcParams['font.size'] = 16
vert_slice_x_f_new =  np.transpose(new_species_grid_f[ : , int((f_x_max-f_x_min)/2), ::-1])
im = plt.imshow(vert_slice_x_f_new,cmap = cmap_comp,extent=(-f_range, f_range, 0, dim[2] * voxel_size),origin='lower',interpolation = "nearest",vmin= 0,vmax = 12)


if class_legend_toggle:
    # Legend
    # get the colors of the values, according to the 
    # colormap used by imshow
    values = np.linspace(0,12,num = 13)
    colors = [ im.cmap(im.norm(value)) for value in values]
    # create a patch (proxy artist) for every color 
    patches = [ mpatches.Patch(color=colors[i], label="Level {l}".format(l=values[i]) ) for i in range(len(values)) ]
    # put those patched as legend-handles into the legend
    L = plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0,prop={'size': 4} )
    #L = plt.legend(handles=patches)
 
    species_labels = ['0 Empty',
    ' 1 PE',
    '2 PE_FSE_D',
    '3  PE_BSE_D',
    '4  PE_BSE_S',
    '5  PE_FSE_DVD',
    '6  PE_BSE_DVD',
    '7  PE_BSE_SVD',
    '8  SE_PE',
    '9  SE_FSE',
    '10  SE_BSE_D',
    '11  SE_BSE_S', 
    '12  Uncategorized']


    for i in range(len(species_labels)):
        L.get_texts()[i].set_text(species_labels[i])

plt.xlabel("y (nm)")
plt.ylabel("z (nm)")
#plt.title("Vertical cross section through center at X = 0")
cs_plot.savefig(cs_vert_x0_str_new, dpi = 200, bbox_inches="tight")

# Vertical cross section at Y=0 (New Classification)-----------------------------------------
cs_plot = plt.figure(fignum)
fignum+=1
plt.rcParams['font.size'] = 16
vert_slice_y_f =  np.transpose(new_species_grid_f[ int((f_x_max-f_x_min)/2) , : , ::-1])
plt.imshow(vert_slice_y_f, cmap=cmap_comp, extent=(-f_range, f_range, 0, dim[2] * voxel_size), origin='lower',vmin=0, vmax = 12)
#plt.ylim(50)
plt.xlabel("x (nm)")
plt.ylabel("z (nm)")
#plt.title("Vertical cross section through center at Y = 0")
cs_plot.savefig(cs_vert_y0_str_new, dpi = 200, bbox_inches="tight")

# Vertical cross section at specific X (New Classification)-----------------------------------------
"""
nm_displacement = 2
plt.rcParams.update({'font.size': 12})
cs_plot = plt.figure(fignum)
fignum+=1
vert_slice_x_f_new =  np.transpose(new_species_grid_f[ : , int((f_x_max-f_x_min)/2+nm_displacement/voxel_size), ::-1])
im = plt.imshow(vert_slice_x_f_new, cmap = cmap_comp,extent=(-f_range, f_range, 0, dim[2] * voxel_size), origin='lower',interpolation = "nearest",vmin= 0,vmax = 12)

if class_legend_toggle:
    # Legend
    # get the colors of the values, according to the 
    # colormap used by imshow
    values = np.linspace(0,12,num = 13)
    colors = [ im.cmap(im.norm(value)) for value in values]
    # create a patch (proxy artist) for every color 
    patches = [ mpatches.Patch(color=colors[i], label="Level {l}".format(l=values[i]) ) for i in range(len(values)) ]
    # put those patched as legend-handles into the legend
    L = plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0,prop={'size': 10} )


    for i in range(len(species_labels)):
        L.get_texts()[i].set_text(species_labels[i])

plt.xlabel("y (nm)")
plt.ylabel("z (nm)")
plt.title("Vertical cross section at X = 2nm")
cs_plot.savefig(cs_vert_x_str_new, dpi = 200, bbox_inches="tight")
"""
# Type Counts -----------------------------------------------------------------------------------------------

unique, counts = np.unique(new_species_grid_t, return_counts=True)
species_dict = dict(zip(unique, counts))
total_deps = np.sum(counts) -counts[0]

# Smith et Al Amount Comparisons (Using Smith definitions of PE,FSE,BSE,SE1,SE2)
S_PE = np.sum(new_species_grid_t == 1)
S_FSE = np.sum(new_species_grid_t == 2) + np.sum(new_species_grid_t == 5)
S_BSE = np.sum(new_species_grid_t == 3) + np.sum(new_species_grid_t == 4) + np.sum(new_species_grid_t == 6) + np.sum(new_species_grid_t == 7)
S_SE_1 = np.sum(new_species_grid_t == 8) + np.sum(new_species_grid_t == 9) #secondary electrons from intial traversal(PE,FSE)
S_SE_2 = np.sum(new_species_grid_t == 10) +np.sum(new_species_grid_t == 11) # secondary electrons from alternate traversal (BSE)
all_deps = np.sum(new_species_grid_t>0)

depo_eff = (total_deps/electron_num)*100
pillar_heights = []
pillar_widths = []
pillar_widths_h1 =[]
pillar_HDW = []
pillar_FWHM = []
pillar_depths = []
pillar_widths_nm_check =[]
incline = str(0)
midpoint = int(dim[1]/2)

#Radial Dist for Point Depositions for New Classification ---------------------------------------------------
spatial_dist = True
if (deposit_type == 'point') and (spatial_dist):
    x = np.linspace(-voxel_size*dim[0]/2,voxel_size*dim[0]/2, dim[0]) 
    y = np.linspace(-voxel_size*dim[1]/2,voxel_size*dim[1]/2, dim[1])

    bin_number = 200

    X, Y, Z = np.meshgrid(x, y, np.ones(dim[2]))

    R = np.sqrt(X**2 + Y**2)

    dist_all = R[new_species_grid != 0]
    dist_PE = R[new_species_grid == 1]
    dist_FSE = R[(new_species_grid == 2) | (new_species_grid == 5)]
    dist_BSE = R[new_species_grid == 3]
    dist_SE_1 = R[(new_species_grid == 8) | (new_species_grid == 9)]
    dist_SE_2 = R[(new_species_grid == 10) | (new_species_grid == 11)]

    plt.rcParams.update({'font.size': 16})


    hist = np.zeros((6, bin_number))

    hist_range = voxel_size*(int)(dim[0]/2)
    #hist_range = f_range
    hist[0, :], bins = np.histogram(dist_all, bins= bin_number, range=(0, hist_range))
    hist[1, :], bins = np.histogram(dist_PE, bins= bin_number, range=(0, hist_range))
    hist[2, :], bins = np.histogram(dist_FSE, bins= bin_number, range=(0, hist_range))
    hist[3, :], bins = np.histogram(dist_BSE, bins= bin_number, range=(0, hist_range))
    hist[4, :], bins = np.histogram(dist_SE_1, bins= bin_number, range=(0, hist_range))
    hist[5, :], bins = np.histogram(dist_SE_2, bins= bin_number, range=(0, hist_range))

    areas = np.zeros(bin_number)

    for i in range(bin_number):
        areas[i] = np.pi*(bins[i+1]**2 - bins[i]**2)

    max_r = voxel_size*(int)(dim[0]/4)

    plt.figure(fignum)
    fignum+=1

    voxel_num_hist, bins = np.histogram(R[:, :, 0].flatten(), bins= bin_number, range=(0, voxel_size*(int)(dim[0]/2)))

    radial_zone = np.linspace(max_r*0.01, max_r*0.99, bin_number)
    plt.plot(radial_zone, hist[0, :]/voxel_num_hist, color = 'black' ) #all depsoitions
    plt.plot(radial_zone, hist[1, :]/voxel_num_hist, color = 'red' )
    plt.plot(radial_zone, hist[2, :]/voxel_num_hist, color = 'green' )
    plt.plot(radial_zone, hist[3, :]/voxel_num_hist, color = 'blue' )
    plt.plot(radial_zone, hist[4, :]/voxel_num_hist, color = '#ffa020' )
    plt.plot(radial_zone, hist[5, :]/voxel_num_hist, color = 'magenta' )

    plt.xlim([0,f_range/2])
    plt.xlabel("Pillar Radius (nm)")
    plt.rcParams.update({'font.size': 16})
    plt.ylabel(r"Deposit Intensity (nm$^{-2}$)")
    labels_distro = ["All","PE","FSE","BSE","SE I","SE II"]
    plt.legend(labels_distro)
    #plt.title("Radial Distribution of Deposition Electrons\n"+energy_str+" Beam")

    plt.savefig(rad_dist_str, bbox_inches="tight")

#Lateral Dist for Point Depositions for New Classification
if (deposit_type == 'line') and (spatial_dist):
    x = np.linspace(-voxel_size*dim[0]/2,voxel_size*dim[0]/2, dim[0]) 
    y = np.linspace(-voxel_size*dim[1]/2,voxel_size*dim[1]/2, dim[1])

    bin_number = 200

    X, Y, Z = np.meshgrid(x, y, np.ones(dim[2]))

    L = Y

    dist_all = L[new_species_grid != 0]
    dist_PE = L[new_species_grid == 1]
    dist_FSE = L[(new_species_grid == 2) | (new_species_grid == 5)]
    dist_BSE = L[new_species_grid == 3]
    dist_SE_1 = L[(new_species_grid == 8) | (new_species_grid == 9)]
    dist_SE_2 = L[(new_species_grid == 10) | (new_species_grid == 11)]

    plt.rcParams.update({'font.size': 16})


    hist = np.zeros((6, bin_number))

    hist_range = voxel_size*(int)(dim[0]/2)
    hist[0, :], bins = np.histogram(dist_all, bins= bin_number, range=(-hist_range, hist_range))
    hist[1, :], bins = np.histogram(dist_PE, bins= bin_number, range=(-hist_range, hist_range))
    hist[2, :], bins = np.histogram(dist_FSE, bins= bin_number, range=(-hist_range, hist_range))
    hist[3, :], bins = np.histogram(dist_BSE, bins= bin_number, range=(-hist_range, hist_range))
    hist[4, :], bins = np.histogram(dist_SE_1, bins= bin_number, range=(-hist_range, hist_range))
    hist[5, :], bins = np.histogram(dist_SE_2, bins= bin_number, range=(-hist_range, hist_range))

    areas = np.zeros(bin_number)

    for i in range(bin_number):
        areas[i] = np.pi*(bins[i+1]**2 - bins[i]**2)

    max_r = voxel_size*(int)(dim[0]/4)

    plt.figure(fignum)
    fignum+=1

    voxel_num_hist, bins = np.histogram(L[:, :, 0].flatten(), bins= bin_number, range=(-hist_range, hist_range))

    lateral_zone = np.linspace(-hist_range, hist_range, bin_number)
    plt.plot(lateral_zone, hist[0, :]/voxel_num_hist, color = 'black' ) #all depsoitions
    plt.plot(lateral_zone, hist[1, :]/voxel_num_hist, color = 'red' )
    plt.plot(lateral_zone, hist[2, :]/voxel_num_hist, color = 'green' )
    plt.plot(lateral_zone, hist[3, :]/voxel_num_hist, color = 'blue' )
    plt.plot(lateral_zone, hist[4, :]/voxel_num_hist, color = '#ffa020' )
    plt.plot(lateral_zone, hist[5, :]/voxel_num_hist, color = 'magenta' )

    x_range = 10
    plt.xlim([-x_range,x_range])
    plt.xlabel("Lateral Shift (nm)")
    plt.rcParams.update({'font.size': 16})
    plt.ylabel(r"Deposit Intensity (nm$^{-2}$)")
    labels_distro = ["All","PE","FSE","BSE","SE I","SE II"]
    plt.legend(labels_distro)
    plt.title("Lateral Distribution of Deposition Electrons\n"+energy_str+" Beam")

    plt.savefig(lat_dist_str, bbox_inches="tight")

# Parameter Studies Universal-------------------------------------------------------------------------------

#Function that finds first non-sequential index in list of indices and pillar values
def find_lower_bound_index(pillar_vals,zero_index):
    zero_indices = zero_index[0]
    for i in range(len(zero_indices)-1):
        if zero_indices[i+1] != zero_indices[i]+1:
            if zero_indices[i+1]+1 == len(pillar_vals):
                return 0
            if pillar_vals[zero_indices[i+1]+1] != 0:
                lsi = int(zero_indices[i+1])
                return lsi
    return 0


# Define HDW - Half Diameter Width  -- 
vert_slice_x =  np.transpose(tag_grid_t[ : , midpoint, ::-1]) #takes vertical_slice at pillar position along x = 0
if deposit_type == 'point':
    sep = pitch/voxel_size
    i = 1
    j = 0
    pillar_vals = vert_slice_x[:,midpoint]# takes single column at x = pillar_pos , y = 0 
    zero_index = np.nonzero(pillar_vals)
    if zero_index[0].size <= 3:
        pillar_heights.append(0)
        pillar_FWHM.append(0)
        pillar_widths_h1.append(0)
    height_index = np.max(zero_index) 
    pillar_height_i = height_index*voxel_size
    if max_toggle:
        pillar_heights.append(np.max(h_map_t))
    else:
        pillar_heights.append(pillar_height_i)
    lower_bound_index = 0
    half_diameter_index = int((height_index-lower_bound_index)/2 + lower_bound_index) #index of half diameter
    row_at_1nm = vert_slice_x[int(1/voxel_size)]
    row_at_nm_check = vert_slice_x[int(nm_check/voxel_size)]
    row_half = vert_slice_x[half_diameter_index]
    zero_index_row_nm_check = np.nonzero(row_at_nm_check)
    zero_index_row_1nm = np.nonzero(row_at_1nm)
    zero_index_row_half =  np.nonzero(row_half)
    if pillar_height_i/2 > 1:
        pillar_FWHM.append((np.max(zero_index_row_half)-np.min(zero_index_row_half))*voxel_size)
    else:
        pillar_FWHM.append(0)
    if zero_index_row_1nm[0].size == 0:
        pillar_widths_h1.append(0)
    else:
        pillar_widths_h1.append((np.max(zero_index_row_1nm)-np.min(zero_index_row_1nm))*voxel_size)
    if zero_index_row_nm_check[0].size == 0:
        pillar_widths_nm_check.append(0)
    else:
        pillar_widths_nm_check.append((np.max(zero_index_row_nm_check)-np.min(zero_index_row_nm_check))*voxel_size)
        
if deposit_type == 'line':
    sep = pitch/voxel_size
    for i in range(pillar_num):
        j = (-pillar_num+1)/2+i
        new_pos = int(midpoint +j*sep)
        pillar_vals = vert_slice_x[:,new_pos]# takes single column at x = pillar_pos , y = 0 
        zero_index = np.nonzero(pillar_vals)
        pillar_height_i = h_map_t[new_pos,midpoint]
        if zero_index[0].size <= 3:
            pillar_HDW.append(0)
            continue
        pillar_heights.append(pillar_height_i)
        height_index = int(pillar_height_i/voxel_size)
        lower_bound_index = find_lower_bound_index(pillar_vals,zero_index)
        depth_nm = pillar_height_i-lower_bound_index*voxel_size
        pillar_depths.append(depth_nm)
        half_diameter_index = int((height_index-lower_bound_index)/2 + lower_bound_index) #index of half diameter
        vert_slice_y_i =  np.transpose(tag_grid_t[ new_pos , : , ::-1])
        row_at_nm_check = vert_slice_y_i[int(nm_check/voxel_size)]
        row_vals = vert_slice_y_i[half_diameter_index,:]
        zero_index_row = np.nonzero(row_vals)
        zero_index_row_nm_check = np.nonzero(row_at_nm_check)
        if zero_index_row_nm_check[0].size == 0:
            pillar_widths_nm_check.append(0)
        else:
            pillar_widths_nm_check.append((np.max(zero_index_row_nm_check)-np.min(zero_index_row_nm_check))*voxel_size)
        if zero_index_row[0].size <= 3:
            pillar_HDW.append(0)
            continue
        pillar_HDW.append((np.max(zero_index_row)-np.min(zero_index_row))*voxel_size)
       
if deposit_type == 'wall':
    # Barebones --> Height/Material Sum Maps sufficient for stability study
    # These loops may be useful for someone to build out
    for i in range(pillar_num_x):
        for j in range(pillar_num_y):
            x_ij = int(midpoint+(-pillar_num_x+1)/2+i)
            y_ij = int(midpoint+(-pillar_num_y+1)/2+j)
            pillar_height_ij = h_map_t[y_ij,x_ij]
            pillar_heights.append(pillar_height_ij) #pillar heights
    

#Angle Determination - take the incline of the diagonal based on heights at first and last points
if deposit_type == 'line':
    pillar_sep = (pillar_num-1)*pitch
    height_sep = pillar_heights[-1]-pillar_heights[0]
    if pillar_sep != 0:
        incline = np.arctan(height_sep/pillar_sep)
    else:
        incline = 0
    incline_deg = round(incline*180/np.pi,2)
    incline = str(incline_deg)

#Writing certain values to a txt file

"""
f = open("Deposition_results.txt", "a")
f.write(title+"\n")
f.write("The average pillar height was "+str(round(pillar_height_av,2))+" nm\n")
f.write("The pillar height standard deviation was "+str(round(np.std(pillar_heights),3))+" nm\n")
f.write("The average HDW was "+str(round(pillar_HDW_av,2))+" nm\n")
f.write("The deposition efficiency was "+str(round(depo_eff,2))+"%"+"("+str(int(np.sum(profile)))+"/"+str(electron_num)+")\n")
f.write("The incline was :" + str(incline_deg)+"\n")

f.write("The Smith Comparison Deposition Amounts: PE:"+str(S_PE)+" FSE:"+str(S_FSE)+" BSE:"+str(S_BSE)+" SE1:"+str(S_SE_1)+" SE2:"+str(S_SE_2)+"\n")
f.write("Quick Values:\n" +str(round(pillar_height_av,2))+" "+str(round(pillar_HDW_av,2))+" "+ str(round(depo_eff,2))+"\n")
f.write("Quick Deposition Amounts:\n" +str(S_PE)+" "+str(S_FSE)+" "+str(S_BSE)+" "+ str(S_SE_1)+" "+str(S_SE_2)+"\n")
f.write("The pillar_Heights_were:\n"+str(pillar_heights)+"\n")
f.write("The pillar_HDW_were:\n"+str(pillar_HDW)+"\n")
f.write("The pillar_depth average is:\n"+str(round(np.average(pillar_depths),2))+"\n")
f.write("The pillar depth deviation is:\n"+str(round(np.std(pillar_depths),2))+"\n\n")

#Single Pillars or Lines
if len(pillar_widths_h1) >= 1:
    f.write("The average pillar width 1nm above substrate was "+str(round(pillar_widths_av,2))+" nm\n")
    f.write("The pillar width standard deviation 1nm above substrate was "+str(round(np.std(pillar_widths_av)))+" nm\n")
f.close()

"""

# Metrics to Strs to save
ph_max = str(round(np.max(pillar_heights),2))
ph_av = str(round(np.average(pillar_heights),2))
hdw_av = str(round(np.average(pillar_HDW),2))
hdw_std = str(round(np.std(pillar_HDW),2))
depth_av = str(round(np.average(pillar_depths),2))
depth_std = str(round(np.std(pillar_depths),2))
width_av =  str(round(np.average(pillar_widths_h1),2))
width_std = str(round(np.std(pillar_widths_h1),2))
fwhm_av = str(round(np.average(pillar_FWHM),2))
p_w_check_av =str(round(np.average(pillar_widths_nm_check),2))



# Write importnant metrics to file
f = open("deposition_results.txt", "a")
f.write(title+"\n")
f.write("Quick Values Order: ph_max ph_av hdw_av hdw_std depth_av dpth_std +width_av width_std fwhm_Av p_w_check_av incline\n")
f.write(ph_max+" "+ph_av+" "+hdw_av+" "+ hdw_std + " "+depth_av + " " + depth_std+ " " +width_av + " " + width_std + " " + fwhm_av +" "+p_w_check_av + " " + incline +"\n")
f.write("Quick Depositions Order: PE FSE BSE SE1 SE2 All_deps\n")
f.write(str(S_PE)+" "+str(S_FSE)+" "+str(S_BSE)+" "+ str(S_SE_1)+" "+str(S_SE_2)+" " +str(all_deps) +"\n")
f.write("\n")
f.close()

toc = time.perf_counter()
print(f"Finished  Plots in {toc - tic:0.4f} seconds")


#Pickling Geometry Data -------------------

print("Pickling Data for Cascades...")
tic = time.perf_counter()

# Debugging Tags for Cascade Tracking
tag_set = np.unique(tag_grid_t) #Primary Tags for Deposition
uncategorized_indices = np.where(new_species_grid_t==12) # Default for Uncategorized
debug_tags = tag_grid_t[uncategorized_indices]

# Unified Data Set
filename_geometry_data = path_to_test_files+"output"+date+parameter_summary+material+cross_section_type+"geometry_data.p"
geometry_data = [tag_set,height_map,debug_tags]
pickle.dump( geometry_data, open( filename_geometry_data, "wb" ) )

toc = time.perf_counter()
print(f"Finished  Pickling in {toc - tic:0.4f} seconds")

