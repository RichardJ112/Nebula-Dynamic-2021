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
import time
from struct import *

# This file is for cone depositions

# IO parameters
date = "/August/19_8_2021/"
material = 'tungsten'
#material = 'graphite-phonon'
cross_section_type = '_SM_'
parameter_summary = "1keV_rpit_1700_ppit_4000_1kpp_425p_pass_1_radii_30_4_cone_dep_4_division_shift_401_401_1001_sb_1000_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_"
#parameter_summary = "1keV_rpit_1700_ppit_4000_1kpp_425p_pass_1_radii_30_4_cone_dep_random_401_401_1001_sb_1000_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_"
deposit_type = "cone"
HPC_toggle = False
desktop_toggle = True
enable_e_grid = False
pillar_num = 1

#Plotting Toggles
surface_map_3d_toggle = True

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

path_list = (input_path.split('/')[-1]).split('.')
title =path_list[0]
title_list = title.split('_')
electrons_per_point = int(title_list[5].replace('kpp',''))*1000
r_pitch_nm = int(title_list[2])/1000
p_pitch_nm = int(title_list[4])/1000
passes = int(title_list[8])
point_num = int(title_list[6].replace('p',''))

electron_num = point_num * electrons_per_point

output_path_folder = path_to_test_files+"figures"+date+title
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
cs_hor_surface0_str_new = general_path +"_cs_hor_surface0_new"# horizontal cross section near surface
t_evo_surface = general_path +"_surface_time_evolution"# horizontal cross section near surface
t_evo_halo = general_path +"_surface_halo_time_evolution"# horizontal cross section at z = 0nm
cs_vert_x0_str_new = general_path +"_cs_vert_X0_new" # vertical cross section at X=0(New Classification)
cs_vert_y0_str_new = general_path +"_cs_vert_Y0_new" # vertical cross section at X=0(New Classification)
cs_vert_x_str_new = general_path +"_cs_vert_X_new" # vertical cross section at specific X (New Classification)
rad_dist_str = general_path + "_rad_dist" # radial distribution for pillar
lat_dist_str = general_path + "_lat_dist" # lateral distribution for lines
height_map_str = general_path + "_height_map"

#Experimental path strings
height_map_3d_str = general_path + "_3d_height_map"
surf_map_3d_str = general_path + "_3d_surf_map"

#Text Loading
print("Loading file...")

tic =  time.perf_counter() 
###Change this to allow for bin outputs

with open(input_path, mode='rb') as file: # b is important -> binary
   fileContent = file.read()

#Binary Unpacking
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

#Fixed Width Ranges --------------------------------------------------------
f_range = 40#fixed range in nm
f_x = np.arange(start =-f_range,stop = f_range,step =0.3)
f_y = np.arange(start = -f_range,stop =f_range,step = 0.3)
f_steps = len(f_x)/2
s_dims = [dim[0]/2,dim[1]/2]
f_x_min = int(s_dims[0]-f_steps)
f_x_max = int(s_dims[0]+f_steps)
f_y_min = int(s_dims[1]-f_steps)
f_y_max = int(s_dims[1]+f_steps)
mat_grid_f = mat_grid_t[f_x_min:f_x_max,f_y_min:f_y_max,:]
tag_grid_f =  tag_grid_t[f_x_min:f_x_max,f_y_min:f_y_max,:]
new_species_grid_f = new_species_grid_t[f_x_min:f_x_max,f_y_min:f_y_max,:]

f_x = np.arange(start =-f_range,stop = f_range,step =0.3)
f_y = np.arange(start = -f_range,stop =f_range,step = 0.3)

toc = time.perf_counter()
print(f"Finished Data Manipulation in {toc - tic:0.4f} seconds")

print("Creating Plots...")

tic = time.perf_counter()


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

#Lower Bound
def find_lower_bound_index(pillar_vals,zero_index):
    zero_indices = zero_index[0]
    for i in range(len(zero_indices)-1):
        if zero_indices[i+1] != zero_indices[i]+1:
            if pillar_vals[zero_indices[i+1]+1] != 0:
                lsi = int(zero_indices[i+1])
                return lsi
    return 0

#Height Map
h_map_f = height_map(tag_grid_f,axis=2)
plt.figure(fignum)
fignum += 1
plt.imshow(h_map_f, extent=(-f_range, f_range, -f_range, f_range), origin={'upper', 'upper'})
plt.xlabel("x (nm)")
plt.ylabel("y (nm)")
#plt.title(r"Height Map", pad= 15)
plt.colorbar(shrink = 1,label = "Height(nm)")
plt.savefig(height_map_str, aspect = 'auto', bbox_inches="tight")

#3D Surface Map ---------------------------------------------------------------------------------

if surface_map_3d_toggle:
    plt.figure(fignum)
    fignum += 1
    plt.rcParams['font.size'] = 16
    surface_list = surface_index(tag_grid_t,dim,remove_surface=1,remove_range=f_range)
    ax = plt.axes(projection='3d')
    #plt.title("Deposit 3D Surface Scan", pad= 15)
    ax.view_init(elev=30, azim=60)
    ax.scatter(surface_list[0], surface_list[1], surface_list[2],c=surface_list[2], cmap='viridis', linewidth=0.5,alpha=0.1)
    #ax.plot_trisurf(surface_list[0], surface_list[1], surface_list[2], cmap='viridis', linewidth=0.3)
    #plt.show()
    ax.set_xlim(-f_range,f_range)
    ax.set_ylim(-f_range,f_range)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_zlabel("z (nm)")
    plt.xticks(np.arange(-f_range, f_range+1, 20))
    plt.yticks(np.arange(-f_range, f_range+1, 20))
    plt.tight_layout()
    plt.savefig(surf_map_3d_str, aspect = 'auto',bbox_inches="tight")

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

cmap = LinearSegmentedColormap.from_list('mycmap', ['white', 'red', 'blue', 'cyan', 'green', 'xkcd:olive brown'])
cmap8 = LinearSegmentedColormap.from_list('mycmap', ['white', 'gray', 'red', 'green', 'blue', 'purple', '#FFA020', '#80FF80', 'cyan', 'magenta']) #10 colours for 8 types
cmap12 = LinearSegmentedColormap.from_list('mycmap', ['white', 'gray', 'red', 'green', 'blue', 'purple', '#FFA020', '#80FF80', 'cyan', 'magenta','orange','black','#301713']) #13 colors for 13 types

#ColorMap comparision
"""
    Old:(colours for cmap_comp)
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
cmap_comp = LinearSegmentedColormap.from_list('mycmap', ['white', 'red', 'green','blue','blue','#dc5858','blue','blue','#ffa020','#80ff80','magenta','purple','gray']) #13 colors for 13 types

#Horizontal cross section at surface (new classification)-----------------------------------------------

ul_plot = plt.figure(fignum)
fignum+=1
plt.imshow(new_species_grid_f[:, :, -1], cmap=cmap_comp, extent=(-f_range, f_range, -f_range, f_range),interpolation = "nearest",vmin= 0,vmax = 12)
plt.xlabel("x (nm)")
plt.ylabel("y (nm)")
#plt.title(" Horizontal cross section at surface", pad = 8)
ul_plot.savefig(cs_hor_surface0_str_new, dpi = 200, bbox_inches="tight")

"""
# Time evolution graph at X = 0 -------------------------------------------------------
tag_cs_vert_x = np.transpose(tag_grid_f[ : , int((f_x_max-f_x_min)/2), :])

time_evo_plot = plt.figure(fignum)
fignum+=1

plt.imshow(tag_cs_vert_x, extent=(-f_range, f_range,0, dim[2] * voxel_size))
plt.xlabel("y (nm)")
plt.ylabel("z (nm)")
plt.title("Time evolution at X = 0")
time_evo_plot.savefig(t_evo_x_str, dpi = 200, bbox_inches="tight")
"""

# Time evolution graph at a specific X Contours*-------------------------------------------------------
tag_cs_vert_x = np.transpose(tag_grid_f[ : , int((f_x_max-f_x_min)/2), :])
time_evo_plot_X_C = plt.figure(fignum)
fignum+=1 
plt.rcParams['font.size'] = 16
f_z = np.arange(start =0,stop = (dim[2])*voxel_size,step =voxel_size)
contour_Z,contour_Y = np.meshgrid(f_y,f_z)

electron_increment = int(electron_num/10) 
contour_levels = np.arange(start=0, stop=electron_num+1, step=electron_increment)
plt.contour(contour_Z,contour_Y,np.flipud(tag_cs_vert_x),levels = contour_levels,colors = 'black',linewidths = 0.3)
plt.contourf(contour_Z,contour_Y,np.flipud(tag_cs_vert_x),levels = contour_levels)
plt.colorbar(shrink = 1,label = "Primary Electron ID")
plt.xlabel("y (nm)")
plt.ylabel("z (nm)")
#plt.title("Time evolution at X = 0 Section Countour Plot")
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

"""
# Time evolution graph at a specific Y -------------------------------------------------------
tag_cs_vert_y = np.transpose(tag_grid_f[ int((f_x_max-f_x_min)/2) , : , :])

time_evo_plot = plt.figure(fignum)
fignum+=1
plt.imshow(tag_cs_vert_y,extent=(-f_range, f_range,0, dim[2] * voxel_size))
plt.xlabel("x (nm)")
plt.ylabel("z (nm)")
plt.title("Time evolution at Y = 0")
time_evo_plot.savefig(t_evo_y_str, dpi = 200, bbox_inches="tight")

# Time evolution graph at a specific Z -------------------------------------------------------

ul_plot = plt.figure(fignum)
fignum+=1
plt.imshow(tag_grid_f[:, :, -3], extent=(-f_range, f_range, -f_range, f_range))
plt.xlabel("x (nm)")
plt.ylabel("y (nm)")
plt.title(" Time evolution at Z = 1nm above surface", pad = 8)
ul_plot.savefig(t_evo_surface, dpi = 200,bbox_inches="tight")

"""
# Cross Sections Electron Classification -------------------------------------------------------------------------------

# Vertical cross section at X=0 (New Classification)-----------------------------------------
class_legend_toggle = True
plt.rcParams['font.size'] = 16
cs_plot = plt.figure(fignum)
fignum+=1
vert_slice_x_f_new =  np.transpose(new_species_grid_f[ : , int((f_x_max-f_x_min)/2), ::-1])
im = plt.imshow(vert_slice_x_f_new,cmap = cmap_comp,extent=(-f_range, f_range, 0, dim[2] * voxel_size),origin='lower',interpolation = "nearest",vmin= 0,vmax = 12);

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
plt.xlabel("x (nm)")
plt.ylabel("z (nm)")
#plt.title("Vertical cross section through center at Y = 0")
cs_plot.savefig(cs_vert_y0_str_new, dpi = 200, bbox_inches="tight")

toc = time.perf_counter()
print(f"Finished  Plots in {toc - tic:0.4f} seconds")



