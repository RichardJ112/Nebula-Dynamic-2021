import numpy as np
import sys
import os
import matplotlib
import math
import os
import matplotlib.pyplot as plt
#This file is for cone depositions

#Automation Additions
HPC_toggle = False
multi_run_folder_toggle = False # single or multi folder (True is multi,False is single)
desktop_toggle = True

date = "/27_8_2021/"
month = "/August/"

# Parameters Geom
voxel_size = 0.3; # voxel size in nanometers (0.27 nm is appr. one atom of Si)
voxel_size_pm = int(voxel_size*1000)

size_x= 401 # horizontal size in the x direction in voxels (now for +/- x)
size_y = 401 # horizontal size in the y direction in voxels (now for +/- y)
size_z = 1001 # vertical size in voxels
volume = size_x*size_y*size_z # total voxel volume
sim_depth = 4000 # simulation depth under the voxels at z < 0 for SEM bulk samples in nm
sample_height = 4000 # height of the sample (length between the sample and the top of the simulation domain, in vacuum) in voxels


#Cone Parameters
r_pitch_nm = 1.7 #radius pitch(spacing towards inner part of cone)
r_bot_nm = 30 #bottom radius target in nm
r_top_nm = 4 #top radius target in nm
p_pitch_nm = 4 #pitch between points
r_pitch_pm = round(r_pitch_nm*1000)
p_pitch_pm = round(p_pitch_nm*1000)
N = int(1000)   # Number of electrons per pillar
passes = 1
plot_preview = True
dep_alt = False
dep_alt_str = ""
if dep_alt:
	dep_alt_str="_alt"
divisions = 4



#Changes Deposition Type for Cones - 'random','pi_shift' and 'no_shift'
cone_dep_type = 'random'
#cone_dep_type = str(divisions)+'_division_shift'
cone_dep_str = cone_dep_type+dep_alt_str

# Parameters Pri:
x_0 = size_x*voxel_size/2       # starting x
y_0 = size_y*voxel_size/2        # starting y
z_0 = 1      # starting z

energy = 1000 # Beam energy, in eV
sigma_beam = 1 # Beam standard deviation in nm'
sigma_beam_pm = round(sigma_beam*1000) # Beam standard deviation in pm

#Point Map

pi = math.pi

# This is a numpy datatype that corresponds to pri files
dt = np.dtype([
	('x', np.float32), ('y', np.float32), ('z', np.float32),    # Starting position
	('dx', np.float32), ('dy', np.float32), ('dz', np.float32), # Starting direction
	('K', np.float32),                                          # Starting energy
	('px', np.uint32), ('py', np.uint32)])         	           # Pixel index

#Point Map Implementation - Generate Deposition Positions
x_p = []
y_p = []

n_c = round(abs(r_top_nm-r_bot_nm)/r_pitch_nm) #number of circles
n_c_steps=round(2*pi*r_bot_nm/p_pitch_nm) # steps large circle
phi = 0
phi_inc = 0

if (cone_dep_type == 'pi_shift'):
	phi_inc = pi

if (cone_dep_type == str(divisions)+"_division_shift"):
	phi_inc = 2*pi/divisions
	phi_inc_alt = 2*pi/divisions

alt_phi = 0
#Generate All Points
n_points_circ  = n_c_steps
rad_circ = r_bot_nm #circle radius
for i in range(n_c):
	theta = np.linspace(0, 2*pi,num = n_points_circ)
	theta = theta + phi
	x_circ_points = rad_circ * np.cos(theta)+x_0
	y_circ_points = rad_circ * np.sin(theta)+y_0
	for j in range(passes):
		x_p.append(x_circ_points)
		y_p.append(y_circ_points)
	#Change parameters for next round
	rad_circ -= r_pitch_nm
	n_points_circ = round(2*pi*rad_circ/p_pitch_nm)
	if (cone_dep_type == 'random'):
		phi_inc = pi*np.random.normal(0,1,1)
	#Alternating
	if dep_alt:
		phi = -phi
		if i%2:
			phi = 0
			phi_inc += phi_inc_alt
			phi += phi_inc
	else:
		phi += phi_inc

#Array Flattening
x_p_array = np.hstack(x_p)
y_p_array = np.hstack(y_p)
x_p_flat = x_p_array.flatten()
y_p_flat = y_p_array.flatten()


#Generate Gaussian Deposition of N electrons at deposition positions
n_points = len(x_p_flat)
x_list = []
y_list = []
for i in range(n_points):
	xi = np.random.normal(x_p_flat[i], sigma_beam, N)
	yi = np.random.normal(y_p_flat[i], sigma_beam, N)
	x_list.append(xi)
	y_list.append(yi)

input_x = np.hstack(x_list).flatten()
input_y = np.hstack(y_list).flatten()
electron_array_size = len(input_x)

#title creation -new naming
# Lines - energy #points #e_pp pitch dimensions beam_diameter
# Cones - energy r_pitch pt_pitch #e_pp #points passes dimensions_beam_diameter_radius_bottom_radius_top type
title_pri = str(int(energy/1000))+"keV_rpit_"+str(r_pitch_pm)+"_ppit_"+str(p_pitch_pm)+"_"+str(int(N/1000))+"kpp_"+str(n_points)+"p_pass_"+str(passes)+"_radii_"+str(r_bot_nm)+"_"+str(r_top_nm)+"_cone_dep_"+cone_dep_str+"_"
title_geom = str(size_x)+"_"+str(size_y)+"_"+str(size_z)+"_sb_"+str(sigma_beam_pm)+"_vs_"+str(voxel_size_pm)
title = title_pri+title_geom+".pri"
title_plot = title_pri + title_geom

#Input Paths
if HPC_toggle: #Adjustments for HPC environment
    #HPC path
    path_to_test_files = "/home/richarddejong/nebula_test_files/"
elif desktop_toggle:
    path_to_test_files = "C:/Users/Richard/source/repos/Nebula/nebula_test_files/"
else: 
    # Laptop path
    path_to_test_files= "C:/Users/richa/Documents/repos/nebula_test_files/"

#Folder Location
if multi_run_folder_toggle:
    #Multi-Runs - Run all files in folder (for multi-runs exectuable)
    if HPC_toggle:
        file_path = "/home/richarddejong/nebula_test_files/vox_tri_pri/multi_runs"+date 
    elif desktop_toggle:
        file_path = "C:/Users/Richard/source/repos/Nebula/nebula_test_files/vox_tri_pri/multi_runs"+date
    else:
        file_path = "C:/Users/richa/Documents/repos/nebula_test_files/vox_tri_pri/multi_runs"+date
else:
    #Single Runs - Run a single file (for single-run executable)
    if HPC_toggle:
        file_path = "/home/richarddejong/nebula_test_files/vox_tri_pri/single_runs"+date 
    elif desktop_toggle:
        file_path = "C:/Users/Richard/source/repos/Nebula/nebula_test_files/vox_tri_pri/single_runs"+date
    else:
        file_path = "C:/Users/richa/Documents/repos/nebula_test_files/vox_tri_pri/single_runs"+date


output_path_folder = path_to_test_files+"figures"+month+date+title_plot

dep_strategy_str = output_path_folder+"_dep_strategy" #top down profile of growth

cmap = matplotlib.cm.get_cmap('viridis')

col_points = cmap(np.linspace(0,1,x_p_flat.size))

#Generate Plot for All deposition areas
plt.rcParams['font.size'] = 16
if plot_preview:
	fignum = 0
	cone_test_plot = plt.figure()
	fignum+=1

	for i in range(len(x_p)):
		x_c = x_p[i]
		y_c = y_p[i]
		col_points = cmap(np.linspace(0,1,x_c.size))
		plt.scatter(x_c,y_c, c=col_points)
	dep_title = "Radial Limits: ("+str(r_bot_nm) +","+str(r_top_nm)+")"+"Radial,Point Pitch: ("+str(r_pitch_nm)+ "," +str(p_pitch_nm) +")" + " Deposition Strategy: "+ cone_dep_str
	plt.colorbar(shrink = 1,label = "Concentric Deposition Order")
	#plt.title(dep_title, pad= 15)
	plt.xlabel("x (nm)")
	plt.ylabel("y (nm)")
	plt.savefig(dep_strategy_str, aspect = 'auto', bbox_inches="tight")
	#plt.show()

# Primary Count Warning
print("\n Warning: This file contains "+str(electron_array_size) + " electrons")

if not os.path.exists(file_path):
	os.makedirs(file_path)

# Open file	
with open(file_path+title, 'wb') as file:
	# Allocate numpy buffer
	array = np.empty(electron_array_size, dtype=dt)

	# Fill with data
	array['x'] = input_x
	array['y'] = input_y
	array['z'] = z_0
	array['dx'] = 0
	array['dy'] = 0
	array['dz'] = 1
	array['K'] = energy
	array['px'] = 0
	array['py'] = 0

	# Write buffer to file
	array.tofile(file)

