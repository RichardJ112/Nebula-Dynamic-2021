import numpy as np
import os

# Point exposure

#This version was adjusted for hpc25 ssh server - 10_05_2021

#Automation Additions
HPC_toggle = False
multi_run_folder_toggle = False # single or multi folder (True is multi,False is single)
desktop_toggle = True

date = "/2_9_2021/"

# Parameters Geom
voxel_size = 0.3; # voxel size in nanometers (0.27 nm is appr. one atom of Si)(Smith uses 0.25nm for Tungsten(W))(0.96 silicon movement)
size_x= 251 # horizontal size in the x direction in voxels (now for +/- x)
size_y = 251 # horizontal size in the y direction in voxels (now for +/- y)
size_z = 1001 # vertical size in voxels
volume = size_x*size_y*size_z # total voxel volume
voxel_size_pm = int(voxel_size*1000)

# Parameters Pri
num = 1 # number of pillars(odd)
numrange = np.linspace(-(num-1)/2,(num-1)/2,num)
x = size_x*voxel_size/2       # starting x
y = size_y*voxel_size/2        # starting y
z = 1      # starting z in nm
N = int(200e3)    # Number of electrons per pillar
energy = int(1e3) # Beam energy, in eV
sigma_beam = 1 # Beam standard deviation in nm
sigma_beam_pm = round(sigma_beam*1000) # Beam standard deviation in pm
line_pitch_list = [0] #space between points nm
#line_pitch_list = np.around(np.linspace(1,4,num=40),decimals=1) #space between points nm

# This is a numpy datatype that corresponds to pri files
dt = np.dtype([
	('x', np.float32), ('y', np.float32), ('z', np.float32),    # Starting position
	('dx', np.float32), ('dy', np.float32), ('dz', np.float32), # Starting direction
	('K', np.float32),                                          # Starting energy
	('px', np.uint32), ('py', np.uint32)])                      # Pixel index

total_runs = len(line_pitch_list)
# For now only iterates over line_pitches --> may update to include any possible combination of parameters
for i in range(total_runs):

	#Update iterative parameters
	line_pitch = line_pitch_list[i]
	line_pitch_pm = int(line_pitch_list[i]*1000) #line_pitch_pm this is since 22 March

	electron_num_str = str(int(N/1000))
	if N < 1000:
		electron_num_str = '0-'+str(N)
	#title creation
	title_pri = str(int(energy/1000))+"keV_"+str(num)+"_"+electron_num_str+"kpp_pitch_"+str(line_pitch_pm)+"_"
	title_geom = str(size_x)+"_"+str(size_y)+"_"+str(size_z)+"_sb_"+str(sigma_beam_pm)+"_vs_"+str(voxel_size_pm)+".pri"
	title = title_pri+title_geom

	x_p = []
	y_p = []
	inputx = np.zeros([num*N])
	inputy = np.zeros([num*N])

	for j in numrange:
		starty = j*line_pitch+y
		xj = np.random.normal(x, sigma_beam, N)
		yj = np.random.normal(starty, sigma_beam, N)
		x_p.append(xj)
		y_p.append(yj)

	for j in range(num):
		inputx[j*N:N*(j+1)] = x_p[j]
		inputy[j*N:N*(j+1)] = y_p[j]

	# Open file

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


	if not os.path.exists(file_path):
		os.makedirs(file_path)

	with open(file_path+title, 'wb') as file:
		# Allocate numpy buffer
		array = np.empty(num*N, dtype=dt)

		# Fill with data
		array['x'] = inputx
		array['y'] = inputy
		array['z'] = z
		array['dx'] = 0
		array['dy'] = 0
		array['dz'] = 1
		array['K'] = energy
		array['px'] = 0
		array['py'] = 0

		# Write buffer to file
		array.tofile(file)

	#Progress Tracker
	if i%10 == 0:
		print("Creating Primary Files :" + str(round(i/total_runs*100)) + "%")
