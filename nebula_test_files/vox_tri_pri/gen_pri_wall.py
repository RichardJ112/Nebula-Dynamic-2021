import numpy as np
import os

# Point exposure for Wall Depositions

#Automation Additions
HPC_toggle = False
multi_run_folder_toggle = False# single or multi folder (True is multi,False is single)
desktop_toggle = True

date = "/20_8_2021/"

# Parameters Geom
voxel_size = 0.3; # voxel size in nanometers (0.27 nm is appr. one atom of Si)(Smith uses 0.25nm for Tungsten(W))
size_x= 401 # horizontal size in the x direction in voxels (now for +/- x)
size_y = 401 # horizontal size in the y direction in voxels (now for +/- y)
size_z = 1001 # vertical size in voxels
volume = size_x*size_y*size_z # total voxel volume
voxel_size_pm = int(voxel_size*1000)

# Parameters Pri
num_x = int(11) # number of pillars in the x direction
num_y = int(21) #number of pillars in the y direction
num_x_range = np.linspace(-(num_x-1)/2,(num_x-1)/2,num_x)
num_y_range = np.linspace(-(num_y-1)/2,(num_y-1)/2,num_y)
x = size_x*voxel_size/2       # starting x
y = size_y*voxel_size/2        # starting y
z = 1      # starting z in nm
N = int(1e3)    # Number of electrons per pillar
seq_lines = 100
energy = 1000 # Beam energy, in eV
sigma_beam = 1 # Beam standard deviation in nm
sigma_beam_pm = round(sigma_beam*1000) # Beam standard deviation in pm
tot_e = int(num_x*num_y*N)


#line_pitch_x_list = [6,6.5,7] 
#line_pitch_x_list = [3,3.5,4] 
line_pitch_x_list = [4.5] 
#line_pitch_x_list = np.around(np.linspace(3,4,num=10),decimals=1) 

#line_pitch_y_list = [1.7,1.7,1.7] #space between points nm
line_pitch_y_list = np.ones(len(line_pitch_x_list))*1.7

#line_pitch_list = np.around(np.linspace(1,5,num=20),decimals=1) #space between points nm
dep_strat = "seq_ll" #seq_ll (sequential layered lines), seq_d (sequential diagonals) 

# This is a numpy datatype that corresponds to pri files
dt = np.dtype([
	('x', np.float32), ('y', np.float32), ('z', np.float32),    # Starting position
	('dx', np.float32), ('dy', np.float32), ('dz', np.float32), # Starting direction
	('K', np.float32),                                          # Starting energy
	('px', np.uint32), ('py', np.uint32)])                      # Pixel index

total_runs = len(line_pitch_y_list)
# For now only iterates over line_pitches --> may update to include any possible combination of parameters 

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

if dep_strat == "seq_d":
    for i in range(total_runs):
        # Update Iterative Parameters
        line_pitch_x = line_pitch_x_list[i]
        line_pitch_x_pm = int(line_pitch_x*1000)
        line_pitch_y = line_pitch_y_list[i]
        line_pitch_y_pm = int(line_pitch_y*1000) 

        #title creation
        title_pri = str(int(energy/1000))+"keV_"+str(num_x)+"_"+str(num_y)+"_"+str(int(N/1000))+"kpp_pitchx_"+str(line_pitch_x_pm)+"_pitchy_"+str(line_pitch_y_pm)+"_"+dep_strat+"_"
        title_geom = str(size_x)+"_"+str(size_y)+"_"+str(size_z)+"_sb_"+str(sigma_beam_pm)+"_vs_"+str(voxel_size_pm)
        title = title_pri+title_geom+".pri"
        
        x_p = []
        y_p = []
        inputx = np.zeros(tot_e)
        inputy = np.zeros(tot_e)

        for k in num_x_range:
            for j in num_y_range:
                startx = k*line_pitch_x+x
                starty = j*line_pitch_y+y
                xj = np.random.normal(startx, sigma_beam, N)
                yj = np.random.normal(starty, sigma_beam, N)
                x_p.append(xj)
                y_p.append(yj)

        for j in range(num_y*num_x):
                inputx[j*N:N*(j+1)] = x_p[j]
                inputy[j*N:N*(j+1)] = y_p[j]

        # Open file

        if not os.path.exists(file_path):
            os.makedirs(file_path)

        with open(file_path+title, 'wb') as file:
            # Allocate numpy buffer
            array = np.empty(tot_e, dtype=dt)

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

if dep_strat == "seq_ll":
    for i in range(total_runs):
        # Update Iterative Parameters
        line_pitch_x = line_pitch_x_list[i]
        line_pitch_x_pm = int(line_pitch_x*1000)
        line_pitch_y = line_pitch_y_list[i]
        line_pitch_y_pm = int(line_pitch_y*1000) 

        #title creation
        title_pri = str(int(energy/1000))+"keV_"+str(num_x)+"_"+str(num_y)+"_"+str(int(N/1000))+"kpp_pitchx_"+str(line_pitch_x_pm)+"_pitchy_"+str(line_pitch_y_pm)+"_"+dep_strat+"_"
        title_geom = str(size_x)+"_"+str(size_y)+"_"+str(size_z)+"_sb_"+str(sigma_beam_pm)+"_vs_"+str(voxel_size_pm)
        title = title_pri+title_geom+".pri"

        x_p = []
        y_p = []
        inputx = np.zeros(tot_e)
        inputy = np.zeros(tot_e)
        N_seq = int(N/seq_lines)
        num_line = num_y*seq_lines
        num_sum = 0
        for k in num_x_range: #We switch x and y
            starty = k*line_pitch_y+y
            for i in range(seq_lines): 
                for j in num_y_range:
                    startx = j*line_pitch_x+x
                    xj = np.random.normal(startx, sigma_beam, N_seq)
                    yj = np.random.normal(starty, sigma_beam, N_seq)
                    x_p.append(xj)
                    y_p.append(yj)
            for j in range(num_line):
                inputx[(j+num_sum)*N_seq:N_seq*(j+1+num_sum)] = x_p[j+num_sum]
                inputy[(j+num_sum)*N_seq:N_seq*(j+1+num_sum)] = y_p[j+num_sum]
            num_sum+=num_line
        # Open file

        if not os.path.exists(file_path):
            os.makedirs(file_path)

        with open(file_path+title, 'wb') as file:
            # Allocate numpy buffer
            array = np.empty(num_y*N*num_x, dtype=dt)

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
