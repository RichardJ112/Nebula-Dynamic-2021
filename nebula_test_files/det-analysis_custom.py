import sys
import os
import numpy as np
import pickle
import matplotlib.pyplot as plt

HPC_toggle = False
desktop_toggle = True

date = "/August/21_8_2021/"
material = 'tungsten'
cross_section_type = '_SM_'


#parameter_summary = "sem_simul_0_120_0_120_1_sd_34000_sh_996_mirrors_vs_300_RC_1keV_1_10kpp_"
#parameter_summary = "sem_simul_24_96_24_96_1_sd_34000_sh_996_mirrors_RCon_1keV_1__"
#parameter_summary = "sem_simul_24_96_24_96_1_sd_34000_sh_996_RC_sd_"
#parameter_summary = "sem_simul_36_84_36_84_1_sd_34000_sh_996_mirrors_RC_1keV_21_1_"
#parameter_summary = "sem_simul_24_96_24_96_1_sd_34000_sh_996_mirrors_RC_1keV_21_1_"
#parameter_summary = "sem_simul_0_120_0_120_1_sd_34000_sh_996_mirrors_vs_300_RC_no_shift_" #cone
#parameter_summary = "sem_simul_0_120_0_120_1_sd_34000_sh_996_mirrors_vs_300_RC_1keV_11_21_1kpp_SLL_" #wall
parameter_summary = "sem_simul_0_120_0_120_1_sd_34000_sh_996_mirrors_vs_300_RC_1keV_21_1kpp_pitch_2400_seq_p_"
#parameter_summary = "1keV_1_100kpp_pitch_0_401_401_1001_sb_1000_sd_34000_sh_996_detect_dome_RC_5keV_"
#parameter_summary = "5keV_1_100kpp_pitch_0_401_401_1001_sb_1000_sd_34000_sh_996_detect_dome_RC_5keV_"
#parameter_summary = "1keV_1_100kpp_pitch_0_401_401_1001_sb_1000_vs_300_sd_34000_sh_996_detect_dome_vs_300_"
#parameter_summary = "1keV_1_10kpp_pitch_0_401_401_1001_sb_1000_vs_300_sd_34000_sh_996_detect_dome_mirror_vs_300_"
#extra = "deposition_tungsten_yields/"
extra = ""

fignum = 1

#Input Paths
if HPC_toggle: #Adjustments for HPC environment
    #HPC path
    path_to_test_files = "/home/richarddejong/nebula_test_files/"
elif desktop_toggle:
    path_to_test_files = "C:/Users/Richard/source/repos/Nebula/nebula_test_files/"
else: 
    # Laptop path
    path_to_test_files= "C:/Users/richa/Documents/repos/nebula_test_files/"


input_path = path_to_test_files+"output"+date+extra+parameter_summary+material+cross_section_type+"detector.det"

path_list = (input_path.split('/')[-1]).split('.')
title =path_list[0] 
title_list = title.split('_')
electron_str = title_list[0]
#electron_energy = int(electron_str.replace('keV',''))*1000
electron_energy = 1000

output_path_folder = path_to_test_files+"figures"+date+title
general_path = output_path_folder+"/"+title

SEM_str = general_path +"_SEM_simul"
REELS_str_all_energy = general_path +"_REELS_simul_all_energy"
REELS_str_low_energy = general_path +"_REELS_simul_low_energy"
REELS_str_high_energy = general_path +"_REELS_simul_high_energy"
pickle_str = path_to_test_files +electron_str +"_energy_dist"

# This is a numpy datatype that corresponds to output files
electron_dtype = np.dtype([
    ('x',  '=f'), ('y',  '=f'), ('z',  '=f'), # Position
    ('dx', '=f'), ('dy', '=f'), ('dz', '=f'), # Direction
    ('E',  '=f'),                             # Energy
    ('px', '=i'), ('py', '=i')])              # Pixel index

# Open the output file
data = np.fromfile(input_path, dtype=electron_dtype)
print("Number of electrons detected: {}".format(len(data)))

if not os.path.exists(output_path_folder):
    os.makedirs(output_path_folder)

num_primary_electrons = 1e5
print("Number of electrons detected: {}".format(len(data)))

title_all = "NBD21 " + electron_str + " Electron Full Energy Distribution for Tungsten"
title_low = "NBD21 " + electron_str + " Electron Low Energy Distribution for Tungsten"
title_high = "NBD21 " + electron_str + " Electron High Energy Distribution for Tungsten"

low_energies = data[data['E'] < 50]
high_energies = data[data['E'] > 50]

N_SE = np.sum(data['E'] < 50)
N_BSE = np.sum(data['E'] > 50)
Yield_SE = round(N_SE/num_primary_electrons,2)
Yield_BSE = round(N_BSE/num_primary_electrons,2)

print("The SE Yield was "+str(Yield_SE))
print("The BSE Yield was "+str(Yield_BSE))


def energy_dist_plot(energies,title,save_str):
    max_energy = int(np.max(energies['E']))
    min_energy = int(np.min(energies['E']))
    N_bins = 4*max_energy
    spectrum, bin_edges = np.histogram(energies['E'], bins=N_bins)
    energy_intensity_max = bin_edges[np.argmax(spectrum)]

    plt.figure(fignum)
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
    plt.plot(bin_centers, spectrum, linewidth=2)
    #plt.title(title)
    plt.xlabel('Energy (eV)')

    if (title == title_high):
        x_ext = 0.025
        #plt.xlim(max_energy-20,max_energy+20)
        plt.xlim(electron_energy-20,electron_energy+20)
        plt.ylim(0,max(spectrum))


    plt.ylabel('Intensity (a.u.)')
    plt.savefig(save_str,bbox_inches = "tight")

plt.rcParams['font.size'] = 16
energy_dist_plot(data,title_all,REELS_str_all_energy)
fignum+=1
energy_dist_plot(low_energies,title_low,REELS_str_low_energy)
fignum+=1
energy_dist_plot(high_energies,title_high,REELS_str_high_energy)
fignum+=1

#pickle energies
pickle.dump(data,open( pickle_str, "wb" ))

# Make a histogram of pixel indices
data = low_energies #only use secondaries
xmin = data['px'].min()
xmax = data['px'].max()
ymin = data['py'].min()
ymax = data['py'].max()
H, xedges, yedges = np.histogram2d(data['px'], data['py'],
	bins = [
		np.linspace(xmin-.5, xmax+.5, xmax-xmin+2),
		np.linspace(ymin-.5, ymax+.5, ymax-ymin+2)
	])

pixel_size = 0.5
xmin_nm = -xmin*pixel_size
xmax_nm = xmax*pixel_size
ymin_nm = -ymin*pixel_size
ymax_nm = ymax*pixel_size

range_x = xmax_nm-xmin_nm
range_y = ymax_nm-ymin_nm

xm = range_x/2
ym = range_y/2

# Make a plot
#f_range = int(title_list[3])-int(title_list[2])
plt.rcParams['font.size'] = 16
plt.figure(fignum)
plt.imshow(H.T, cmap='gray', vmin=0,extent=(-xm,xm,-ym,ym))
plt.colorbar()
plt.xlim(-30,30)
plt.ylim(-30,30)
plt.xlabel('x (nm)')
plt.ylabel('y (nm)')
plt.savefig(SEM_str,bbox_inches = "tight")
#plt.savefig("Test")

print("Plots Done")