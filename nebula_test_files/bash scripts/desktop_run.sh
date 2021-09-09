#!/bin/sh

#Shell script to run simulations using hpc 25 with access to nebula_vsc and nebula_test_files repositories.

#Needed Input Paths
geometry_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/vox_tri_pri/161_161_1001_sd_34000_sh_996_detect_dome_mirror_vs_250.bin"
#geometry_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/vox_tri_pri/81_81_101_sd_34000_sh_96_detect_dome_vs_250_RC_1keV_1_5kpp.bin"
#geometry_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/vox_tri_pri/401_401_1001_sd_34000_sh_996_mirrors_RC_1keV_21_1kpp_pitch_1700_seq_p.bin"
#geometry_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/vox_tri_pri/401_401_1001_sd_34000_sh_996_mirrors_RC_1keV_21_1kpp_pitchy_1700_pitchx_5500.bin"
#electron_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/vox_tri_pri/single_runs/3_8_2021/1keV_1_10kpp_pitch_0_401_401_101_sb_1000.pri"
#electron_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/vox_tri_pri/single_runs/4_8_2021/1keV_1_5kpp_pitch_0_401_401_101_sb_1000.pri"
electron_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/vox_tri_pri/single_runs/6_8_2021/1keV_1_10kpp_pitch_0_161_161_1001_sb_1000.pri"
#electron_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/vox_tri_pri/single_runs/12_7_2021/sem_simul_24_96_24_96_1.pri"
#electron_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/vox_tri_pri/single_runs/15_8_2021/1keV_rpit_1700_ppit_4000_1kpp_388p_pass_1_radii_30_10_cone_dep_random_401_401_101_sb_1000.pri"
#material_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/graphite-phonon.mat.hdf5"
#material_path_2="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/graphite-phonon_adsorbate.mat.hdf5"
material_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/tungsten.mat.hdf5"
material_path_2="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/tungsten_adsorbate.mat.hdf5"

#Executable Path
nebula_name="nebula_cpu_vsc"

#nebula_path geometry_path electron_path material_path
$nebula_name $geometry_path $electron_path $material_path $material_path_2
