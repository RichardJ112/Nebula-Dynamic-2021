#!/bin/sh

#Shell script to run desktop

cd ..
cd nebula_test_files

#Shell script to run simulations using desktop with access to nebula_vsc and nebula_test_files repositories.

#Needed Input Paths
geometry_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_dynamic_2021/nebula_test_files/vox_tri_pri/401_401_1001_sd_34000_sh_996_detect_dome_vs_300.bin"

electron_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_dynamic_2021/nebula_test_files/vox_tri_pri/single_runs/28_10_2021/1keV_1_10kpp_pitch_0_401_401_1001_sb_1000_vs_300.pri"

material_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_dynamic_2021/nebula_test_files/material_files/silicon.mat.hdf5"
material_path_2="/mnt/c/Users/Richard/source/repos/Nebula/nebula_dynamic_2021/nebula_test_files/material_files/silicon_adsorbate.mat.hdf5"

#nebula_path geometry_path electron_path material_paths
nebula_cpu_vsc $geometry_path $electron_path $material_path $material_path_2
