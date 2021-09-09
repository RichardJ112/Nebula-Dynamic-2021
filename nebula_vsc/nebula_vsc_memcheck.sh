#!/bin/sh

#Shell script to build nebula_vsc

#Needed Input Paths
geometry_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/vox_tri_pri/161_161_1001_sd_34000_sh_996_detect_dome_mirror_vs_250.bin"
electron_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/vox_tri_pri/single_runs/6_8_2021/1keV_1_1kpp_pitch_0_161_161_1001_sb_1000.pri"
material_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/tungsten.mat.hdf5"
material_path_2="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/tungsten_adsorbate.mat.hdf5"

test_name="memcheck_test_log.txt" 


valgrind --leak-check=full --show-leak-kinds=all --log-file=$test_name nebula_cpu_vsc $geometry_path $electron_path $material_path $material_path_2
