#!/bin/sh

#Shell script to run simulations using hpc 25 with access to nebula_vsc and nebula_test_files repositories.

#Needed Input Paths
geometry_path="/mnt/c/Users/richa/Documents/repos/nebula_test_files/vox_tri_pri/401_401_35_sd_1000_sh_30_mirrors.bin"
electron_path="/mnt/c/Users/richa/Documents/repos/nebula_test_files/vox_tri_pri/single_runs/7_7_2021/1keV_1_100kpp_pitch_0_401_401_35_sb_0.pri"
material_path="/mnt/c/Users/richa/Documents/repos/nebula_test_files/tungsten.mat.hdf5"
material_path_2="/mnt/c/Users/richa/Documents/repos/nebula_test_files/tungsten_adsorbate.mat.hdf5"
#material_path="/home/richarddejong/nebula_test_files/tungsten.mat.hdf5"
#Executable Path
nebula_name="/mnt/c/Users/richa/Documents/repos/nebula_vsc/build/bin/nebula_cpu_vsc"

#nebula_path geometry_path electron_path material_path
$nebula_name $geometry_path $electron_path $material_path $material_path_2
