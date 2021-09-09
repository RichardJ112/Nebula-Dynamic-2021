#!/bin/sh

#Shell script to run simulations using hpc 25 with access to nebula_vsc and nebula_test_files repositories.

#Needed Input Paths
geometry_path="/home/richarddejong/nebula_test_files/vox_tri_pri/401_401_801_sd_4000.bin"
electron_path="/home/richarddejong/nebula_test_files/vox_tri_pri/single_runs/17_5_21/10keV_21_50kpp_pitch_4500_401_401_801_sb_1000.pri"
material_path="/home/richarddejong/nebula_test_files/graphite-phonon.mat.hdf5"
#material_path="/home/richarddejong/nebula_test_files/tungsten.mat.hdf5"
#Executable Path
nebula_name="/home/richarddejong/nebula_vsc/build/bin/nebula_cpu_vsc"

clear
cd ~
#nebula_path geometry_path electron_path material_path
$nebula_name $geometry_path $electron_path $material_path
