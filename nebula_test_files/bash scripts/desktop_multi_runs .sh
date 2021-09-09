#!/bin/sh

#Shell script to run simulations using desktop with access to nebula_vsc and nebula_test_files repositories.

#Needed Input Paths
geometry_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/vox_tri_pri/401_401_35_sd_1000_sh_30_mirrors.bin"
electron_path_folder="/home/richarddejong/nebula_test_files/vox_tri_pri/multi_runs/7_7_21/*"
material_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/tungsten.mat.hdf5"
material_path_2="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/tungsten_adsorbate.mat.hdf5"

#Executable Path
nebula_name="nebula_cpu_vsc"


declare -i num=1
clear
cd ~
#nebula_path geometry_path electron_path material_path
for file in $electron_path_folder
do
  echo "Processing $file file... $num"
  $nebula_name $geometry_path $file $material_path $material_path_2
  num+=1
done

