#!/bin/sh

#Shell script to run simulations using hpc 25 with access to nebula_vsc and nebula_test_files repositories.

#Needed Input Paths
geometry_path="/home/richarddejong/nebula_test_files/vox_tri_pri/401_401_1001_sd_4000.bin"
electron_path_folder="/home/richarddejong/nebula_test_files/vox_tri_pri/multi_runs/10_05_21/*"
material_path="/home/richarddejong/nebula_test_files/graphite-phonon.mat.hdf5"
declare -i num=1
#Executable Path
nebula_name="/home/richarddejong/nebula_vsc/build/bin/nebula_cpu_vsc"

clear
cd ~
#nebula_path geometry_path electron_path material_path
for file in $electron_path_folder
do
  echo "Processing $file file... $num"
  $nebula_name $geometry_path $file $material_path
  num+=1
done

