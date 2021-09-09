#!/bin/sh

#Shell script to run simulations using desktop with access to nebula_vsc and nebula_test_files repositories.

cd ..
cd nebula_test_files

#Needed Input Paths
geometry_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/vox_tri_pri/401_401_1001_sd_34000_sh_996_detect_dome_mirror_vs_300.bin"
electron_path_folder="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/vox_tri_pri/multi_runs/21_8_2021_helios_c1/*"
#material_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/tungsten.mat.hdf5"
#material_path_2="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/tungsten_adsorbate.mat.hdf5"
material_path="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/graphite-phonon.mat.hdf5"
material_path_2="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/graphite-phonon_adsorbate.mat.hdf5"

declare -i num=1

#nebula_path geometry_path electron_path material_path material_path_2
for file in $electron_path_folder
do
  echo "Processing $file file... $num"
  nebula_cpu_vsc $geometry_path $file $material_path $material_path_2
  num+=1
done

