#!/bin/sh

# Shell script which reads and plots the data found in a folders output file sequentially
electron_path_folder="/home/richarddejong/nebula_test_files/vox_tri_pri/multi_runs/10_05_21/*"

declare -i num=1
read_name="/home/richarddejong/nebula_vsc/build/bin/nebula_cpu_vsc"

clear
cd ~
#nebula_path geometry_path electron_path material_path
for file in $electron_path_folder
do
  echo "Processing $file file... $num"
  python $read_file $file
  num+=1
done