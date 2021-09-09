#!/bin/sh

# Shell script which reads and plots the data found in a folders output file sequentially
electron_path_output_folder="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/vox_tri_pri/multi_runs/9_7_2021/*"

declare -i num=1
read_script="/mnt/c/Users/Richard/source/repos/Nebula/nebula_test_files/read_geometry_bin.py"

#Single file test
#echo "Test single file $electron_path_output_file"
#python3 $read_script $electron_path_output_file

clear
#nebula_path geometry_path electron_path material_path
for file in $electron_path_output_folder
do
  echo "Processing $file file... $num"
  python3 $read_script $file
  num+=1
done
