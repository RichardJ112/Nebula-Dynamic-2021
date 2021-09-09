#!/bin/sh

# Shell script which reads and plots the data found in a folders output file sequentially
electron_path_output_folder="/mnt/c/Users/richa/Documents/repos/nebula_test_files/output/May/20_5_2021/*"
electron_path_output_file="/mnt/c/Users/richa/Documents/repos/nebula_test_files/output/May/20_5_2021/10keV_1_50kpp_pitch_0_401_401_801_sb_1000_sd_4000_tungsten_SM_output.bin"


declare -i num=1
read_script="/mnt/c/Users/richa/Documents/repos/nebula_test_files/read_geometry_bin.py"

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
