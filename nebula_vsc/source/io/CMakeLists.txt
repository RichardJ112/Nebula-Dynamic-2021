list(APPEND SOURCE_FILES hdf5_file.cpp)
list(APPEND SOURCE_FILES load_pri_file.cpp)
list(APPEND SOURCE_FILES load_vox_file.cpp)
add_library(io OBJECT ${SOURCE_FILES})
