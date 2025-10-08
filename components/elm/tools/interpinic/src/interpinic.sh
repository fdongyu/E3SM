module load intel/19.0.4
#module load netcdf/4.6.3
module load netcdf/4.7.4

#export NETCDF_DIR=/share/apps/netcdf/4.6.3/intel/19.0.4
export NETCDF_DIR=/share/apps/netcdf/4.7.4/intel/19.0.4
export LIB_NETCDF=${NETCDF_DIR}/lib
export INC_NETCDF=${NETCDF_DIR}/include
export USER_FC=ifort
export USER_CC=icc
#export USER_LDFLAGS="-L${LIB_NETCDF} -lnetcdf -lnetcdff -lnetcdf_intel"
export USER_LDFLAGS="-L${LIB_NETCDF} -lnetcdf -lnetcdff"
make
