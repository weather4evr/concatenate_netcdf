#!/bin/csh
rm *.o *.mod *.x
mpif90 -L${NETCDF}/lib -lnetcdf -lnetcdff -I${NETCDF}/include -c kinds.f90
mpif90 -L${NETCDF}/lib -lnetcdf -lnetcdff -I${NETCDF}/include -c mpisetup.f90
mpif90 -L${NETCDF}/lib -lnetcdf -lnetcdff -I${NETCDF}/include -c netcdf_mod.f90

mpif90 -L${NETCDF}/lib -lnetcdf -lnetcdff -I${NETCDF}/include kinds.o mpisetup.o netcdf_mod.o concatenate_netcdf_files.f90 -o concatenate_netcdf_files.x

exit
