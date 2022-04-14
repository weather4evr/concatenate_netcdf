program concatenate_netcdf_files

use netcdf
use netcdf_mod, only: open_netcdf, close_netcdf, get_netcdf_dims, define_output_file_from_template, &
                      get_and_output_netcdf_var, get_netcdf_info, get_and_output_netcdf_var_2d, &
                      setup_groups, large_variable_support, numgrps, ncgroup_ids_out, ncgroup_ids_in
use kinds, only : i_kind, r_kind
use mpisetup, only : mpi_initialize, mpi_cleanup, mpi_datacounts, mype, npe, stop2, mype_out, pe_name, &
                       make_mpi_subcommunicator, new_comm

implicit none

integer(i_kind), parameter :: nmax = 300
character(len=500) obsfile, varname, obspath, output_path, outname, fname_prefx, group_name
character(len=100) :: obs_platforms(nmax) = ''
character(len=100) :: ob_platform, suffix
character(len=NF90_MAX_NAME) :: dim_names(4)
integer(i_kind) :: i, itype, ipe, j, nobs_curr, v, nobs_tot, iret, g
integer(i_kind) :: ndims, nvars, num_obs_platforms, nvars_group
integer(i_kind) :: ncfileid, ncidout, ncvarid, rcode
integer(i_kind) :: xtype, ndims_var, natts, dimids(10), dims(4), char_len
integer(i_kind) :: iunit = 10
integer(i_kind) :: varids(1000)
logical :: fexist, just_copy
logical, allocatable, dimension(:) :: fexist_all
integer, allocatable, dimension(:) :: nobs_all

namelist /share/ obspath, output_path, obs_platforms, fname_prefx, large_variable_support

! mpi definitions.
include 'mpif.h'

!!!!!

! Initialize mpi (sets mype, npe, pe_name)
call mpi_initialize

! All processors read the namelist
open(file='input.nml',unit=iunit, status = 'old')
read(iunit,nml=share)
close(iunit)

! All processors figure out the number of ob platforms
num_obs_platforms = 0
do i = 1,nmax
   if ( obs_platforms(i).eq.'' ) cycle
   num_obs_platforms = num_obs_platforms + 1
enddo

allocate(fexist_all(npe))
allocate(nobs_all(0:npe-1))

do itype=1, num_obs_platforms ! loop over the obs platforms you specified

   ob_platform = obs_platforms(itype) ! 'e.g., sondes'

   if ( mype == mype_out ) then 
      write(*,*)' '
      write(*,*)'-----------------------------------'
      write(*,*)'Processing ',trim(adjustl(ob_platform))
      write(*,*)'-----------------------------------'
      write(*,*)' '
   endif
   
   ! each processor reads its own file
  !obsfile = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(ob_platform))//'_'//pe_name//'.nc4'
   obsfile = trim(adjustl(obspath))//'/'//trim(adjustl(fname_prefx))//'_'//trim(adjustl(ob_platform))//'_'//pe_name//'.nc4'
   suffix = 'nc4'
   inquire(file=obsfile,exist=fexist)
   if (.not.fexist) then
      obsfile = trim(adjustl(obspath))//'/'//trim(adjustl(fname_prefx))//'_'//trim(adjustl(ob_platform))//'_'//pe_name//'.h5'
      suffix = 'h5'
   endif

   ! Make sure file exists; all processors skip if file is missing on any processor
   inquire(file=obsfile,exist=fexist)
   call mpi_allgather(fexist, 1, mpi_logical, fexist_all(:), 1, mpi_logical, mpi_comm_world,iret)
   if (.not.all(fexist_all)) then
      write(*,*)'missing file: ',trim(adjustl(obsfile))
      call mpi_barrier(mpi_comm_world,iret)  ! sync up before going to next ob type
      cycle
   endif

   ! Define output file name; output from mype_out
  !outname = trim(adjustl(output_path))//'/'//trim(adjustl(fname_prefx))//'_'//trim(adjustl(ob_platform))//'_all.nc4'
   outname = trim(adjustl(output_path))//'/'//trim(adjustl(fname_prefx))//'_'//trim(adjustl(ob_platform))//'_all.'//suffix

   ! File exists.  Open it, figure out how many obs this file has.
   call open_netcdf(obsfile,ncfileid)
   call get_netcdf_dims(ncfileid,'nlocs',nobs_curr)
  !write(*,fmt='(a7,a4,a12,i8,1x,a3)')' on pe ',pe_name,' there are ',nobs_curr, 'obs'

   ! Get number of dimensions and number of variables in file
   call get_netcdf_info(obsfile,ncfileid,ndims,nvars) ! ndims not used

   ! Figure out if this file has "groups". If so, we need to force large_variable_support = true (done in subroutine).
   !  This will output a netCDF4 file, and "groups" are only allowed in netCDF4 files.
   !  Also, if there are indeed "groups", we can assume we are reading netCDF4 files.
   !  See here: https://www.unidata.ucar.edu/software/netcdf/docs-fortran/f90_groups.html
   ! numgrps is the number of groups, and the group ids in input files are in "ncgroup_ids_in"...we can refer to the group by its id
   call setup_groups(ncfileid) ! sets numgrps, ncgroup_ids_out, ncgroup_ids_in (from netcdf_mod)

   ! Now figure out the total number of observations across all files/processors
   call mpi_allreduce(nobs_curr, nobs_tot, 1, mpi_integer, mpi_sum, mpi_comm_world, iret)

   ! Need to make sure mype_out has nobs_curr > 0....
    call mpi_allgather(nobs_curr, 1, mpi_integer, nobs_all(:), 1, mpi_integer, mpi_comm_world, iret)
    do i = 0,npe-1
       if ( nobs_all(i) > 0 ) then
          mype_out = i
          exit
       endif
    enddo

   ! Make subcommunicator...we only want to do operations where nobs_curr > 0
   ! Basically, if nobs_curr <= 0, we don't want it to participate in the concatenating
   !   Easiest way to do that is to just put those processors on a different mpi communicator
   !   ! The new communicator is new_comm, from mpisetup_mod
   call make_mpi_subcommunicator(nobs_curr)

   if ( nobs_curr <= 0 ) numgrps = -9999 ! This will bypass the loop over variables, and basically skip everything

   ! Create the output netcdf file that will have all concatenated obs over all processors.
   ! Define all dimensions and variables, but don't yet fill output file with data
   ! Also fills ncgroup_ids_out (from netcdf_mod)
   if ( mype == mype_out ) then
      call define_output_file_from_template(ncfileid,outname,nobs_tot,ncidout) ! ncidout only on mype_out...seems to be ok.
      write(*,fmt='(a20,i8)')' total num obs = ',nobs_tot
   endif

   numgrps = numgrps + 1 ! Need to also loop over the "root" group, which is basically just the main file.  
                         !   For g = 1, we are effectively outputting traditional variables not in a group
                         !   It's okay to re-use this variable at this point in the code.
                         ! Note that ncgroup_ids_in(1) = ncfileid and ncgroup_ids_out(1) = ncidout

   ! Each processor loops over all groups and variables in the file
   do g = 1,numgrps

   rcode = nf90_inq_varids(ncgroup_ids_in(g), nvars_group, varids) ! get all the variables in the ith group from input file
   rcode = nf90_inq_grpname(ncgroup_ids_in(g), group_name)
   if ( mype == mype_out ) write(*,fmt='(a,i3)')'number of vars in '//trim(adjustl(group_name))//' group is ',nvars_group

   do v = 1,nvars_group

      dims(:) = 1 ! reset each time through loop
      dim_names(:) = ''
      rcode = nf90_Inquire_Variable(ncgroup_ids_in(g), v, varname, xtype, ndims_var, dimids, natts) ! Output is varname, xtype,ndims_var, dimids
      do j = 1,ndims_var
         rcode = nf90_inquire_dimension( ncgroup_ids_in(g), dimids(j), name=dim_names(j), len=dims(j) )
         dim_names(j) = trim(adjustl(dim_names(j))) ! can't do trim on an array, so do it here
      enddo
      if ( mype == mype_out ) write(*,fmt='(a24, a30,i2,4i8)')'variable, ndims, dims = ',trim(adjustl(varname)),ndims_var,dims(1:ndims_var)

      ! collect information about amount of data (total number of points) on each processor
      !  we will send product(dims) from each processor to account for 
      !  possibility of multi-dimensional variables
      !  after we have done that (and filled $scount), we can compute the displacement ($displs)
      !  needed for mpi_gatherv
      call mpi_datacounts(dims,nobs_curr) ! fills $scount, $displs in mpisetup

      ! Assume the variable needs to be concatenated across all processors (just_copy = .false.)
      ! But, if variable does not have size nobs_curr, just make a copy of it; nothing to concatenate (just_copy = .true.)
      !   A better check is to look at the dimension name--if it's "nlocs" then we need to concaenate, otherwise, just_copy = .true.
      just_copy = .false. ! assume the variable needs to be concatenated across all processors
      if ( ndims_var == 1 ) then
!        if ( dims(1) /= nobs_curr ) just_copy = .true.
         if ( dim_names(1) /= 'nlocs' ) just_copy = .true.
         call get_and_output_netcdf_var(ncgroup_ids_in(g),varname,dims,ncgroup_ids_out(g),nobs_tot,xtype,just_copy)
      else if ( ndims_var == 2 ) then
         ! if a character variable, it's really a 1-d array, with 2nd dimension equal to the number of characters
         if ( xtype == nf90_char ) then
!           if ( dims(2) /= nobs_curr ) just_copy = .true.
            if ( dim_names(2) /= 'nlocs' ) just_copy = .true.
            char_len = dims(1)  ! array of variable is (/ char_len, nobs_curr /)
            call get_and_output_netcdf_var(ncgroup_ids_in(g),varname,dims,ncgroup_ids_out(g),nobs_tot,xtype,just_copy,char_len)
         else if ( xtype == nf90_float .or. xtype == nf90_int .or. xtype == nf90_int64 ) then
           !if ( dims(2) /= nobs_curr ) just_copy = .true. ! array of variable is (/ dim1, nobs_curr /)
            if ( dim_names(2) /= 'nlocs' ) just_copy = .true.
            call get_and_output_netcdf_var_2d(ncgroup_ids_in(g),varname,dims,ncgroup_ids_out(g),nobs_tot,xtype,just_copy)
         endif
      endif
      if ( just_copy .and. mype == mype_out ) write(*,*)'copying variable = ',trim(adjustl(varname))

   enddo ! end loop over variables
   enddo ! end loop over groups

   call close_netcdf(obsfile,ncfileid)
   if ( mype == mype_out ) call close_netcdf(outname,ncidout) ! only open on mype_out

   call mpi_barrier(mpi_comm_world,iret)  ! sync up before going to next ob type
   call mpi_comm_free(new_comm,iret) ! free/clear the sub-communicator

enddo ! end loop over "itype"

! clean up
call mpi_cleanup
deallocate(fexist_all,nobs_all)
stop

end program concatenate_netcdf_files
