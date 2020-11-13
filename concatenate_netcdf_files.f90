program concatenate_netcdf_files

use netcdf
use netcdf_mod, only: open_netcdf, close_netcdf, get_netcdf_dims, define_output_file_from_template, &
                      get_and_output_netcdf_var, get_netcdf_info
use kinds, only : i_kind, r_kind
use mpisetup, only : mpi_initialize, mpi_cleanup, mpi_datacounts, mype, npe, stop2, mype_out, pe_name

implicit none

integer(i_kind), parameter :: nmax = 300
character(len=500) obsfile, varname, obspath, output_path, outname
character(len=100) :: obs_platforms(nmax) = ''
character(len=100) :: ob_platform
integer(i_kind) :: i, itype, ipe, j, nobs_curr, v, nobs_tot, iret
integer(i_kind) :: ndims, nvars, num_obs_platforms
integer(i_kind) :: ncfileid, ncidout, ncvarid, rcode
integer(i_kind) :: xtype, ndims_var, natts, dimids(10), dims(4), char_len
integer(i_kind) :: iunit = 10
logical :: fexist, just_copy
logical, allocatable, dimension(:) :: fexist_all

namelist /share/ obspath, output_path, obs_platforms

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

do itype=1, num_obs_platforms ! loop over the obs platforms you specified

   ob_platform = obs_platforms(itype) ! 'sondes'

   ! output file name; output from mype_out
   outname = trim(adjustl(output_path))//'/obsout_omb_'//trim(adjustl(ob_platform))//'_all.nc4'

   if ( mype == mype_out ) then 
      write(*,*)' '
      write(*,*)'Processing ',trim(adjustl(ob_platform))
      write(*,*)' -----------------------------------'
   endif
   
   ! each processor reads its own file
   obsfile = trim(adjustl(obspath))//'/obsout_omb_'//trim(adjustl(ob_platform))//'_'//pe_name//'.nc4'

   ! Make sure file exists; all processors skip if file is missing on any processor
   inquire(file=obsfile,exist=fexist)
   call mpi_allgather(fexist, 1, mpi_logical, fexist_all(:), 1, mpi_logical, mpi_comm_world,iret)
   if (.not.all(fexist_all)) then
      write(*,*)'missing file: ',trim(adjustl(obsfile))
      call mpi_barrier(mpi_comm_world,iret)  ! sync up before going to next ob type
      cycle
   endif

   ! File exists.  Open it, figure out how many obs this file has.
   call open_netcdf(obsfile,ncfileid)
   call get_netcdf_dims(ncfileid,'nlocs',nobs_curr)
   write(*,fmt='(a7,a4,a12,i8,1x,a3)')' on pe ',pe_name,' there are ',nobs_curr, 'obs'

   ! Get number of dimensions and number of variables in file
   call get_netcdf_info(obsfile,ncfileid,ndims,nvars)

   ! Now figure out the total number of observations across all files/processors
   call mpi_allreduce(nobs_curr, nobs_tot, 1, mpi_integer, mpi_sum, mpi_comm_world, iret)

   ! Create the output netcdf file that will have all concatenated obs over all processors.
   ! Define all dimensions and variables, but don't yet fill with data
   if ( mype == mype_out ) then
      call define_output_file_from_template(ncfileid,outname,nobs_tot,ncidout)
      write(*,fmt='(a20,i8)')' total num obs = ',nobs_tot
   endif

   ! Each processor loops over all variables in the file
   do v = 1,nvars

      dims(:) = 1 ! reset each time through loop
      rcode = nf90_Inquire_Variable(ncfileid, v, varname, xtype, ndims_var, dimids, natts) ! Output is varname, xtype,ndims_var, dimids
      do j = 1,ndims_var
         rcode = nf90_inquire_dimension( ncfileid, dimids(j), len=dims(j) )
      enddo
!     if ( mype == mype_out ) write(*,*)'variable, ndims, dims = ',trim(adjustl(varname)),ndims_var,dims(1:ndims_var)

      ! collect information about amount of data (total number of points) on each processors
      !  we will send product(dims) from each processor to account for 
      !  possibility of multi-dimensional variables
      !  after we have done that (and filled $scount), we can compute the displacement ($displs)
      !  needed for mpi_gatherv
      call mpi_datacounts(dims) ! fills $scount, $displs in mpisetup

      ! Assume the variable needs to be concatenated across all processors (just_copy = .false.)
      ! But, if variable does not have size nobs_curr, just make a copy of it; nothing to concatenate (just_copy = .true.)
      just_copy = .false. ! assume the variable needs to be concatenated across all processors
      if ( ndims_var == 1 ) then
         if ( dims(1) /= nobs_curr ) just_copy = .true.
         call get_and_output_netcdf_var(ncfileid,varname,dims,ncidout,nobs_tot,xtype,just_copy)
      else if ( ndims_var == 2 ) then
         ! if a character variable, it's really a 1-d array, with 2nd dimension equal to the number of characters
         if ( xtype == nf90_char ) then
            if ( dims(2) /= nobs_curr ) just_copy = .true.
            char_len = dims(1)  ! array of variable is (/ char_len, nobs_curr /)
            call get_and_output_netcdf_var(ncfileid,varname,dims,ncidout,nobs_tot,xtype,just_copy,char_len)
         else
           !call get_and_output_netcdf_var(ncfileid,varname,dims,output2d)
            ! mpi_gather all obs on mype_out
          ! call mpi_gatherv(output2d, scount(mype), mpi_char, output2d_allData, scount(:), displs(:), mpi_char, mype_out ,mpi_comm_world,iret)
          ! ! output the variable from mype_out
          !  if ( mype == mype_out ) then
          !     rcode = nf90_inq_varid(ncidout,trim(adjustl(varname)),ncvarid)
          !     rcode = nf90_put_var(ncidout,ncvarid,output2d_allData)
          !  endif
         endif
      endif
      if ( just_copy .and. mype == mype_out ) write(*,*)'copying variable = ',trim(adjustl(varname))

   enddo ! end loop over variables

   call close_netcdf(obsfile,ncfileid)
   if ( mype == mype_out ) call close_netcdf(outname,ncidout) ! only open on mype_out

   call mpi_barrier(mpi_comm_world,iret)  ! sync up before going to next ob type

enddo ! end loop over "itype"

! clean up
call mpi_cleanup
deallocate(fexist_all)
stop

end program concatenate_netcdf_files
