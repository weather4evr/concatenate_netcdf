!ifort  -L${NETCDF}/lib -lnetcdf -lnetcdff -I${NETCDF}/include test.f90 -o a.out
module netcdf_mod

use netcdf
use kinds, only : i_kind
use mpisetup, only: stop2, npe, mype, mype_out, displs, scount, displs1d, scount1d, new_comm

implicit none

!include 'netcdf.inc'

private

! public subroutines
public :: open_netcdf, close_netcdf, get_netcdf_dims, define_output_file_from_template, get_and_output_netcdf_var
public :: get_netcdf_info, get_and_output_netcdf_var_2d, setup_groups

! parameter variable visible to this module only
integer(i_kind), parameter :: nmax_groups = 500

! public variables
logical, public :: large_variable_support = .false. 
integer(i_kind), public :: numgrps, ncgroup_ids_in(nmax_groups), ncgroup_ids_out(nmax_groups)

! variables visible to this module only
character(len=100) :: DIMSNAME,VARSNAME,ATTSNAME
integer(i_kind) :: ncidin,ncidout,ndims,nvars,ngatts,unlimdimid,idims,dimsval,ivars,idims2,ivars2,&
           varstype,varsndims,varsdimids(4),varsnatts, ivarsnatts,igatts
integer(i_kind) :: ncvarid, ierr, iret, rcode, ncstatus, i
integer(i_kind) :: ncgroup_ids(nmax_groups)

! mpi definitions.
include 'mpif.h'

!interface get_and_output_netcdf_var_2d
!  module procedure get_and_output_netcdf_var_2d_real
!   module procedure get_and_output_netcdf_var_2d_integer
!end interface

!interface get_and_output_netcdf_var
!   module procedure get_netcdf_var_1d_real
!   module procedure get_netcdf_var_1d_integer
!   module procedure get_netcdf_var_1d_char
!end interface

contains

subroutine get_and_output_netcdf_var_2d(fileid,variable,dims,ncidout,nobs_tot,xtype,just_copy)
   integer(i_kind),  intent(in) :: fileid, ncidout,nobs_tot, xtype
   character(len=*), intent(in) :: variable
   integer(i_kind),  intent(in), dimension(4) :: dims
   logical,          intent(in)           :: just_copy
   if ( xtype == nf90_float ) then
      call get_and_output_netcdf_var_2d_real(fileid,variable,dims,ncidout,nobs_tot,just_copy)
   elseif ( xtype == nf90_int .or. xtype == nf90_int64 ) then
      call get_and_output_netcdf_var_2d_integer(fileid,variable,dims,ncidout,nobs_tot,just_copy)
   else
      write(*,*)'Unsure what to do with '//trim(adjustl(variable))
      write(*,*)'What type of variable is it?'
      call stop2(53)
   endif
end subroutine get_and_output_netcdf_var_2d

subroutine get_and_output_netcdf_var(fileid,variable,dims,ncidout,nobs_tot,xtype,just_copy,char_len)
   integer(i_kind),  intent(in) :: fileid, ncidout,nobs_tot, xtype
   character(len=*), intent(in) :: variable
   integer(i_kind),  intent(in), dimension(4) :: dims
   logical,          intent(in)           :: just_copy
   integer(i_kind),  intent(in), optional :: char_len
   if ( xtype == nf90_float ) then
      call get_netcdf_var_1d_real(fileid,variable,dims,ncidout,nobs_tot,just_copy)
   elseif ( xtype == nf90_int .or. xtype == nf90_int64 ) then
      call get_netcdf_var_1d_integer(fileid,variable,dims,ncidout,nobs_tot,just_copy)
   elseif ( xtype == nf90_char) then
      call get_netcdf_var_1d_char(fileid,variable,dims,ncidout,nobs_tot,just_copy,char_len)
   elseif ( xtype == nf90_string) then
      if ( mype == mype_out ) write(*,*)'Variable '//trim(adjustl(variable))//' is a string; skipping output'
   else
      write(*,*)'Unsure what to do with '//trim(adjustl(variable))
      write(*,*)'What type of variable is it?'
      call stop2(54)
   endif
end subroutine get_and_output_netcdf_var

subroutine open_netcdf(fname,ncfileid)
   character(len=*), intent(in) :: fname
   integer(i_kind), intent(out) :: ncfileid

   ncstatus = nf90_open(path=trim(adjustl(fname)),mode=nf90_nowrite,ncid=ncfileid)  ! open file
   if ( ncstatus .eq. 0 ) then
!     write(*,fmt='(a)') 'opened '//trim(adjustl(fname))//' for reading'
!     write(*,fmt='(a,i8)') 'fileid = ',ncfileid
   else
      write(*,fmt='(a)') 'error reading '//trim(adjustl(fname))
      call stop2(31) ! stop
   endif
end subroutine open_netcdf

subroutine close_netcdf(fname,ncfileid)
   character(len=*), intent(in) :: fname
   integer(i_kind), intent(in) :: ncfileid
   ncstatus = nf90_close(ncfileid) ! close file
   if ( ncstatus .ne. 0 ) then
      write(*,fmt='(a)') 'error closing '//trim(adjustl(fname))
      call stop2(32) ! stop
   endif
end subroutine close_netcdf

subroutine get_netcdf_dims(ncfileid,variable,output)
   integer(i_kind), intent(in) :: ncfileid
   character(len=*), intent(in) :: variable
   integer(i_kind), intent(out) :: output

   integer(i_kind) :: ncdimid, ierr

   ierr = 0
   ncstatus = nf90_inq_dimid(ncfileid,trim(adjustl(variable)),ncdimid) ; ierr = ierr + ncstatus
   ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,len=output)      ; ierr = ierr + ncstatus
   if ( ierr /= 0 ) then
      write(0,*) 'Error reading dimension for '//trim(adjustl(variable))
      write(0,*) 'ierr = ',ierr
      call stop2(33) ! stop
   endif

end subroutine get_netcdf_dims

subroutine get_netcdf_info(fname,ncfileid,ndims,nvars)
   character(len=*), intent(in) :: fname
   integer(i_kind), intent(in) :: ncfileid
   integer(i_kind), intent(inout) :: ndims,nvars
   ncstatus=nf90_inquire(ncfileid,ndims,nvars,ngatts,unlimdimid)
   if ( ncstatus /= 0 ) then
      write(0,*) 'Error getting data for for '//trim(adjustl(fname))
      write(0,*) 'ncstatus = ',ncstatus
      call stop2(34) ! stop
   endif
end subroutine get_netcdf_info

subroutine setup_groups(ncfileid)
   integer(i_kind), intent(in)  :: ncfileid
   ncstatus = nf90_inq_grps(ncfileid, numgrps, ncgroup_ids)
   ncgroup_ids_in(1) = ncfileid
   ncgroup_ids_in(2:numgrps) = ncgroup_ids
   ! If there are groups, we can assume we read a netCDF4 file, and so we need
   !   to output a netCDF4 file, too.  This can be done by simply setting
   !   large_variable_support = true.
   if ( numgrps > 0 ) then
      if ( mype == mype_out ) write(*,*)'There are ',numgrps,' groups in the file'
      large_variable_support = .true.
   endif
end subroutine setup_groups

subroutine define_output_file_from_template(ncidin,fout,nobs_tot,ncidout)
   integer(i_kind), intent(in)    :: ncidin
   character(len=*), intent(in) :: fout
   integer(i_kind), intent(in)    :: nobs_tot
   integer(i_kind), intent(inout) :: ncidout

   integer(i_kind) :: netcdf_file_type, group_ncid
   integer(i_kind) :: varids(1000)
   character(len=500) :: group_name

   ! this variable is read in from namelist
   if ( large_variable_support ) then
      netcdf_file_type = NF90_NETCDF4
   else
      netcdf_file_type = NF90_CLOBBER
   endif

!  rcode=nf90_open(path=trim(adjustl(fin)),mode=nf90_nowrite,ncid=ncidin)
!  rcode=nf90_create(path=trim(adjustl(fout)),cmode=nf90_clobber,ncid=ncidout)
!  rcode=nf90_create(path=trim(adjustl(fout)),cmode=NF90_NETCDF4,ncid=ncidout) ! make netcdf4 format...this should be able to overwrite existing files
   rcode=nf90_create(path=trim(adjustl(fout)),cmode=netcdf_file_type,ncid=ncidout)

   ! dimensions
   rcode=nf90_inquire(ncidin,ndims,nvars,ngatts,unlimdimid)
   do idims=1,ndims
      rcode=nf90_inquire_dimension(ncidin,idims,DIMSNAME,dimsval)
      if ( trim(adjustl(DIMSNAME)) == 'nlocs') then
         rcode=nf90_def_dim(ncidout,trim(adjustl(DIMSNAME)),nobs_tot,idims2)
      else
         rcode=nf90_def_dim(ncidout,trim(adjustl(DIMSNAME)),dimsval,idims2)
      endif
   end do

   ! variables
   do ivars=1,nvars
      rcode=nf90_inquire_variable(ncidin,ivars,VARSNAME,varstype,varsndims,varsdimids,varsnatts)
      rcode=nf90_def_var(ncidout,trim(adjustl(VARSNAME)),varstype,varsdimids(1:varsndims),ivars2)
     !rcode=nf90_def_var(ncidout,VARSNAME,varstype,varsndims,varsdimids,ivars)
      do ivarsnatts=1,varsnatts
         rcode=nf90_inq_attname(ncidin,ivars,ivarsnatts,ATTSNAME)
         rcode=nf90_copy_att(ncidin,ivars,ATTSNAME,ncidout,ivars2)
      end do
   end do

   ! groups...numgrps could be 0, in which case this block doesn't get executed, which is ok.
  !rcode = nf90_inq_grps(ncidin, numgrps, ncgroup_ids)
   ncgroup_ids_out(1) = ncidout ! Need to keep the ID for the output file as a whole (also the "root" group)
   do i = 1,numgrps
      rcode = nf90_inq_grpname(ncgroup_ids(i), group_name)
      rcode = nf90_def_grp(ncidout, group_name, group_ncid) ! define group in output file; refer to group by $group_nicd
      rcode = nf90_inq_varids(ncgroup_ids(i), nvars, varids) ! get all the variables in the ith group from input file
      do ivars=1,nvars ! loop over the number of variables in the group
         rcode=nf90_inquire_variable(ncgroup_ids(i),ivars,VARSNAME,varstype,varsndims,varsdimids,varsnatts)
         rcode=nf90_def_var(group_ncid,trim(adjustl(VARSNAME)),varstype,varsdimids(1:varsndims),ivars2)
        !rcode=nf90_def_var(ncidout,VARSNAME,varstype,varsndims,varsdimids,ivars)
         do ivarsnatts=1,varsnatts
            rcode=nf90_inq_attname(ncgroup_ids(i),ivars,ivarsnatts,ATTSNAME)
            rcode=nf90_copy_att(ncgroup_ids(i),ivars,ATTSNAME,group_ncid,ivars2)
         end do
      end do
      ncgroup_ids_out(i+1) = group_ncid ! Add 1 because element 1 is for ncidout
   end do

   ! global attributes
   do igatts=1,ngatts
      rcode=nf90_inq_attname(ncidin,nf90_global,igatts,ATTSNAME)
      rcode=nf90_copy_att(ncidin,nf90_global,ATTSNAME,ncidout,nf90_global)
   end do
   rcode=nf90_enddef(ncidout)
end subroutine define_output_file_from_template

subroutine get_netcdf_var_1d_real(fileid,variable,dims,ncidout,nobs_tot,just_copy)

   integer(i_kind), intent(in) :: fileid, ncidout,nobs_tot
   character(len=*), intent(in) :: variable
   integer(i_kind), intent(in), dimension(4) :: dims
   logical, intent(in) :: just_copy

   real, dimension(dims(1))         :: output
   real, allocatable, dimension(:)  :: output_allData

   ierr = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid) ; ierr = ierr + ncstatus
   ncstatus = nf90_get_var(fileid,ncvarid,output)                    ; ierr = ierr + ncstatus

   if ( ierr /= 0 ) then
      write(0,*) 'Error reading data for '//trim(adjustl(variable))
      write(0,*) 'ierr = ',ierr
      call stop2(35) ! stop
   endif

   if ( just_copy ) then
      allocate(output_allData(size(output)))
      output_allData = output
   else
      allocate(output_allData(nobs_tot))
      call mpi_gatherv(output, scount(mype), mpi_real, output_allData, scount(:), displs(:), mpi_real, mype_out ,new_comm,iret)
   endif

   ! output the variable from mype_out
    if ( mype == mype_out ) then
       rcode = nf90_inq_varid(ncidout,trim(adjustl(variable)),ncvarid)
       rcode = nf90_put_var(ncidout,ncvarid,output_allData)
    endif
   deallocate(output_allData)

end subroutine get_netcdf_var_1d_real

subroutine get_netcdf_var_1d_integer(fileid,variable,dims,ncidout,nobs_tot,just_copy)

   integer(i_kind), intent(in) :: fileid, ncidout,nobs_tot
   character(len=*), intent(in) :: variable
   integer(i_kind), intent(in), dimension(4) :: dims
   logical, intent(in) :: just_copy

   integer(i_kind), dimension(dims(1))         :: output
   integer(i_kind), allocatable, dimension(:)  :: output_allData

   ierr = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid) ; ierr = ierr + ncstatus
   ncstatus = nf90_get_var(fileid,ncvarid,output)                    ; ierr = ierr + ncstatus

   if ( ierr /= 0 ) then
      write(0,*) 'Error reading data for '//trim(adjustl(variable))
      write(0,*) 'ierr = ',ierr
      call stop2(36) ! stop
   endif

   if ( just_copy ) then
      allocate(output_allData(size(output)))
      output_allData = output
   else
      allocate(output_allData(nobs_tot))
      call mpi_gatherv(output, scount(mype), mpi_integer, output_allData, scount(:), displs(:), mpi_integer, mype_out ,new_comm,iret)
   endif

   ! output the variable from mype_out
    if ( mype == mype_out ) then
       rcode = nf90_inq_varid(ncidout,trim(adjustl(variable)),ncvarid)
       rcode = nf90_put_var(ncidout,ncvarid,output_allData)
    endif
   deallocate(output_allData)

end subroutine get_netcdf_var_1d_integer

subroutine get_netcdf_var_1d_char(fileid,variable,dims,ncidout,nobs_tot,just_copy,char_len)

   integer(i_kind), intent(in) :: fileid, ncidout,nobs_tot,char_len
   character(len=*), intent(in) :: variable
   integer(i_kind), intent(in), dimension(4) :: dims
   logical, intent(in) :: just_copy

   character(len=char_len), dimension(dims(2)) :: output ! char_len == dims(1) is number of characters, dims(2) is number of obs
   character(len=char_len), allocatable, dimension(:)        :: output_allData

   ierr = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid) ; ierr = ierr + ncstatus
   ncstatus = nf90_get_var(fileid,ncvarid,output)                    ; ierr = ierr + ncstatus

   if ( ierr /= 0 ) then
      write(0,*) 'Error reading data for '//trim(adjustl(variable))
      write(0,*) 'ierr = ',ierr
      call stop2(37) ! stop
   endif

!  write(*,*)'char_len = ',char_len
!  write(*,*) ' scount = ',scount
!  write(*,*) ' size, sizeof = ',size(output),sizeof(output)
!  write(*,*)output

   if ( just_copy ) then
      allocate(output_allData(size(output)))
      output_allData = output
   else
      allocate(output_allData(nobs_tot))
      call mpi_gatherv(output, scount(mype), mpi_char, output_allData, scount(:), displs(:), mpi_char, mype_out ,new_comm,iret)
   endif
!  write(*,*) ' size1, sizeof1 = ',size(output_allData),sizeof(output_allData)

   ! output the variable from mype_out
    if ( mype == mype_out ) then
       rcode = nf90_inq_varid(ncidout,trim(adjustl(variable)),ncvarid)
       rcode = nf90_put_var(ncidout,ncvarid,output_allData)
    endif
   deallocate(output_allData)

end subroutine get_netcdf_var_1d_char

subroutine get_and_output_netcdf_var_2d_real(fileid,variable,dims,ncidout,nobs_tot,just_copy)

   integer(i_kind), intent(in) :: fileid, ncidout,nobs_tot
   character(len=*), intent(in) :: variable
   integer(i_kind), intent(in), dimension(4) :: dims
   logical, intent(in) :: just_copy

   real, dimension(dims(1),dims(2))         :: output
   real, allocatable, dimension(:,:)  :: output_allData

   ierr = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid) ; ierr = ierr + ncstatus
   ncstatus = nf90_get_var(fileid,ncvarid,output)                    ; ierr = ierr + ncstatus

   if ( ierr /= 0 ) then
      write(0,*) 'Error reading data for '//trim(adjustl(variable))
      write(0,*) 'ierr = ',ierr
      call stop2(38) ! stop
   endif

   if ( just_copy ) then
      allocate( output_allData(dims(1),dims(2)) ) ! nlevs, nobs_curr
      output_allData = output
   else
      allocate(output_allData(dims(1),nobs_tot)) ! nlevs, nobs
     ! Handle one vertical level at a time
      do i = 1,dims(1)
         call mpi_gatherv(output(i,:), scount1d(mype), mpi_real, output_allData(i,:), scount1d(:), displs1d(:), mpi_real, mype_out ,new_comm,iret)
         if ( iret /= 0 ) then
            if ( mype == mype_out ) write(*,*)'mpi_gatherv problem for '//trim(adjustl(variable))
         endif
      enddo
   endif

   ! output the variable from mype_out
    if ( mype == mype_out ) then
       rcode = nf90_inq_varid(ncidout,trim(adjustl(variable)),ncvarid)
       rcode = nf90_put_var(ncidout,ncvarid,output_allData)
    endif
   deallocate(output_allData)

end subroutine get_and_output_netcdf_var_2d_real

subroutine get_and_output_netcdf_var_2d_integer(fileid,variable,dims,ncidout,nobs_tot,just_copy)

   integer(i_kind), intent(in) :: fileid, ncidout,nobs_tot
   character(len=*), intent(in) :: variable
   integer(i_kind), intent(in), dimension(4) :: dims
   logical, intent(in) :: just_copy

   integer, dimension(dims(1),dims(2))         :: output
   integer, allocatable, dimension(:,:)  :: output_allData

   ierr = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid) ; ierr = ierr + ncstatus
   ncstatus = nf90_get_var(fileid,ncvarid,output)                    ; ierr = ierr + ncstatus

   if ( ierr /= 0 ) then
      write(0,*) 'Error reading data for '//trim(adjustl(variable))
      write(0,*) 'ierr = ',ierr
      call stop2(39) ! stop
   endif

   if ( just_copy ) then
      allocate( output_allData(dims(1),dims(2)) ) ! nlevs, nobs_curr
      output_allData = output
   else
      allocate(output_allData(dims(1),nobs_tot)) ! nlevs, nobs
     ! Handle one vertical level at a time
      do i = 1,dims(1)
         call mpi_gatherv(output(i,:), scount1d(mype), mpi_real, output_allData(i,:), scount1d(:), displs1d(:), mpi_real, mype_out ,new_comm,iret)
         if ( iret /= 0 ) then
            if ( mype == mype_out ) write(*,*)'mpi_gatherv problem for '//trim(adjustl(variable))
         endif
      enddo
   endif

   ! output the variable from mype_out
    if ( mype == mype_out ) then
       rcode = nf90_inq_varid(ncidout,trim(adjustl(variable)),ncvarid)
       rcode = nf90_put_var(ncidout,ncvarid,output_allData)
    endif
   deallocate(output_allData)

end subroutine get_and_output_netcdf_var_2d_integer

end module netcdf_mod
