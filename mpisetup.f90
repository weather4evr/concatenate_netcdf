module mpisetup

! Adapted from .../GSI/src/enkf/mpisetup.f90

use kinds, only: r_kind, r_single, r_double, i_kind

implicit none

! mpi definitions.
include 'mpif.h'

private

public :: mpi_initialize, mpi_cleanup, mpi_datacounts, stop2
public :: make_mpi_subcommunicator
public :: displs, scount, displs1d, scount1d

integer(i_kind), public :: npe, mype
integer(i_kind), public :: mype_out = 0
integer(i_kind), public :: mpi_status(mpi_status_size)
integer(i_kind), public :: mpi_realkind
integer(i_kind), public :: new_comm, npe_new, mype_new
character(len=4), public :: pe_name

integer(i_kind),allocatable,dimension(:) :: displs, scount
integer(i_kind),allocatable,dimension(:) :: displs1d, scount1d

contains

subroutine mpi_initialize()
integer(i_kind) :: ierr
call mpi_init(ierr)
! mype is process number, npe is total number of processes.
call mpi_comm_rank(mpi_comm_world,mype,ierr)
call mpi_comm_size(mpi_comm_world,npe,ierr)
write(pe_name,'(i4.4)') mype

if (mype == mype_out) print *,'running on ',npe,' processors ...'
if (r_kind == r_single) then
   mpi_realkind = mpi_real4
else if (r_kind == r_double) then
   mpi_realkind = mpi_real8
else
   print *,'illegal r_kind (must be single or double)'
   call mpi_cleanup()
endif

end subroutine mpi_initialize

subroutine mpi_cleanup()
integer(i_kind) :: ierr
flush(6,err=10)
flush(0,err=10)
10 continue
call mpi_barrier(mpi_comm_world,ierr)
if (mype == mype_out ) write(6,*) 'all done!'
call mpi_finalize(ierr)
if (ierr /= 0) then
 print *, 'MPI_Finalize error status = ',ierr
end if
if( allocated(displs)) deallocate(displs)
if( allocated(scount)) deallocate(scount)
if( allocated(displs1d)) deallocate(displs1d)
if( allocated(scount1d)) deallocate(scount1d)
end subroutine mpi_cleanup

subroutine mpi_datacounts(dims,nobs_curr)
   integer(i_kind), intent(in) :: dims(4)
   integer(i_kind), intent(in) :: nobs_curr
   integer(i_kind) :: ierr, i

   call mpi_allgather(product(dims), 1, mpi_integer, scount(0:npe_new-1), 1, mpi_integer, new_comm ,ierr)
   displs(0) = 0
   do i = 1, npe_new-1
      displs(i) = scount(i-1) + displs(i-1)
   enddo

!  if ( mype == mype_out ) write(*,*) ' scount = ',scount
   ! For true 2d variables (not 1-d character vars)
   !  process 1 level at a time, so we need to know the
   !  correct counts for nobs_curr, and not product(dims)
   call mpi_allgather(nobs_curr, 1, mpi_integer, scount1d(0:npe_new-1), 1, mpi_integer, new_comm,ierr)
   displs1d(0) = 0
   do i = 1, npe_new-1
      displs1d(i) = scount1d(i-1) + displs1d(i-1)
   enddo
!  if ( mype == mype_out ) write(*,*) ' scount1d = ',scount1d
end subroutine mpi_datacounts

subroutine make_mpi_subcommunicator(nobs)
   integer(i_kind), intent(in) :: nobs
   integer(i_kind) :: color, ierr

   if ( nobs <= 0 ) then  ! these processers don't have obs, so. don't do anything with them
      color = -1 !MPI_UNDEFINED ! don't add to the new communicator (new_comm)
   else
      color = 0
   endif
   ! communicator new_comm will refer to the processors that actually do stuff
   ! important to note that mype_out will be part of new_comm
   call mpi_comm_split(mpi_comm_world, color, mype, new_comm, ierr) ! use new_comm in later mpi_gather
   call mpi_comm_rank(new_comm,mype_new,ierr)
   call mpi_comm_size(new_comm,npe_new,ierr)
  !write(*,*)'mype, new_comm rank/size,color ',mype,mype_new,npe_new,color
   if ( mype == mype_out ) write(*,*)'number of procs with nobs > 0 is ',npe_new

   ! This gets called for each variable, so need to deallocate each time, before reallocating
   if( allocated(displs)) deallocate(displs)
   if( allocated(scount)) deallocate(scount)
   if( allocated(displs1d)) deallocate(displs1d)
   if( allocated(scount1d)) deallocate(scount1d)
   allocate(displs(0:npe_new-1), scount(0:npe_new-1))
   allocate(displs1d(0:npe_new-1), scount1d(0:npe_new-1))

  !if ( mype == mype_out ) then
  !   if ( nobs <= 0 ) mype_out = mype_out + 1
  !endif
  !if ( mype == mype_out ) write(*,*)'outputting files from pe = ',mype_out

end subroutine make_mpi_subcommunicator

subroutine stop2(ierror_code)
! adapted from GSI/src/main/stop1.f90

  integer(i_kind), intent(in) :: ierror_code
  integer(i_kind)             :: ierr

  write(6,*)'****STOP2****  ABORTING EXECUTION w/code=',ierror_code
  write(0,*)'****STOP2****  ABORTING EXECUTION w/code=',ierror_code
  call mpi_abort(mpi_comm_world,ierror_code,ierr)
  stop
  return
end subroutine stop2

end module mpisetup
