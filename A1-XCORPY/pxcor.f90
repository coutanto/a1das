!
!   pxcor.f90
!   
!
!   Created by coutanto on 11/08/2021.
!
!   This module describe the PXC class that
!   manages the pre-processing to apply on signals to compute
!   cross correlation
!
module pxcor_mod
use iso_fortran_env
implicit none

!
!--------- DEFINITION AND DESCRIPTION OF THE VARIOUS CORRELATION PROCESSING -----------
!
type description
  character (len=10) :: label
  integer :: narg
end type

type(description), parameter :: process_list(*) = &
   (/ &
     description("onebit",   0), &
     description("clip",     1), &
     description("white",    2), &
     description("bandpass", 2), &
     description("lowpass",  1) &
   /)
!--------- END DEFINITION AND DESCRIPTION OF THE VARIOUS CORRELATION PROCESSING -----------

! Definition of the various processing applied in the session
type process
  character(10) :: label
  integer :: narg
  real(kind=real64),allocatable :: args(:)
  logical :: flag
end type process


type pxcor
  type(process), allocatable :: proc(:)
contains
   procedure :: init
   procedure :: flag_raised
end type pxcor

type(pxcor) :: pxc

contains
! ===================
! init
! ===================
subroutine init(this, d)
   class(pxcor) :: this
   type(description) :: d(:)

   integer :: nproc, i

   nproc=size(d)
   allocate(this%proc(nproc))
   do i=1,nproc
      this%proc(i)%label = trim(d(i)%label)
      this%proc(i)%narg  = d(i)%narg
      this%proc(i)%flag  = .false.
   enddo
end subroutine

! ===================
! process_info
! ===================
subroutine pxc_print_info()
   integer :: i

   if (.not. allocated(pxc%proc))  call pxc%init(process_list)
   write(0,*) 'processing available for cross-correlation:'
   do i=1,size(pxc%proc)
      write(0,*) ' <',trim(pxc%proc(i)%label),'> taking ',pxc%proc(i)%narg,' parameters'
   enddo
end subroutine

! ===================
! register_proc
   !
   ! register a process to be performed during correlation
   ! input:
   !     proc = name of the process
   !     val = array of parameters
   !     nval = number of parameters
   ! output:
   !     ier = 0 success, 1 wrong process name, 2 wrong argument number
! ===================
!TODO message d'erreur si pas bon nombre de parametre ou si mauvais
! processing
subroutine pxc_register_proc(proc, val, nval, ier)
   character (len=10) :: proc
   integer      :: nval, ier
   real(kind=8) :: val(nval)

   integer :: i

   ier = 0
   ! first call to register_proc, initialize
   if (.not. allocated(pxc%proc))  call pxc%init(process_list)
   ! look for process in the list and check parameter number
   do i=1,size(pxc%proc)
      if (trim(pxc%proc(i)%label) /= trim(proc)) cycle
      if (nval /= pxc%proc(i)%narg) then
         ier = 2
         return
      end if
      if (allocated(pxc%proc(i)%args)) deallocate(pxc%proc(i)%args)
      allocate(pxc%proc(i)%args(nval))
      pxc%proc(i)%args = val
      pxc%proc(i)%flag = .true.
      exit
   enddo
   ! test if we didn't find the process
   if (i>size(pxc%proc)) ier = 1

end subroutine
! ===================
! flag_raised
! ===================
logical function flag_raised(this, proc) result(flag)
   class(pxcor) :: this
   character(len=*) :: proc

   integer :: i

   if (.not. allocated(pxc%proc))  call pxc%init(process_list)
   flag = .false.
   do i=1,size(this%proc)
      if (trim(this%proc(i)%label) /= trim(proc)) cycle
      flag = this%proc(i)%flag
      exit
   enddo
   if (i==size(this%proc)+1) write(0,*) 'Warning: processing ',proc,' unknown'

end function


end module pxcor_mod
