!
!   a1xcor.f90
!   
!
!   Created by coutanto on 09/08/2021.
!   Copyright 2021 ___ORGANIZATIONNAME___. All rights reserved.
!
!
!   This module contains the interface to python
!   register_couple()
!   register_par()
!   trace_index_list()
!   get_ntrace()
!   reset_couple()
!   compute_xcorr()
!
!
module a1xcor

!use iso_fortran_env
use cxcor_mod       ! configuration of cxcor
use pxcor_mod
implicit none


contains

! =========================================================
! display_cxc
! =========================================================
subroutine display_cxc()
integer ::i,j
write(0,*) 'ntrace=',cxc%ntrace
write(0,*) 'ncouple=',cxc%nx
write(0,*) 'matrix of correlations:'
    if (cxc%ntrace<=10) then
        do i=1,cxc%ntrace
            write(0,*) ('    *',j=1,i-1),(cxc%x(cxc%k(i,j)),j=i,cxc%ntrace)
        enddo
    else
        do i=1,5
            write(0,*) ('    *',j=1,i-1),(cxc%x(cxc%k(i,j)),j=i,5),'   ...',(cxc%x(cxc%k(i,j)),j=cxc%ntrace-5,cxc%ntrace)
        end do
        write(0,*) '                     ...       ...'
        write(0,*) '                     ...       ...'
        do i=cxc%ntrace-5,cxc%ntrace
            !write(0,*) ('    *',j=6,i+15-cxc%ntrace),(cxc%x(cxc%k(i,j)),j=i,5),'   ...', &
             write(0,*) ('    *',j=1,i-cxc%ntrace+10),(cxc%x(cxc%k(i,j)),j=i,5),'   ... ', &
                        (cxc%x(cxc%k(i,j)),j=i,cxc%ntrace)
        end do
    end if
end subroutine display_cxc

! =========================================================
! display_pxc
! =========================================================
subroutine display_pxc()
   call pxc_print_info()
end subroutine display_pxc

! =========================================================
! register_couple
! =========================================================
subroutine register_couple(couple, ntraces, status, two, ncouple)
    integer, intent(in) :: couple(two, ncouple)
    integer, intent(in) :: ntraces
    integer, intent(out) :: status
    integer, intent(in) :: ncouple, two

!f2py intent(in, c) :: couple
!f1py integer(out) :: status
!f2py integer intent(hide), depend(couple) :: two=shape(couple,1)
!f2py integer intent(hide), depend(couple) :: ncouple=shape(couple,0)

    integer :: i

    status =0
    call cxc%init(ntraces)
    do i=1,ncouple
        call cxc%register_ij(couple(1,i), couple(2,i), status)
        if (status >0) exit
    enddo
end subroutine

! =========================================================
! register_par
! =========================================================
subroutine register_par(proc, ier, val, nval)
    character (len=10) :: proc
    integer      :: nval, ier
    real(kind=8) :: val(nval)
!f2py intent(in)    :: proc
!f2py intent(out)   :: ier
!f2py intent(in, c) :: val
!f2py integer intent(hide), depend(val) :: nval=shape(val,0)


    write(0,*) 'processing <',proc,'>'
    write(0,*) 'received ',nval,'parameters'
    write(0,*) 'as ',val

    call pxc_register_proc(proc, val, nval, ier)
end subroutine

! =========================================================
! trace_index_list
!
! retourne la liste des traces utilisees pour le calcul des correlations
! =========================================================
subroutine trace_index_list(list, nx, n)
    integer, intent(out) :: list(n), nx
    integer,intent(in)   :: n
!f2py intent(in, out) :: list
!f2py intent(out) :: nx
!f2py integer intent(hide), depend(list) :: n=shape(list,0)

    call cxc%list_traces(list,nx)

end subroutine

! =========================================================
! get_ntrace
!
! retourne le nombre de traces déclarees
!
! =========================================================
subroutine get_ntrace(ntrace)
!f2py intent(out) :: ntrace
     integer, intent(out) :: ntrace

!write(0,*) 'ntrace =',ntrace
     ntrace = cxc%get_ntrace()
!write(0,*) 'ntrace =',ntrace
end subroutine

! =========================================================
! get_nxcorr
!
! retourne le nombre de correlation à calculer
!
! =========================================================
subroutine get_nxcorr(nxcorr)
!f2py intent(out) :: nxcorr
     integer, intent(out) :: nxcorr

!write(0,*) 'ntrace =',nxcorr
     nxcorr = cxc%get_nxcorr()
!write(0,*) 'ntrace =',nxcorr
end subroutine

! =========================================================
! reset_couple
! =========================================================
subroutine reset_couple()
     call cxc%cleanUp()
end subroutine reset_couple

! =========================================================
! compute_xcorr
!
! interface to python for compute_xcorr_mkl() fortran routine
! input:
!   arrayIn(m1,m2) = section(space_dim, time_dim)  IN FORTRAN, in PYTHON IT IS arrayIn[m2, m1]
!   arrayOut(n1,n2) = xcor(time_dim, nxcor)
!   ijx(nn1,nn2)    = indexing of xcorr second dimension with respect to trace ordering in arrayIn
!   lag = xcorr lag in samples
!   stack = flag to stack FFT of length 2*lag+1 along the time_im dimension
!   m1,m2,n1,n2,nn1,nn2 = input array dimension managed by f2py
! output:
!   ier= error code
! =========================================================
subroutine compute_xcorr(arrayIn, arrayOut, ijx, lag, stack, transposed, ier, m1, m2, n1, n2, nn1, nn2)

    use xcorr_mod
    !use iso_fortran_env
    implicit none

!f2py intent(inout, overwrite, c) :: arrayIn
!f2py intent(in, out, overwrite, c) ::  ijx
!f2py intent(in, out, overwrite, c) :: arrayOut
!f2py intent(in) lag, stack, transposed
!f2py intent(out) ier
!f2py integer intent(hide), depend(arrayIn) :: m1=shape(arrayIn,1)   ! m1 is dim 1 in python
!f2py integer intent(hide), depend(arrayIn) :: m2=shape(arrayIn,0)   ! m2 is dim 0 in python
!f2py integer intent(hide), depend(arrayOut) :: n1=shape(arrayOut,1)
!f2py integer intent(hide), depend(arrayOut) :: n2=shape(arrayOut,0)
!f2py integer intent(hide), depend(ijx)     :: nn1=shape(ijx,1)
!f2py integer intent(hide), depend(ijx)     :: nn2=shape(ijx,0)

    integer, intent(in)    :: m1,m2,n1,n2
    integer, intent(in)    :: lag, stack
    logical, intent(in)    :: transposed
    integer                :: nn1,nn2
    integer, intent(out)   :: ier
    integer, intent(inout) :: ijx(nn1,nn2)

    real(kind=8) :: arrayIn(m1,m2), arrayOut(n1,n2)
    integer      :: ntime,nspace

write(0,*) 'arrayOut(a,b)',n1,n2,'ijx(a,b) (2,*)',nn1,nn2
write(0,*) 'arrayIn(a,b) (space, time)',m1,m2
write(0,*) 'arrayIn',arrayIn(2,1:3)
write(0,*) 'transposed',transposed
    !write(0,*) 'stack',stack
    if (transposed) then
       ntime=m1; nspace=m2
    else
       nspace=m1; ntime=m2
    endif
    call compute_xcorr_core(arrayIn, nspace, ntime, n1, lag, stack, arrayOut, nn2, ijx, ier )


end subroutine
end module a1xcor
