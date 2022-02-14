#ifdef MAIN
program main
implicit none

integer, parameter :: nt=1000000,ns=96
real(kind=8) :: sos(6,3),g
real(kind=8),allocatable :: sig(:,:)
integer :: i

allocate(sig(nt,ns))

do i=1,3
sos(1,i)=1.d0
sos(2,i)=-2.d0
sos(3,i)=1.d0
sos(4,i)=1.d0
enddo
sos(1,1)=0.0295882236386606d0
sos(2,1)=-0.0591764472773213d0
sos(3,1)=0.0295882236386606d0


sos(5,1)=1.87350135405495d-16
sos(5,2)=1.06858966120171d-15
sos(5,3)=-4.23272528138341d-16
sos(6,1)=0.0173323801209992d0
sos(6,2)=0.17157287525381d0
sos(6,3)=0.588790706480864d0

g=0.0295882236386606d0

sig=0.d0
do i=1,ns
sig(1,i)=1.d0
enddo

call sosfilt(sos,3,sig,nt,ns)

write(6,*) sig(1:10,2)



end program
#endif


! ===============================  Double precision version ====================
! sos filtering
! reference: https://en.wikipedia.org/wiki/Digital_biquad_filter, matlab zp2sos
!
! The filter is defined as a cascade of digital biquad filters
! H(z) = (b0 + b1*z^-1 + b2*z^-2) / (1 + a1*z^-1 + a2*z^-2)
!
! given by the SOS matrix (second order section)
!
! | b0 b1 b2 1 a1 a2 |  1st filter
! | b0 b1 b2 1 a1 a2 |  2nd filter
! | b0 b1 b2 1 a1 a2 |  ...
!
! The number of biquad filters is the lowest integer larger or equal than (npoles/2, nzeros/2)
! A sixth order butterworh filter yields thus a 3 x 6 sos matrix
!
! After different testings the formula below is slightly slower than 
! the very standard formula
! w(n) =    x(n) - a1*w(n-1) + a2*w(n-2)
! y(n) = b0*w(n) + b1*w(n-1) + b2*w(n-2)
! standard formula:
! y(n) = b0x(n)+b1*x(n-1)+b2*x(n-2)-a1*y(n-1)-a2*y(n-2)
!
! after different testings, the fastest algo is 1)
! 1) each thread write/read directly from sig(:,:) shared array
! 2) copy a small part of sig(:,:) in the stack would be faster but
!    the stack isn't large enough to hold long time series (nt<10000)
! 3) copy a small part of sig(:,:) in allocated memory for each thread is fine
!    for any time length but is slower than 1)
!
subroutine sosfilt_d(sos,sig,j1,j2,jj,order,ncof,nt,ns)
   implicit none
   integer, intent(in)        :: j1, j2, jj    ! indices for second dimension of sig, 0 based
   integer,intent(in)         :: order         ! filter order
   integer,intent(in)         :: ncof          ! should be 6
   real(kind=8),intent(in)    :: sos(ncof,order)  ! sos matrix
   integer,intent(in)         :: nt,ns         ! length of signal, number of signal
   real(kind=8),intent(inout) :: sig(nt,ns)      ! signals

!f2py intent(in) j1,j2,jj
!f2py intent(in,out,overwrite,c) sig
!f2py intent(in,c) sos
!f2py intent(hide),depend(sig) :: ns=shape(sig,0), nt=shape(sig,1)
!f2py intent(hide),depend(sos) :: order=shape(sos,0), ncof=shape(sos,1)

   integer :: o,i,j
   real(kind=8) :: a1,a2,b0,b1,b2
   integer :: omp_get_thread_num,omp_get_num_threads

   real(kind=8) :: x0,xm1,xm2,y0,ym1,ym2

   if (ncof /= 6) then
       write(0,*) 'number of sos coefficient should be 6 not ',ncof
       stop
   endif

!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) &
!$OMP SHARED(sig,sos)
#if defined(_OPENMP)
if (omp_get_thread_num()==1) then
!  write(0,*) 'number of threads ',omp_get_num_threads()
! write(0,*) 'omp stack size ',omp_get_stacksize_s()
endif
#endif


#if defined(_OPENMP)
!write(0,*) 'allocation thread ',omp_get_thread_num(),' is ',ok
#endif



!$OMP DO ORDERED, SCHEDULE(DYNAMIC)
!   do j=1,ns             ! loop over signals
   do j=j1+1,j2,jj        ! loop over time
#if defined(_OPENMP)
!    write(0,*) 'running on thread ',omp_get_thread_num()
#endif

     do o=1,order            ! loop over filter cascade
       b0=sos(1,o);b1=sos(2,o);b2=sos(3,o)
       a1=sos(5,o);a2=sos(6,o)
       ! initialize recursion, point 1
       x0 = sig(1,j)
       y0 = b0*x0
       sig(1,j) = y0
       xm1=0.d0
       xm2=0.d0
       ym1=0.d0
       ym2=0.d0
       !  ------             , point 2
       xm1=x0
       x0=sig(2,j)
       ym1=y0
       y0 = b0*x0 + b1*xm1 - a1*ym1
       sig(2,j) = y0
       do i=3,nt             ! loop over time
         xm2 = xm1
         xm1 = x0
         x0 = sig(i,j)
         ym2 = ym1
         ym1 = y0
         y0 = b0*x0 + b1*xm1 +b2*xm2 -a1*ym1 -a2*ym2
         sig(i,j)=y0
       enddo
     enddo
   enddo

!$OMP END PARALLEL
end subroutine


subroutine sosfiltfilt_d(sos,sig,j1,j2,jj,ncof,order,nt,ns)
   implicit none
   integer, intent(in)        :: j1, j2, jj    ! indices for second dimension of sig, 0 based
   integer,intent(in)         :: order         ! filter order
   integer,intent(in)         :: ncof          ! should be 6
   real(kind=8),intent(in)    :: sos(ncof,order)  ! sos matrix
   integer,intent(in)         :: nt,ns         ! length of signal, number of signal
   real(kind=8),intent(inout) :: sig(nt,ns)      ! signals

!f2py intent(in) j1,j2,jj
!f2py intent(in,out,overwrite,c) sig
!f2py intent(in,c) sos
!f2py integer,intent(hide),depend(sig) :: ns=shape(sig,0), nt=shape(sig,1)
!f2py integer,intent(hide),depend(sos) :: order=shape(sos,0), ncof=shape(sos,1)


   integer :: o,i,j
   real(kind=8) :: a1,a2,b0,b1,b2
   integer :: omp_get_thread_num,omp_get_num_threads

   real(kind=8) :: x0,xm1,xm2,y0,ym1,ym2

   if (ncof /= 6) then
       write(0,*) 'number of sos coefficient should be 6 not ',ncof
       stop
   endif

!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) &
!$OMP SHARED(sig,sos)
#if defined(_OPENMP)
if (omp_get_thread_num()==1) then
!  write(0,*) 'number of threads ',omp_get_num_threads()
! write(0,*) 'omp stack size ',omp_get_stacksize_s()
endif
#endif

#if defined(_OPENMP)
!write(0,*) 'allocation thread ',omp_get_thread_num(),' is ',ok
#endif
!$OMP DO ORDERED, SCHEDULE(DYNAMIC)
!   do j=1,ns             ! loop over signals
   do j=j1+1,j2,jj
#if defined(_OPENMP)
!    write(0,*) 'running on thread ',omp_get_thread_num()
#endif

     do o=1,order            ! loop over filter cascade
       b0=sos(1,o);b1=sos(2,o);b2=sos(3,o)
       a1=sos(5,o);a2=sos(6,o)
       ! initialize recursion, point 1
       x0 = sig(1,j)
       y0 = b0*x0
       sig(1,j) = y0
       xm1=0.d0
       xm2=0.d0
       ym1=0.d0
       ym2=0.d0
       !  ------             , point 2
       xm1=x0
       x0=sig(2,j)
       ym1=y0
       y0 = b0*x0 + b1*xm1 - a1*ym1
       sig(2,j) = y0
       do i=3,nt             ! loop over time
         xm2 = xm1
         xm1 = x0
         x0 = sig(i,j)
         ym2 = ym1
         ym1 = y0
         y0 = b0*x0 + b1*xm1 +b2*xm2 -a1*ym1 -a2*ym2
         sig(i,j) = y0
       enddo
     enddo

     do o=1,order            ! loop over filter cascade 
       b0=sos(1,o);b1=sos(2,o);b2=sos(3,o)
       a1=sos(5,o);a2=sos(6,o)
       ! initialize recursion, point 1
       x0 = sig(nt,j)
       y0 = b0*x0
       sig(nt,j) = y0
       xm1=0.d0
       xm2=0.d0
       ym1=0.d0
       ym2=0.d0
       !  ------             , point 2
       xm1=x0
       x0=sig(nt-1,j)
       ym1=y0
       y0 = b0*x0 + b1*xm1 - a1*ym1
       sig(nt-1,j) = y0
       do i=nt-2,1,-1             ! loop over time backward
         xm2 = xm1
         xm1 = x0
         x0 = sig(i,j)
         ym2 = ym1
         ym1 = y0
         y0 = b0*x0 + b1*xm1 +b2*xm2 -a1*ym1 -a2*ym2
         sig(i,j) = y0
       enddo
     enddo     

   enddo


!$OMP END PARALLEL
end subroutine

! ===============================  Single precision version ====================
!
! signal is single precision but computation are performed in double precision
!

subroutine sosfilt_s(sos,sig,j1,j2,jj,order,ncof,nt,ns)
   implicit none
   integer, intent(in)        :: j1, j2, jj    ! indices for second dimension of sig, 0 based
   integer,intent(in)         :: order         ! filter order
   integer,intent(in)         :: ncof          ! should be 6
   real(kind=8),intent(in)    :: sos(ncof,order)  ! sos matrix
   integer,intent(in)         :: nt,ns         ! length of signal, number of signal
   real(kind=4),intent(inout) :: sig(nt,ns)      ! signals

!f2py intent(in) j1,j2,jj
!f2py intent(in,out,overwrite,c) sig
!f2py intent(in,c) sos
!f2py intent(hide),depend(sig) :: ns=shape(sig,0), nt=shape(sig,1)
!f2py intent(hide),depend(sos) :: order=shape(sos,0), ncof=shape(sos,1)

   integer :: o,i,j
   real(kind=8) :: a1,a2,b0,b1,b2
   integer :: omp_get_thread_num,omp_get_num_threads
   logical :: isnan

   real(kind=8) :: x0,xm1,xm2,y0,ym1,ym2

   if (ncof /= 6) then
       write(0,*) 'number of sos coefficient should be 6 not ',ncof
       stop
   endif

!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) &
!$OMP SHARED(sig,sos)
#if defined(_OPENMP)
if (omp_get_thread_num()==1) then
  !write(0,*) 'number of threads ',omp_get_num_threads()
 !write(0,*) 'omp stack size ',omp_get_stacksize_s()
endif
#endif


#if defined(_OPENMP)
!write(0,*) 'allocation thread ',omp_get_thread_num(),' is ',ok
#endif

!$OMP DO ORDERED, SCHEDULE(DYNAMIC)
!   do j=1,ns             ! loop over signals
   do j=j1+1,j2,jj
#if defined(_OPENMP)
!    write(0,*) 'running on thread ',omp_get_thread_num()
#endif

     do o=1,order            ! loop over filter cascade
       b0=sos(1,o);b1=sos(2,o);b2=sos(3,o)
       a1=sos(5,o);a2=sos(6,o)
       ! initialize recursion, point 1
       x0 = sig(1,j)
       y0 = b0*x0
       sig(1,j) = y0
       xm1=0.d0
       xm2=0.d0
       ym1=0.d0
       ym2=0.d0
       !  ------             , point 2
       xm1=x0
       x0=sig(2,j)
       ym1=y0
       y0 = b0*x0 + b1*xm1 - a1*ym1
       sig(2,j) = y0
       do i=3,nt             ! loop over time
         xm2 = xm1
         xm1 = x0
         x0 = sig(i,j)
         !where (ieee_is_nan(x0)) x0 = 0.d0 !better move that to python
         ym2 = ym1
         ym1 = y0
         y0 = b0*x0 + b1*xm1 +b2*xm2 -a1*ym1 -a2*ym2
         sig(i,j)=y0
       enddo
     enddo
   enddo

!$OMP END PARALLEL
end subroutine



subroutine sosfiltfilt_s(sos,sig,j1,j2,jj,ncof,order,nt,ns)
   implicit none
   integer, intent(in)        :: j1, j2, jj    ! indices for second dimension of sig, 0 based
   integer,intent(in)         :: order         ! filter order
   integer,intent(in)         :: ncof          ! should be 6
   real(kind=8),intent(in)    :: sos(ncof,order)  ! sos matrix
   integer,intent(in)         :: nt,ns         ! length of signal, number of signal
   real(kind=4),intent(inout) :: sig(nt,ns)      ! signals

!f2py intent(in)                 :: j1,j2,jj
!f2py real(kind=4), intent(in,out,overwrite,C) :: sig
!f2py intent(in,c)               :: sos
!f2py integer intent(hide),depend(sig) :: ns=shape(sig,0), nt=shape(sig,1)
!f2py integer intent(hide),depend(sos) :: order=shape(sos,0), ncof=shape(sos,1)


   integer :: o,i,j
   real(kind=8) :: a1,a2,b0,b1,b2
   integer :: omp_get_thread_num,omp_get_num_threads

   real(kind=8) :: x0,xm1,xm2,y0,ym1,ym2
   real(kind=8),allocatable :: siglocal(:)

   !write(0,*) 'sosfilt 2',sos,j1,j2,jj,order,ncof,nt,ns
   if (ncof /= 6) then
       write(0,*) 'number of sos coefficient should be 6 not ',ncof
       stop
   endif

!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) &
!$OMP SHARED(sig,sos)
#if defined(_OPENMP)
if (omp_get_thread_num()==1) then
!  write(0,*) 'number of threads ',omp_get_num_threads()
! write(0,*) 'omp stack size ',omp_get_stacksize_s()
endif
#endif
    allocate(siglocal(nt))
#if defined(_OPENMP)
!write(0,*) 'allocation thread ',omp_get_thread_num(),' is ',ok
#endif
!$OMP DO ORDERED, SCHEDULE(DYNAMIC)
!   do j=1,ns             ! loop over chunk of signals
   do j=j1+1,j2,jj
#if defined(_OPENMP)
!    write(0,*) 'running on thread ',omp_get_thread_num()
#endif
     ! copy of data in double precision in the thread stack
     ! it may be necessary to increase thread stack using OMP_THREADSTACK
     ! env variable
     siglocal = sig(:,j)

     !
     ! forward filtering
     !
     do o=1,order            ! loop over filter cascade
       b0=sos(1,o);b1=sos(2,o);b2=sos(3,o)
       a1=sos(5,o);a2=sos(6,o)
       ! initialize recursion, point 1
       x0 = siglocal(1)
       y0 = b0*x0
       siglocal(1) = y0
       xm1=0.d0
       xm2=0.d0
       ym1=0.d0
       ym2=0.d0
       !  ------             , point 2
       xm1=x0
       x0 = siglocal(2)
       ym1=y0
       y0 = b0*x0 + b1*xm1 - a1*ym1
       siglocal(2)=y0
       do i=3,nt             ! loop over time
         xm2 = xm1
         xm1 = x0
         x0 = siglocal(i)
         ym2 = ym1
         ym1 = y0
         y0 = b0*x0 + b1*xm1 +b2*xm2 -a1*ym1 -a2*ym2
         siglocal(i) = y0
       enddo
     enddo

     !
     ! backward filtering
     !
     do o=1,order            ! loop over filter cascade
       b0=sos(1,o);b1=sos(2,o);b2=sos(3,o)
       a1=sos(5,o);a2=sos(6,o)
       ! initialize recursion, point 1
       x0=siglocal(nt)
       y0 = b0*x0
       siglocal(nt) = y0
       xm1=0.d0
       xm2=0.d0
       ym1=0.d0
       ym2=0.d0
       !  ------             , point 2
       xm1=x0
       x0 = siglocal(nt-1)
       ym1=y0
       y0 = b0*x0 + b1*xm1 - a1*ym1
       siglocal(nt-1) = y0
       do i=nt-2,1,-1             ! loop over time backward
         xm2 = xm1
         xm1 = x0
         x0 = siglocal(i)
         ym2 = ym1
         ym1 = y0
         y0 = b0*x0 + b1*xm1 +b2*xm2 -a1*ym1 -a2*ym2
         siglocal(i) = y0
       enddo
     enddo
     sig(:,j) = siglocal
   enddo

   deallocate(siglocal)
!$OMP END PARALLEL
end subroutine

