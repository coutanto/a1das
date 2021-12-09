!
!   cxcorr.f90
!   
!
!   Created by coutanto on 06/08/2021.
!
!   This module describe the CXC class that
!   manages the list of couple (i,j) for which one decide to compute
!   cross correlation
!
! compatibility: fortran 2003 and above
!
module cxcor_mod

!
implicit none

integer(kind=1),parameter :: one=1, two=2
!
! ================= members =================
!
private
public :: cxcor

!
! ================= class definition =================
!

type cxcor

    ! X:  array of cross-correlation indices
    !     It contains the index of a given {i,j} xcorr
    !     in the list of all the computed xcorr
    !
    ! X is a triangular superior matrix stored as a vector
    !   where the index of a given (i,j) couple is stored at
    !  location k(i,j) = i + n*(j-1) - j*(j-1)/2
    !
    ! Ex for 3 traces, X =     x  x  x    k=   1  2  3
    !                             x  x            4  5
    !                                x               6
    !
    !  X(l) with l=k(i,j) is the index in the cross correlation list
    !  of the cross-correlation Cx(Sig_i, Sig_j)
    ! if i>j, we replace by (j,i) and it is the responsability of the
    ! calling program to flip the result along time
    integer(kind=1), allocatable :: X(:)

    ! number of correlation to be performed (or # of non zero entry in X)
    integer :: nx

    ! number of trace, X is ntrace x ntrace
    integer :: ntrace

    ! - list of trace for which one must compute fft
    ! - number of trace
    ! - index in list_trace() for each trace, or 0 if not used
    integer, pointer :: list_trace(:)
    integer :: nt
    integer, pointer :: index_in_list(:)
    
contains
   !
   ! ================= methods =================
   !

    ! register couple or range of couple for which one compute xcorr
    procedure :: register_ij, register_irange
    ! remove couple or range of couple for which one compute xcorr
    procedure :: remove_ij, remove_irange
    ! get number of xcorr to be computed
    procedure :: get_nxcorr
    ! test for couple xcorr to be computed
    procedure :: is_ij_cxcorr, is_i_cxcorr
    ! obtain list of couple given one trace
    procedure :: list_j_cxcorr
    ! obtain list of unique trace used for correlation
    ! usefull to determine the FFT to be computed
    ! or to select initial traces to be kept
    procedure :: list_traces
    !
    procedure :: cleanUp
    !
    procedure :: init
    !
    procedure :: get_ntrace
    !
    procedure :: k

    ! overload register procedure
    generic :: register => register_ij, register_irange
    generic :: remove => remove_ij, remove_irange
#ifdef REMOVE
    ! class destructor
    final :: cleanUp_
#endif
end type cxcor

#ifdef REMOVE
class constructor
interface cxcor
   module procedure Constructor
end interface
#endif

type(cxcor),public :: cxc
   !
   ! ================= contains =================
   !
contains

#ifdef REMOVE
! ===================================
! Constructor
!
! class constructor
! !! noter que <c> n'a pas besoin d'etre declare
! ===================================
function Constructor() result(c)
   type(cxcor) :: c
   !integer, intent(in) :: ntrace

   !allocate(c%X(ntrace,ntrace))
   c%ntrace=0
   c%nx=0
end function
#endif
! ===================================
! init
!
! class constructor
! !! noter que <c> n'a pas besoin d'etre declare
! ===================================
subroutine init(c, ntrace)
   class(cxcor) :: c
   integer, intent(in) :: ntrace

   if (allocated(c%X)) deallocate(c%X)
   allocate(c%X(ntrace*(ntrace+1)/2))
   c%X=0
   c%ntrace=ntrace
   c%nx=0
   c%nt=0
   nullify(c%list_trace)
   allocate(c%index_in_list(ntrace))
end subroutine

! ===================================
! k(i,j)
!
! compute the index of the location of (i,j) in X
! ===================================
function k(c,i,j)
   class(cxcor), intent(inout) :: c
   integer, intent(in) :: i,j
   integer :: k

   if (i>=j) then
      k = i + c%ntrace*(j-1) - j*(j-1)/2
   else
      k = j + c%ntrace*(i-1) - i*(i-1)/2
   endif

end function
! ===================================
! register_ij
!
! add the couple (i,j) to the list of xcorr to be computed
! WARNING: register_ij expect indices to be on base 1 (starting at 1, not 0)
! ===================================
subroutine register_ij(c,i,j, status)
   class(cxcor), intent(inout) :: c
   integer, intent(inout) :: status
   integer, intent(in) :: i,j

   integer :: k
   status = 0

   if (i > c%ntrace .or. j > c%ntrace) status = 1
   if (i <1 .or. j <1) status = 1
   if (status>0) return

   if (c%x(c%k(i,j))==0) c%nx = c%nx+1

   k = c%k(i,j)
   ! set the value of X as a 2 bit number: 00 (=0) no correlation
   !                                       01 (=1) correlate (i,j) where i>j
   !                                       10 (=2) correlate (j,i) where j>i
   !                                       11 (=3) correlate (i,j) and (j,i)
   if (i>j) then
        c%x(k) = or(c%x(k),one)
   else
       c%x(k) = or(c%x(k),two)
   end if
   
end subroutine

! ===================================
! register_irange
!
! add the couples (i,j), ...(i,k) to the list of xcorr to be computed
! ===================================
subroutine register_irange(c,i,j,k,status)

   class(cxcor), intent(inout) :: c
   integer, intent(inout) :: i,j,k,status

   integer :: jj

   status = 0

   if (k<j) then
     jj=j
     j=k
     k=jj
   endif
   if (i > c%ntrace .or. j > c%ntrace .or. k>c%ntrace) status = 1
   if (i <1 .or. j <1 .or. k<1) status = 1
   if (status>0) return

   do jj=j,k
       if (c%k(i,jj) == 0) c%nx=c%nx + 1
      if (i>jj) then
            c%x(c%k(i,jj)) = or(c%x(c%k(i,jj)),one)
      else
            c%x(c%k(i,jj)) = or(c%x(c%k(i,jj)),two)
      end if
   enddo

end subroutine

! ===================================
! remove_ij
!
! remove the couple (i,j) from the list of xcorr to be computed
! ===================================
subroutine remove_ij(c,i,j)

   class(cxcor), intent(inout) :: c
   integer, intent(inout) :: i,j

   c%x(c%k(i,j)) = 0
   c%nx = c%nx-1

   deallocate(c%list_trace)
   nullify(c%list_trace)
   c%nt=0
end subroutine

! ===================================
! remove_irange
!
! remove the couples (i,j), ...(i,k) to the list of xcorr to be computed
! ===================================
subroutine remove_irange(c,i,j,k)

   class(cxcor), intent(inout) :: c
   integer, intent(inout) :: i,j,k

   integer :: jj

   if (k<j) then
     jj=j
     j=k
     k=jj
   endif
   do jj=j,k
     c%x(c%k(i,jj))=0
   enddo
   c%nx=c%nx - j+i-1

   deallocate(c%list_trace)
   nullify(c%list_trace)
   c%nt=0
end subroutine


! ===================================
! get_nxcorr
!
! return the number of xcorr to be computed
! ===================================
function get_nxcorr(c) result(nx)
   class(cxcor), intent(inout) :: c

   integer nx
   nx=c%nx
end function

! ===================================
! get_ntrace
!
! return the number of traces declared
! ===================================
function get_ntrace(c) result(ntrace)
   class(cxcor), intent(inout) :: c

   integer ntrace
   ntrace = c%ntrace
end function
! ===================================
! is_ij_cxcorr
!
! return true if the xcorr of couple (i,j) has to be computed
! ===================================
function is_ij_cxcorr(c,i,j) result(t)
   class(cxcor), intent(inout) :: c
   integer, intent(in)  :: i,j

   logical  :: t

   t=.false.
   if (c%x(c%k(i,j)) /= 0) t=.true.
end function

! ===================================
! is_i_cxcorr
!
! return true if the xcorr of couple(s) (i,?) have to be computed
! or in other words, if trace i is used or not
! ===================================
function is_i_cxcorr(c,i) result(t)
   class(cxcor), intent(inout) :: c
   integer, intent(in)  :: i
   logical  :: t

   integer :: j
   t=.false.

   do j=1,c%ntrace
     if (c%x(c%k(i,j)) == 0) cycle
     t=.true.
     return
   enddo
end function

! ===================================
! list_j_cxcorr
!
! return the list of indices j for which one must compute (i,j) xcorr
! ===================================
subroutine list_j_cxcorr(c,i,jlist,nj)
   class(cxcor), intent(inout) :: c
   integer, intent(in)  :: i
   integer, intent(inout) :: jlist(*)
   integer, intent(out) :: nj

   integer :: j
   nj=0

   do j=1,c%ntrace
     if (c%x(c%k(i,j)) == 0) cycle
     nj=nj+1
     jlist(nj)=j
   enddo
end subroutine

! ===================================
! list_traces
!
! return the sorted list of indices j for which one must compute a fft
!
! This list is stored in CXC class as <list_trace>
! together with the <index_in_list> array
! ===================================
subroutine list_traces(c,list,nt)
   class(cxcor), intent(inout) :: c
   integer, intent(inout),target :: list(*)
   integer, intent(out) :: nt

   integer :: i,j
   integer, allocatable :: tlist(:)

   ! test si on a deja fait ce travail, alors on le recupere
   if (associated(c%list_trace)) then
       nt=c%nt
       list(1:nt) = c%list_trace(1:nt)
       return
   endif

   !
   ! temporary array to keep track of trace_i, trace_j
   ! used in xcorrs
   !
   allocate (tlist(c%ntrace))
   tlist=0
   do i=1,c%ntrace
      do j=1,c%ntrace
         if (c%x(c%k(i,j)) == 0) cycle
         tlist(j)=1
         tlist(i)=1
      enddo
   enddo

   !
   ! count and sort the list of trace used in xcorrs
   !
   nt=0
   do i=1,c%ntrace
      if (tlist(i)==0) cycle
      nt=nt+1
      list(nt)=i
      if (nt==1) cycle
      if (list(nt)<list(nt-1)) then
          j=list(nt-1)
          list(nt-1)=list(nt)
          list(nt)=j
      endif
   enddo
   deallocate(tlist)

   !
   ! store list_trace array in CXC class
   ! and create the reverse array index_in_list
   allocate(c%list_trace(nt))
   allocate(c%index_in_list(c%ntrace))
   c%index_in_list = 0
   c%list_trace(1:nt) = list(1:nt)
   c%nt = nt
   do i=1,nt
     c%index_in_list(c%list_trace(i))=i
   enddo

end subroutine

! ===================================
! cleanUp: destructor
!
! ===================================
subroutine cleanUp(c)
   class(cxcor) :: c  ! pourquoi faut-il type et non class ??
   deallocate(c%x)
   c%nx=0
   if (associated(c%index_in_list)) deallocate(c%index_in_list)
   if (associated(c%list_trace)) deallocate(c%list_trace)
   c%nt=0
   c%ntrace=0
end subroutine

subroutine cleanUp_(c)
   type(cxcor) :: c  ! pourquoi faut-il type et non class ??
   deallocate(c%x)
   c%nx=0
   if (associated(c%index_in_list)) deallocate(c%index_in_list)
   if (associated(c%list_trace)) deallocate(c%list_trace)
   c%nt=0
   c%ntrace=0
end subroutine



end module cxcor_mod
