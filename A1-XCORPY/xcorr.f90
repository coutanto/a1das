!
!   xcorr.f90
!   
!
!   Created by coutanto on 11/08/2021.
!
module xcorr_mod

#define MKL

use iso_fortran_env
use cxcor_mod
use pxcor_mod

#ifdef MKL
use mkl_dfti
#endif

implicit none

#ifdef FFTW
include 'fftw_f77.i
#endif

type fft_handle
#ifdef MKL
    type(dfti_descriptor), pointer :: h
#endif
#ifdef FFTW
    integer :: h
#endif
end type fft_handle

contains
include 'whitening.f90'
include 'filtering.f90'

! =========================================================
! compute_xcorr_core
!
! input:
!    arrayIn(nspace, ntime) = bloc de données
!    nspace, ntime = dimension du bloc
!    len_xcor = longueur de la cross-correlation (dépend du lag et de la parité de ntime)
!    lag = lag pour la xcorrelation
!    nx = dimension du tableau ijx(2,nx) et arrayOut(ntime, nx)
! output:
!    xcorr(ntime, *) = cross correlation pour les couples {i,j}
!                         specifies (uniquement pour j>=i)
!    ijx(2,nx) = indice des traces {i,j} dans xcorr
!    nx = nbre de xcorr exactement calculees
!    ier = code de sortie
!          0 ok
!          1 sous dimensionnement de nx
! =========================================================
subroutine compute_xcorr_core(arrayIn, nspace, ntime, len_xcor, lag, stack, array_xcorr, nx, ijx, ier )

use cxcor_mod
use pxcor_mod

real(kind=real64), intent(inout) :: arrayIn(:,:), array_xcorr(len_xcor,nx)
integer, intent(in) :: ntime, nspace, len_xcor, lag, stack
integer, intent(inout) :: nx, ijx(2,*)
integer, intent(out) :: ier
real(kind=real64):: fmin,fmax,r,fs

integer, allocatable :: list_trace(:)
integer :: nt,stats, fft_len, i0
real(kind=real64), allocatable :: taper(:), arrayTmp(:,:)
real(kind=real64) :: perc
complex(kind=real64), allocatable :: arrayFFT(:,:),xspectra(:,:)
type (fft_handle) :: fft_h
integer :: i,l
logical :: transposed
ier=0

!
! test if array_in is (space x time) or transposed into (time x space)
!
if (size(arrayIn,1) == nspace) then
   transposed = .false.
else
   transposed = .true.
endif
!
! test si la dimension de xcorr et ijx est suffisante
!
if (nx <cxc%nx) then
   write(0,*) 'nx=',nx,'cxc%nx=',cxc%nx
   ier=1
   return
else
   nx = cxc%nx
endif
!
! check value of lag time
!
if (lag>ntime/2) then
    write(0,*) 'wrong value for lag time, too large >',ntime/2
    return
end if
!
! recuperation de la liste des traces a TF transformer
!
allocate(list_trace(cxc%ntrace))
call cxc%list_traces(list_trace,nt)


!
! transposition out of place
!
! ?? est-ce qu'on reduit le tableau uniquement aux traces à transformer
! pour avoir une transposition + petite, ou non?
! est-ce que cela est fait dans python?
allocate(arrayTmp(ntime,nspace))
if (transposed) then
    arrayTmp = arrayIn
else
    call transpose(arrayIn, arrayTmp, nspace, ntime)
endif


!
! preprocessing temporel
!
if (pxc%flag_raised("onebit")) then
    arrayTmp=sign(1.d0,arrayTmp)
endif

if (pxc%flag_raised("clip")) then
   ! compute std
   ! and clip with a mask to a factor of std
endif

if (stack == 1) then
    fft_len = len_xcor
else
    fft_len = ntime
end if

if (pxc%flag_raised("taper")) then
   ! compute taper window to be used with fft
   perc=pxc%proc(pxc%index("taper"))%args(1)
   allocate(taper(fft_len))
   call create_taper(fft_len, perc, taper)
endif

allocate(arrayFFT(fft_len,nt))
call init_fft(fft_len, fft_h)


allocate(xspectra(fft_len,nx))
xspectra = 0.d0

do i0 = 1,ntime,fft_len
	
!
! compute FFT
!
    call compute_FFT(arrayTmp,arrayFFT, fft_len, i0, list_trace, nt, fft_h, taper)
!
! processing spectraux
!
    if (pxc%flag_raised("white")) then
		fmin=pxc%proc(pxc%index("white"))%args(1)
		fmax=pxc%proc(pxc%index("white"))%args(2)
      	call whitening(arrayFFT,fft_len,nt,fmin,fmax)
    endif
    if (pxc%flag_raised("bandpass")) then
		fmin=pxc%proc(pxc%index("bandpass"))%args(1)
		fmax=pxc%proc(pxc%index("bandpass"))%args(2)
      	call bandpass(arrayFFT,fft_len,nt,fmin,fmax)
    endif

    if (pxc%flag_raised("lowpass")) then
		fmax=pxc%proc(pxc%index("lowpass"))%args(1)
      	call lowpass(arrayFFT, fft_len, nt, fmax)
    endif

!
! compute cross spectra
!
    call compute_cross_spectra(arrayFFT, xspectra, fft_len, nt, cxc%nx, ijx, ier)

end do
!
! back in time, et correlation
!
call compute_cross_corr(xspectra, array_xcorr, fft_len, cxc%nx, len_xcor, lag, fft_h, ier)

deallocate(arrayTmp)
deallocate(arrayFFT)
deallocate(xspectra)
if (allocated(taper)) deallocate(taper)

call release_fft(fft_h)

end subroutine


! =========================================================
! transpose
!
! out of place transposition
! =========================================================
subroutine transpose(arrayi, arrayo, nrow, ncol)

    include 'mkl.fi'

    real(kind=real64),intent(in) :: arrayi(nrow,ncol)
    real(kind=real64),intent(out) :: arrayo(ncol,nrow)
    integer,intent(in)              :: nrow,ncol
    real(kind=real64)               :: alpha=1.d0

    call mkl_domatcopy('C','T',nrow,ncol,alpha,arrayi,nrow,arrayo,ncol)
end subroutine

! =========================================================
! INIT_FFT
! =========================================================
subroutine init_fft(nfft, fft_h)
    use mkl_dfti
    integer , intent(in) :: nfft
    type(fft_handle), intent(inout) :: fft_h

#ifdef MKL
    integer :: stats

    stats = DftiCreateDescriptor(fft_h%h, DFTI_DOUBLE, DFTI_REAL, 1, nfft)
    stats = DftiSetValue(fft_h%h, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    stats = DftiSetValue(fft_h%h, DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX)
    stats = DftiSetValue(fft_h%h, DFTI_PACKED_FORMAT, DFTI_CCE_FORMAT);
    stats = DftiSetValue(fft_h%h, DFTI_INPUT_STRIDES,1);
    stats = DftiSetValue(fft_h%h, DFTI_OUTPUT_STRIDES, 1);
    stats = DftiSetValue(fft_h%h, DFTI_BACKWARD_SCALE, 1.d0/real(nfft))
    stats = DftiCommitDescriptor(fft_h%h)

#endif
#ifdef FFTW
    call create_fft_plan_my_c_wrapper(fft_h%h,???)
#endif
end subroutine init_fft

! =========================================================
! RELEASE_FFT
! =========================================================
subroutine release_fft(fft_h)
    type(fft_handle) :: fft_h

    integer :: stats
#ifdef MKL
    stats = DftiFreeDescriptor (fft_h%h)
#endif
#ifdef FFTW
    call destroy_fft_plan_my_c_wrapper(fft_h%h)
#endif
end subroutine release_fft
! =========================================================
! compute_fft
! =========================================================
subroutine compute_FFT(array,c_array, nfft, i0, list, ntrace, fft_h, taper)
use mkl_dfti

real(kind=real64),intent(inout)   :: array(:,:)
real(kind=real64),intent(in), allocatable   :: taper(:)
complex(kind=real64),intent(out):: c_array(nfft,ntrace)
integer , intent(in) :: nfft, i0,ntrace,list(ntrace)
type(fft_handle), intent(in) :: fft_h

#ifdef MKL
integer :: stats, i, ii

!$OMP PARALLEL default(shared), private(i,ii,stats)
!$OMP DO
do i=1,ntrace
  ii=list(i)
  if (allocated(taper)) array(i0:i0+nfft-1,ii) = array(i0:i0+nfft-1,ii)*taper
  Stats = DftiComputeForward (fft_h%h, array(i0:i0+nfft-1,ii), c_array(:,i))
enddo
!$OMP END PARALLEL
#endif

#ifdef FFTW
    call fftw_multiple_my_c_wrapper(fft_h%h,)
#endif

end subroutine

! =========================================================
! compute_cross_spectra
!
! input:
!    arrayFFT = array of nt FFT of length nfft (nfft x nt)
!    xspectra = array of nx cross-spectra (nfft x nx)
!    nfft = fft length
!    nt = number of traces in FFT array
!    nx = number of cross-spectra to compute
! output:
!    ijx = array of cxc%nx cross-correlation
!    status = 0 if ok, 2, 3 if errors
!
! global:
!    cxc and pxc class/structure instance
! =========================================================
!
! Rule to determine the cross-correlation length depending on lag <lag> and parity of signal length <ntime>
!
!                    xcor_length    lag            mkl_fft_shift        copy_indices
!  ntime even(2*k)
!   lag = 0             ntime       ntime/2           ntime/2 (=lag)      [1:ntime]
!   lag # 0            2*lag+1        --              lag                 [1:2*lag+1]
! ntime odd (2*k+1)
!   lag = 0              ntime      (ntime-1)/2      (ntime-1)/2(=lag)    [1:ntime](=[1:2*lag+1])
!   lag # 0             2*lag+1        --            lag                  [1:2*lag+1]
!
!    cross-correlation zero lag is thus always at position <lag> in the xcorr array (xcor[lag-1])
!
!
subroutine compute_cross_spectra(arrayFFT, xspectra, nfft, nt, nx, ijx, status)

complex(kind=real64),intent(in)  :: arrayFFT(nfft,nt)
complex(kind=real64),intent(inout)  :: xspectra(nfft,nx)
integer, intent(inout) :: status
integer, intent(in)  :: nfft, nt, nx
integer, intent(out) :: ijx(2,nx)


integer :: i,j,k,l,stats, i_index, j_index
integer, allocatable :: ii(:),jj(:)
complex(kind=real64),allocatable :: cxf(:)
real(kind=real64),allocatable :: cxt(:)

status = 0
!
! count and set indices shared by thread
!
allocate(ii(nx),jj(nx))
l=1
do i=1,cxc%ntrace
  do j=i,cxc%ntrace
     k = cxc%k(i,j)
     if (cxc%x(k) /= 0) then
       ii(l) = i
       jj(l) = j
       ijx(1,l) = i
       ijx(2,l) = j
       l=l+1
     endif
  enddo
enddo

! check number of xcorr
if (l-1 /= nx) then
   write(0,*) 'compute_correlation: incoherent number of xcorr between nx and l',nx,l-1
   status = 1
   return
endif

! check indexing
do l=1,nx
    i_index = cxc%index_in_list(ii(l))
    j_index = cxc%index_in_list(jj(l))
    if (i_index==0 .or. j_index==0) then
        write(0,*) 'compute_correlation: Abort: bug: trace i ou j non enregistre pour la fft'
        status = 2
        return
    endif
end do

! compute cross-spectra
!$OMP PARALLEL default(shared), firstprivate(l, i_index, j_index)
!$OMP DO
do l=1,nx  ! loop over all x-corr
  ! get the index in the list of FFTs for trace_i and trace_j
  i_index = cxc%index_in_list(ii(l))
  j_index = cxc%index_in_list(jj(l))

  xspectra(:,l) = xspectra(:,l) + arrayFFT(:,i_index)*conjg(arrayFFT(:,j_index))

enddo
!$OMP END PARALLEL

deallocate(ii,jj)
end subroutine

! =========================================================!
! compute cross-correlation from cross-spectra by inverse Fourier transform + time shift
! input:
!    xspectra = array of nx cross-spectra (nfft x nx)
!    xcorr = array of cross-correlation (len_xcor x nx)
!    nfft = fft length
!    len_xcor = cross-correlation length
!    nx = number of cross-spectra to compute
!    h = handle on FFT MKL parameters
    !
subroutine compute_cross_corr(xspectra, xcorr, nfft, nx, len_xcor, lag, fft_h, status)
use mkl_dfti
integer, intent(in) :: nx, len_xcor, lag, nfft, status
complex(kind=real64),intent(in)  :: xspectra(nfft,nx)
real(kind=real64),intent(out) :: xcorr(len_xcor,nx)
type(fft_handle), intent(in) :: fft_h

integer :: l, Stats
real(kind=real64), allocatable :: cxt(:)

#ifdef MKL
!$OMP PARALLEL default(shared), firstprivate(l, Stats, cxt)
allocate(cxt(nfft))
!$OMP DO
do l=1,nx
    Stats = DftiComputeBackward (fft_h%h, xspectra(:,l), cxt)
    ! zero lag is at index 1, shift to the right to have it at index
    ! lag+1 consdering that we keep 2*lag+1 samples
    ! cshift can yield a stack overflow for large arraies
    ! see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=12533
    !cxt = cshift(cxt,-lag)
    !xcorr(1:len_xcor,l) = cxt(1:len_xcor)
    xcorr(1:lag,l) = cxt(nfft-lag+1:nfft)
    xcorr(lag+1:len_xcor,l) = cxt(1:len_xcor-lag)
end do
deallocate(cxt)
!$OMP END PARALLEL
#endif

end subroutine compute_cross_corr

subroutine create_taper(len, perc, taper)
implicit none
integer, intent(in) :: len
real(kind=real64), intent(in) :: perc
real(kind=real64), intent(inout), allocatable :: taper(:)

#define pi 3.14159265d0
integer :: i,l
l = floor(len*perc)
taper=1.d0
do i=1,l
  taper(i) = (1. - cos(pi*(i-1)/l))/2.0
  taper(len-l+i) = (1. - cos(pi*(l-i)/l))/2.0
enddo

end subroutine
end module xcorr_mod

