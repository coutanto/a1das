!   
!
!   Created by nziengud on 08/09/2021.
!
!



! =========================================================
!                 BANDPASS FILTERING
! =========================================================
!
!
! input:
!    arrayIn(nf,ntrace)= array of FFT | nf=fft_len (number of frequencies)
!	 fmin= min frequency (Hz)
!	 fmax= max frequency (Hz)
!
!    RQ: frequencies should be normalized by fs ...
!	 
! output: 
!	arrayIn bandpass filtered
!
! =========================================================

subroutine bandpass(arrayIn,nf,ntrace,fmin,fmax)

use iso_fortran_env
implicit none

complex(kind=real64),intent(inout)	:: arrayIn(nf,ntrace)
integer,intent(in)                  :: nf,ntrace
real(kind=real64), intent(in)       :: fmin,fmax

integer								:: k 
real(kind=real64)					:: window(nf)




call taper_bandpass(window,nf,fmin,fmax)


!$OMP PARALLEL default(shared), firstprivate(k)
!$OMP DO

do k=1,ntrace
	arrayIn(:,k) = arrayIn(:,k)*window          ! taper (bandpass)	
enddo

!$OMP END PARALLEL


end subroutine bandpass



! =========================================================
!             TAPER  ~ BANDPASS 
! =========================================================
!
!
! input:
!    win(nf)= window | 
!	 nf=number of frequencies in FFT
!	 fmin= min frequency (Hz)
!	 fmax= max frequency (Hz)
!	 fs= frequency sampling (Hz)
!
! output: 
!	win(nf)
!
! =========================================================

subroutine taper_bandpass(win,nf,fmin,fmax) 

use iso_fortran_env
implicit none


real(kind=real64), intent(inout)  :: win(nf)
real(kind=real64), intent(in)     :: fmin,fmax
integer, intent(in)               :: nf

real(kind=real64) :: s,r
real(kind=real64), parameter :: pi = 3.1415926535897931d0 
integer :: j,I1,I2,I3,I4,Napod

! les fréquences sont normalisées 
!I1,I2,I3,I4=fmin,fmin+df,fmax-df,fmax 
!k=nf*(fk/fs): indice k lié à la fréquence fk. si fk est normalisée, k=nf*fk


r=0.2d0                ! cosine fraction, df=(fmax-fmin)*r

Napod=(fmax-fmin)*r*nf ! taper rise width in pts


I1= nf*fmin    		   ! fmin index
I2= I1+Napod	       ! fmin+df index

I4= nf*fmax	           ! fmax index
I3= I4-Napod           ! fmax-df index


win=0.d0
do j=I1,I2
		s=REAL(j-I1)/REAL(Napod) ! 0 <= s <= 1
		s=s*pi/2                 ! 0 <= s <= pi/2
		win(j)=SIN(s)**2
enddo
win(I2+1:I3)=1.d0

do j=I3+1,I4

		s=REAL(j-I3)/REAL(Napod) ! 0 <= s <= 1
		s=(s+1.d0)*pi/2          ! pi/2 <= s <= pi
		win(j)=SIN(s)**2
enddo


end subroutine taper_bandpass


! =========================================================
!                 LOWPASS FILTERING
! =========================================================
!
!
! input:
!    arrayIn(nf,ntrace)= array of FFT | nf=fft_len (number of frequencies)
!	 fmax= max frequency (Hz)
!	 
! output: 
!	arrayIn lowpass filtered
!
! =========================================================

subroutine lowpass(arrayIn,nf,ntrace,fmax)

use iso_fortran_env
implicit none

complex(kind=real64),intent(inout)	:: arrayIn(nf,ntrace)
integer,intent(in)                  :: nf,ntrace
real(kind=real64), intent(in)       :: fmax

integer								:: k 
real(kind=real64)					:: window(nf)




call taper_lowpass(window,nf,fmax)


!$OMP PARALLEL default(shared), firstprivate(k)
!$OMP DO

do k=1,ntrace
	arrayIn(:,k) = arrayIn(:,k)*window          ! taper	(lowpass)
enddo

!$OMP END PARALLEL


end subroutine lowpass


! =========================================================
!             TAPER  ~ BANDPASS 
! =========================================================
!
!
! input:
!    win(nf)= window | 
!	 nf=number of frequencies
!	 fmax= max frequency (Hz)
!
! output: 
!	win(nf)
!
! =========================================================

subroutine taper_lowpass(win,nf,fmax) 

use iso_fortran_env
implicit none


real(kind=real64), intent(inout)  :: win(nf)
real(kind=real64), intent(in)     :: fmax
integer, intent(in)               :: nf

real(kind=real64) :: s,r
real(kind=real64), parameter :: pi = 3.1415926535897931d0 
integer :: j,I1,I2,Napod


!I1,I2=fmax-df,fmax 
!k=nf*(fk/fs): indice k lié à la fréquence fk. si fk est normalisée, k=nf*fk


r=0.2                ! cosine fraction, df=(fmax-0)*r

Napod=(fmax)*r*nf




I2= nf*fmax	         ! fmax index
I1= I2-Napod         ! fmax-df index


win=0.d0

win(1:I1)=1.d0

do j=I1+1,I2

		s=REAL(j-i1)/REAL(Napod)  ! 0 <= s <= 1
		s=(s+1.d0)*pi/2.          ! pi/2 <= s <= pi
		win(j)=SIN(s)**2
enddo


end subroutine taper_lowpass





