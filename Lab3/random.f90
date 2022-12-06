module randommod
use precision, only: pr
implicit none

!  INTEGER,PARAMETER :: K4B=selected_int_kind(9)
!  INTEGER,PARAMETER :: DP =KIND(1.0D0)

CONTAINS

  function ran0(idum)
    implicit none
    integer             :: idum
    real(pr)            :: ran0
    integer, parameter  :: IA=16807, IM=2147483647, IQ=127773, IR=2836, MASK=123459876
    real(pr), parameter :: AM=1._pr/IM
    integer             :: k

    idum = ieor(idum,MASK)
    k    = idum/IQ
    idum = IA*(idum-k*IQ) - IR*k

    if (idum.lt.0) then
      idum = idum + IM
    end if

    ran0 = AM*idum
    idum = ieor(idum,MASK)

    return

  end function ran0

  function ran2(idum)
    implicit none
    integer             :: idum
    real(pr)            :: ran2
    integer, parameter  :: IM1=2147483563, IM2=2147483399, IMM1=IM1-1
    integer, parameter  :: IA1=40014, IA2=40692, IQ1=53668, IQ2=52774
    integer, parameter  :: IR1=12211, IR2=3791, NTAB=32, NDIV=1+67108861
    real(pr), parameter :: AM=1._pr/IM1, EPS=1.2_pr/10000000._pr, RNMX=1._pr-EPS
    integer             :: idum2,j,k,iv(NTAB),iy

    SAVE iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/

    if (idum.le.0) then
      idum  = max(-idum,1)
      idum2 = idum
      do j = NTAB+8, 1, -1
        k = idum/IQ1
        idum = IA1*(idum-k*IQ1) - k*IR1
        if (idum.lt.0) idum=idum+IM1
        if (j.le.NTAB) iv(j)=idum
      end do
      iy=iv(1)
    endif

    k     = idum/IQ1
    idum  = IA1*(idum-k*IQ1) - k*IR1
    k     = idum2/IQ2
    idum2 = IA2*(idum2-k*IQ2) - k*IR2

    if (idum2.lt.0) idum2 = idum2 + IM2
      j  = 1 + iy/NDIV
      iy = iv(j) - idum2
    if(iy.lt.1) iy = iy + IMM1
      ran2 = min(AM*iy,RNMX)
    return

  end function ran2

!  function rmzran()
!        implicit none
!        integer  ::  mzran,n,i,j,k,is,js,ks,ns,mzranset
!        real(pr) ::  rmzran
!
!        save i,j,k,n
!        data i,j,k,n/521288629,362436069,16163801,1131199299/
!        mzran = i-k
!        if(mzran.lt.0) mzran = mzran + 2147483579
!        i = j
!        j = k
!        k = mzran
!        n = 69069*n + 1013904243
!!        n = ishft(3533*ishft(n,-16)+iand(n,65535),16)+3533*iand(n,65535))
!!       n = n + 1013904243
!        mzran = mzran + n
!!       For random reals on  (0,1): UNI() = .5+.2328306E-9*mzran()
!!       For random reals on (-1,1): VNI() = .4656613E-9*mzran()
!
!        rmzran = real(mzran,kind=DP)*0.2328306*10**(-9) + 0.5_pr
!
!        return
!        entry mzranset(is,js,ks,ns)
!        i = 1+iabs(is)
!        j = 1+iabs(js)
!        k = 1+iabs(ks)
!        n = ns
!        mzranset = n
!        return
!  end function rmzran

end module randommod
