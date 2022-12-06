program  p1
use, intrinsic :: iso_c_binding
implicit none
include "/usr/include/fftw3.f03"
  type(C_ptr)                            :: plan_rc,plan_cr
  real(C_double), allocatable            :: in(:)
  complex(C_double_complex), allocatable :: out(:)
  integer                                :: i, j, N
  real(C_double)                         :: t, ti, tf, factor, tt0, tt1, tt2, tt3

  open(111,file='fft.dat')

  N = 1024

  allocate( in(N), out(N/2 + 1))

  call cpu_time(tt0)
  plan_rc = fftw_plan_dft_r2c_1d(n, in, out, FFTW_MEASURE) !FFTW_ESTIMATE
  call cpu_time(tt1)
  write(*,*) 'tiempo plan', (tt1-tt0)

  ti = 0.d0
  tf = 4.d0

  do i=1,N
    t = dble(i-1)*(tf-ti)/dble(N)+ti
    in(i) = f(t)
  end do

  call cpu_time(tt2)
  call fftw_execute_dft_r2c(plan_rc, in, out)
  call cpu_time(tt3)
  write(*,*) 'tiempo transformando', tt3-tt2

  factor = 1.d0/(tf-ti)
  do i = -N/2, N/2
    if ( i < 0 ) then
      write(111,*) factor*real(i,kind(1.d0)),real(conjg( out(-i+1) )/dble(N)), aimag(conjg( out(-i+1) )/dble(N))
    else
      write(111,*) factor*real(i,kind(1.d0)), real(out(i+1) /dble(N)), aimag(out(i+1) /dble(N))
    endif
  enddo

  close(111)

CONTAINS

  function f(x)
  implicit none
    real(kind(1.d0)), intent(in) :: x
    real(kind(1.d0))             :: f
    real(kind(1.d0)), parameter  :: pi=4.d0*atan(-1.d0)

    f = sin(pi/2.d0*x) + cos (20.d0*pi*x)

  end function f

end program p1
