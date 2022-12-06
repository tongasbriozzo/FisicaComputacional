program p3b
use precision, only: pr
use, intrinsic :: iso_c_binding
use ODEs
implicit none
include "/usr/include/fftw3.f03"

! Este programa calcula el espectro de potencias de los ángulos del pémduilo doble

  type(C_ptr)                            :: plan_rc
  complex(C_double_complex), allocatable :: out1(:), out2(:)
  real(C_double)                         :: factor

  integer                  :: i, n
  real(pr)                 :: a, b, fr, p1, p2
  real(pr), dimension(1:4) :: x0
  real(pr), allocatable    :: x(:,:), t(:)

  ! Inicializamos los parámetros
  a = 0._pr    ! Tiempo inicial
  b = 1200._pr ! Teimpo final
  n = 800000   ! Número de divisiones del intervalo de integración

  ! Damos dimension a los vectores x (estado del sistema), t (tiempo) y out (transformada de Fourier)
  allocate(x(1:4,0:n))
  allocate(t(0:n))
  allocate(out1(n/2 + 1),out2(n/2 + 1))
  plan_rc = fftw_plan_dft_r2c_1d(n, x(1,:), out1, FFTW_MEASURE)  !  También calculamos el plan de trabajo de la fftw

  ! Sistema 1 con E=0

  ! Inicializamos las condiciones iniciales
  x0(1) = 0._pr           ! tita1
  x0(2) = 0._pr           ! tita2
  x0(3) = sqrt(1.125_pr)  ! omega1
  x0(4) = 0._pr           ! omega2

  ! Hacemos evolucionar el sistema en el tiempo y tomamos la tranformada de fourier de x1 (tita1) y x2 (tita2)
  call RK4(a,b,n,x0,x,t)
  call fftw_execute_dft_r2c(plan_rc, x(1,:), out1)
  call fftw_execute_dft_r2c(plan_rc, x(2,:), out2)
  factor = 1._pr/(t(n)-t(0))

  ! Abrimos el documento de texto y anotamos el espectro de potencias
  open(321,file='ep1.dat')
  do i = 0, n/2
    fr = factor*real(i,kind(1.d0))                                            ! Esta es la frecuencia en la que estamos evaluando la potencia
    p1 = sqrt((real(out1(i+1)/dble(N)))**2 + (aimag(out1(i+1) /dble(N)))**2)  ! Potencia para tita1
    p2 = sqrt((real(out2(i+1)/dble(N)))**2 + (aimag(out2(i+1) /dble(N)))**2)  ! potencia para tita2
    write(321,*) fr, p1, p2
  end do
  close(321)

  ! Sistema 2 con E=-0.745, idem al anterior

  x0(1) = 0._pr
  x0(2) = 0._pr
  x0(3) = sqrt(0.0075_pr)
  x0(4) = 0._pr

  call RK4(a,b,n,x0,x,t)
  call fftw_execute_dft_r2c(plan_rc, x(1,:), out1)
  call fftw_execute_dft_r2c(plan_rc, x(2,:), out2)
  factor = 1._pr/(t(n)-t(0))

  open(322,file='ep2.dat')
  do i = 0, n/2
    fr = factor*real(i,kind(1.d0))
    p1 = sqrt((real(out1(i+1)/dble(N)))**2 + (aimag(out1(i+1) /dble(N)))**2)
    p2 = sqrt((real(out2(i+1)/dble(N)))**2 + (aimag(out2(i+1) /dble(N)))**2)
    write(322,*) fr, p1, p2
  end do
  close(322)

end program p3b
