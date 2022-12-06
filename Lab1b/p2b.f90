program p2b

! Este programa calcula el espectro de potencias para las trayectorias obtenidas en el problema anterior
! para esto, hacemos uso de la libreria fftw3 y su subrutina

use precision, only: pr
use, intrinsic :: iso_c_binding
implicit none
include "/usr/include/fftw3.f03"
  type(C_ptr)                            :: plan_rc
  complex(C_double_complex), allocatable :: out(:)
  integer                                :: i, n
  real(C_double)                         :: factor
  real(pr), dimension(1:500)             :: x1
  real(pr)                               :: r

  ! Elejimos el número de pasos y el valor de r
  n  = 500
  r  = 1.5_pr
  ! para evitar trabajar con vectores y manejar demasiados datos a la vez, r debe elejirse manualmente

  ! Damos dimensión al vector out, que contendrá los valores del espectro de potencias
  ! También calculamos el plan de acción para la subrutina fftw3
  allocate( out(n/2 + 1))
  plan_rc = fftw_plan_dft_r2c_1d(n, x1, out, FFTW_MEASURE)

  ! Inicializamos el punto de partida de la trayectoria y la hacemos evolucionar en el tiempo
  x1(1) = 0.6_pr
  do i = 1, n-1
    x1(i+1) = r*x1(i)*(1._pr-x1(i))
  end do
  ! Descartamos un transitorio de 500 pasos y volvemos a calcular la trayectoria para los 500 pasos siguientes
  x1(1) = x1(500)
  do i = 1, n-1
    x1(i+1) = r*x1(i)*(1._pr-x1(i))
  end do

  ! Llamamo a la subrutina fftw para obtener el espectro de potencia de la trayectoria
  call fftw_execute_dft_r2c(plan_rc, x1, out)
  ! Definimos la frecuencia mìnima del espectro
  factor = 1.d0/500._pr

  ! abrimos el documento de texto y anotamos el modulo de la potencia para cada frecuencia
  open(23,file='p1.5.dat')
  do i = 0, 300  ! se anota hasta 300 y no hasta n/2 para evitar problemas de borde en el gráfico, no cambia el resultado
    write(23,*) factor*real(i,kind(1.d0)), sqrt((real(out(i+1) /dble(N)))**2 + (aimag(out(i+1) /dble(N)))**2)
  end do
  close(23)

end program p2b
