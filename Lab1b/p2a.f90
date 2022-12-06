program p2a

! Este porograma calcula las trayectorias de la ec logística para distintas condiciones iniciales, a un mismo r

use precision, only: pr
implicit none
  integer                                :: i, n
  real(pr), dimension(1:500)             :: x1, x2, x3, x4
  real(pr)                               :: r

  ! Primero, definimos los formatos que usaremos en los documentos de texto
  1      format(A6,    18x, A7,    17x, A6,    18x, A6,    18x, A6)
  2      format(I3,     2x, E22.16, 2x, E22.16, 2x, E22.16, 2x, E22.16)

  ! elegimos n y r, donde n es el nùmero de pasos y r el paràmetro variable
  n  = 500
  r  = 4._pr
  ! Para evitar trabajar con vectores y manejar muchos datos a la vez, se prefirio cambiar de r manualmente

  ! Elejimos los valores iniciales con los que trabajaremos
  x1(1) = 0.06_pr
  x2(1) = 0.3_pr
  x3(1) = 0.6_pr
  x4(1) = 0.9_pr

  ! Abrimos el documento de texto
  open(21,file='xt.dat')
  write(21,1)  'Tiempo',   'x0=0.06',  'x0=0.3',   'x0=0.6',   'x0=0.9'
  write(21,1)   1,          x1(1),      x2(1),      x3(1),      x4(1)

  ! Hacemos evolucionar el sistema hasta el instante n
  do i = 1, n-1

    x1(i+1) = r*x1(i)*(1._pr-x1(i))
    x2(i+1) = r*x2(i)*(1._pr-x2(i))
    x3(i+1) = r*x3(i)*(1._pr-x3(i))
    x4(i+1) = r*x4(i)*(1._pr-x4(i))

    ! Escribimos las trayectorias obtenidas
    write(21,2) i+1,        x1(i+1),    x2(i+1),    x3(i+1),    x4(i+1)

  end do

  ! Cerramos el archibo de texto
  close(21)

end program p2a
