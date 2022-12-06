program p3a
use precision, only: pr
use ODEs
implicit none

! Este programa calcula las trayectorias del péndulo doble, asi como sus àngulos y velocidades angulares en función del tiempo

  integer                  :: i, n
  real(pr)                 :: a, b, L1, L2, x1, y1, x2, y2
  real(pr), dimension(1:4) :: x0
  real(pr), allocatable    :: x(:,:), t(:)

  ! Primero, definimos los formatos que vamos a utilizar en los documentos de texto
  1       format(A6,    18x, A2,    22x, A2,    22x, A2,    22x, A2)
  2       format(A6,    18x, A5,    19x, A5,    19x, A2,    22x, A2)
  3       format(E22.16, 2x, E22.16, 2x, E22.16, 2x, E22.16, 2x, E22.16)

  ! Inicializamos los parámetros,
  a = 0._pr     !  tiempo inicial
  b = 150._pr   !  tiempo final
  n = 100000    ! numero de diviciones del intervalo de integración

  L1 = 19.6_pr  ! Longitud del péndulo 1
  L2 = 9.8_pr   ! Longitud del péndulo 2

  ! Damos dimensión al vector x, compuesto por tita1, tita2, omega1 y omega2. Tambien dimensionamos el tiempo
  allocate(x(1:4,0:n))
  allocate(t(0:n))

  ! Inicializamos las condiciones iniciales
  x0(1) = 0._pr
  x0(2) = 0._pr
  x0(3) = sqrt(1.125_pr)
  x0(4) = 0._pr

  ! Llamamos a la subrutina rk4 para obtener la evolucion del sistema en el tiempo
  call RK4(a,b,n,x0,x,t)

  ! Abrimos los documento de texto y titulamos las columnas
  open(311,file='pos.dat')
  open(312,file='ang.dat')

  write(311,1)  'tiempo',   'x1',       'y1',       'x2',       'y2'
  write(312,2)  'tiempo',   'tita1',    'tita2',    'w1',       'w2'

  ! Escribimos las posiciones, angulos y velocidad angulares en función del tiempo
  do i = 0, n

    x1 = -L1*cos(x(1,i))       ! coordenada x de la masa 1
    y1 = -L1*sin(x(1,i))       ! coordenada y de la masa 1
    x2 = -L2*cos(x(2,i)) + x1  ! coordenada x de la masa 2
    y2 = -L2*sin(x(2,i)) + y1  ! coordenada y de la masa 2

    write(311,3) t(i),       x1,         y1,         x2,         y2
    write(312,3) t(i),       x(1,i),     x(2,i),     x(3,i),     x(4,i)

  end do

  close(311)
  close(312)

end program p3a
