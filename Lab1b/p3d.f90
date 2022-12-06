program p3d
use precision, only: pr
implicit none

! Este programa calcula el tiempo de flip (tf) para distintas condiciones iniciales de tita1 y tita2, partiendo con w1=w2=0

  integer                          :: i, j, k, n
  real(pr)                         :: a, b, t, h, det
  real(pr), dimension(1:4)         :: x0, xi, xf
  real(pr), parameter              :: pi = acos(-1._pr)

  ! Inicializamos los parámetros
  a = 0._pr      ! Tiempo inicial
  b = 10000._pr  ! Tiempo final
  n = 100000000  ! Número de divisiones del intervalo de integración. Si, llevó como 20 horas

  h = (b-a)/real(n,pr)  ! step size

  x0(:) = 0._pr  ! x0 nos da la condicion inicial que utilizaremos para calcular el tf
  xi(:) = 0._pr  ! xi es el estado inicial para calcular el siguiente paso por rk4
  xf(:) = 0._pr  ! xf es el siguiente estado, calculamos la evolución paso a paso para no guardar vectores que ocupen demasiada memoria

  ! Abrimos el documento de texto
  open(33,file='tf.dat')

  !  Tomamos los angulos iniciales de una grilla de 601x601 valores equiespaciados en el rango de -3 a 3 radianes (notar que tf(i,j)=tf(-i,-j))
  do i =    0, 300
  do j = -300, 300

    ! Inicializamos el estado inicial del sistema
    x0(1) = real(i,pr)/100._pr  ! tita1 inicial
    x0(2) = real(j,pr)/100._pr  ! tita2 inicial

    det = 2._pr*cos(x0(1))+cos(x0(2))      ! El determinante nos indica si el sistema tiene energía suficiente para realizar un flip

    if (det>1._pr) then                    ! Si el determinante es mayor a 1, no hay flip y seteamos tf en 100 000
      write(33,*)  x0(1),  x0(2), 100000
      write(33,*) -x0(1), -x0(2), 100000
      cycle
    end if

    ! Si det<1, inicializamos el estado inicial del sistema xi
    xi(:) = x0(:)

    ! Hacemos que el sistema evolucione en el tiempo
    do k = 0, n

      t = a + h*real(k,pr)  ! Tiempo

      call RK4(h,t,xi,xf)   ! Llamamos a rk4 para alcular el estado del sistema en el intante siguiente

      if (abs(xf(1))>pi) then           ! Si x1 (tita1) es mayor a pi, escribimos el tf
        write(33,*)  x0(1),  x0(2), t
        write(33,*) -x0(1), -x0(2), t
        exit                            ! Salimos para no desperdiciar tiempo
      end if
      if (abs(xf(2))>pi) then           ! Si x2 (tita2) es mayor a pi, escribimos el tf
        write(33,*)  x0(1),  x0(2), t
        write(33,*) -x0(1), -x0(2), t
        exit
      end if
      if (k==n) then                     ! Si agotamos el tiempo, definimos el tf como 100 000
        write(33,*)  x0(1),  x0(2), 100000
        write(33,*) -x0(1), -x0(2), 100000
      end if

      xi(:) = xf(:)                      ! Definimos el estado inicial del siguiente paso como el final de este

    end do

  end do
  end do

  close(33)

contains

! Estas funciones y subutinas se encuentran explicadas en el archivo ODEs.f90
! Únicamente se altero la subrutina rk4 para que calcule la evolución paso a paso, en vez de para todo el intervalo

function f(t,x)
implicit none

  real(pr), intent(in)                 :: t
  real(pr), intent(in), dimension(1:4) :: x
  real(pr), dimension(1:4)             :: f
  real(pr)                             :: a, b, g
  real(pr)                             :: s, s1, s2, c

  a = 1._pr
  b = 1._pr
  g = 1._pr

  s  = sin(x(1)-x(2))
  s1 = sin(x(1))
  s2 = sin(x(2))
  c  = cos (x(1)-x(2))

  f(1) = x(3)
  f(2) = x(4)
  f(3) = (-(1._pr+a)*g*s1-a*b*(x(4)**2)*s-a*c*((x(3)**2)*s-g*s2))/(1._pr+a*s**2)
  f(4) = ((1._pr+a)*((x(3)**2)*s-g*s2)+c*((1._pr+a)*g*s1+a*b*(x(4)**2)*s))/(b*(1._pr+a*s**2))

end function f

subroutine RK4(h,t,xi,xf)
implicit none

  real(pr), intent(in)                  :: h, t
  real(pr), intent(in), dimension (1:4) :: xi
  real(pr), dimension(1:4)              :: xf
  real(pr), dimension(1:4)              :: K1, K2, K3, K4

    K1(:) = h*f(t,        xi(:))
    K2(:) = h*f(t+h/2._pr,xi(:)+K1(:)/2._pr)
    K3(:) = h*f(t+h/2._pr,xi(:)+K2(:)/2._pr)
    K4(:) = h*f(t+h,      xi(:)+K3(:))

    xf(:) = xi(:)+(K1(:)+K2(:)*2._pr+K3(:)*2._pr+K4(:))/6._pr

end subroutine RK4

end program p3d
