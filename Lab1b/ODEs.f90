module ODEs
use precision, only: pr
implicit none

! Este módulo contiene la subrutina rk4, empleada para obtener la evolucion de ODEs en el tiempo

contains

function f(t,x)
implicit none

! Esta función contiene las expresiones de la derivada de cada una de las cuatro variables que nos interesan

  real(pr), intent(in)                 :: t
  real(pr), intent(in), dimension(1:4) :: x
  real(pr), dimension(1:4)             :: f
  real(pr)                             :: a, b, g
  real(pr)                             :: s, s1, s2, c

  ! Inicializamos los parámetros a=m2/m1, b=L2/L1 y g=g/L1
  a = 1._pr/3._pr
  b = 0.5_pr
  g = 0.5_pr

  ! Definimos algunas variables que nos serán útiles
  s  = sin(x(1)-x(2))
  s1 = sin(x(1))
  s2 = sin(x(2))
  c  = cos(x(1)-x(2))

  ! Definimos la función f
  f(1) = x(3)                                                                                  ! La derivada de x1 (tita1) es x3 (omega1)
  f(2) = x(4)                                                                                  ! La derivada de x2 (tita2) es x4 (omega2)
  f(3) = (-(1._pr+a)*g*s1-a*b*(x(4)**2)*s-a*c*((x(3)**2)*s-g*s2))/(1._pr+a*s**2)               ! La derivada de x3 (omega1)
  f(4) = ((1._pr+a)*((x(3)**2)*s-g*s2)+c*((1._pr+a)*g*s1+a*b*(x(4)**2)*s))/(b*(1._pr+a*s**2))  ! la derivada de x4 (omega2)

end function f

subroutine RK4(a,b,n,x0,x,t)
implicit none

! Esta subrutina calcula la evolucion temporal de un sistema de ODEs

  integer, intent(in)                   :: n
  integer                               :: i
  real(pr), dimension(0:n)              :: t
  real(pr), intent(in)                  :: a, b
  real(pr), intent(in), dimension (1:4) :: x0
  real(pr), dimension(1:4,0:n)          :: x
  real(pr), dimension(1:4)              :: K1, K2, K3, K4
  real(pr)                              :: h

  ! Inicializamos las condiciones iniciales y el step size
  x(:,0) = x0
  h = (b-a)/real(n,pr)

  ! Definimos los tiempos a evaluar
  do i = 0, n
    t(i) = a+h*real(i,pr)
  end do

  ! Obtenemos el vextor x para cada tiempo
  do i = 0, n-1

    K1 = h*f(t(i),x(:,i))
    K2 = h*f(t(i)+h/2._pr,x(:,i)+K1(:)/2._pr)
    K3 = h*f(t(i)+h/2._pr,x(:,i)+K2(:)/2._pr)
    K4 = h*f(t(i)+h,x(:,i)+K3(:))

    x(:,i+1) = x(:,i)+(K1(:)+K2(:)*2._pr+K3(:)*2._pr+K4(:))/6._pr

  end do

end subroutine RK4

end module ODEs
