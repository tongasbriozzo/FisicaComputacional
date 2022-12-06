module ODEsso
use precision, only: pr
implicit none

contains

function f(t,x)
implicit none

  real(pr), intent(in) :: t, x
  real(pr)             :: f, k

  k = 1._pr

  f = -k*x

end function f

function g(t,y)
implicit none

   real(pr), intent(in) :: t, y
   real(pr)             :: g

   g = y

end function g

subroutine Eulerso(a,b,n,x0,y0,x,y)
implicit none

  integer, intent(in)      :: n
  integer                  :: i
  real(pr), dimension(0:n) :: t
  real(pr), intent(in)     :: a, b, x0, y0
  real(pr), dimension(0:n) :: x, y
  real(pr)                 :: h

  x(0) = x0
  y(0) = y0
  h = (b-a)/real(n,pr)

  do i = 0, n

    t(i) = a+h*real(i,pr)

  end do

  do i = 0, n-1

    x(i+1) = x(i)+h*g(t(i),y(i))
	  y(i+1) = y(i)+h*f(t(i),x(i))

  end do

end subroutine Eulerso

subroutine RK2so(a,b,n,x0,y0,x,y)
implicit none

integer, intent(in)      :: n
integer                  :: i
real(pr), dimension(0:n) :: t
real(pr), intent(in)     :: a, b, x0, y0
real(pr), dimension(0:n) :: x, y
real(pr)                 :: K1, K2, L1, L2, h

x(0) = x0
y(0) = y0
h = (b-a)/real(n,pr)

do i = 0, n

  t(i) = a+h*real(i,pr)

end do

do i = 0, n-1

  K1 = h*g(t(i),y(i))
  L1 = h*f(t(i),x(i))
  K2 = h*g(t(i)+h,y(i)+L1)
  L2 = h*f(t(i)+h,x(i)+K1)

  x(i+1) = x(i)+(K1+K2)/2._pr
  y(i+1) = y(i)+(L1+L2)/2._pr

  end do

end subroutine RK2so

subroutine RK4so(a,b,n,x0,y0,x,y)
implicit none

integer, intent(in)      :: n
integer                  :: i
real(pr), dimension(0:n) :: t
real(pr), intent(in)     :: a, b, x0, y0
real(pr), dimension(0:n) :: x, y
real(pr)                 :: K1, K2, K3, K4, L1, L2, L3, L4, h

x(0) = x0
y(0) = y0
h = (b-a)/real(n,pr)

do i = 0, n

  t(i) = a+h*real(i,pr)

end do

do i = 0, n-1

  K1 = h*g(t(i),y(i))
  L1 = h*f(t(i),x(i))
  K2 = h*g(t(i)+h/2._pr,y(i)+L1/2._pr)
  L2 = h*f(t(i)+h/2._pr,x(i)+K1/2._pr)
  K3 = h*g(t(i)+h/2._pr,y(i)+L2/2._pr)
  L3 = h*f(t(i)+h/2._pr,x(i)+K2/2._pr)
  K4 = h*g(t(i)+h,y(i)+L3)
  L4 = h*f(t(i)+h,x(i)+K3)

  x(i+1) = x(i)+(K1+K2*2._pr+K3*2._pr+K4)/6._pr
  y(i+1) = y(i)+(L1+L2*2._pr+L3*2._pr+L4)/6._pr

  end do

end subroutine RK4so

end module ODEsso
