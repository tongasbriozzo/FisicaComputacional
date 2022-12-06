module ODEssopr
use precision, only: pr
implicit none

contains

function f(x)
implicit none

  real(pr), intent(in) :: x
  real(pr)             :: f, k

  k = 1._pr

  f = -k*x

end function f

function g(y)
implicit none

   real(pr), intent(in) :: y
   real(pr)             :: g

   g = y

end function g

subroutine Eulerso(a,b,n,x0,y0,xen,yen)
implicit none

  integer(4), intent(in)   :: n
  integer(4)               :: i
  real(pr), intent(in)     :: a, b, x0, y0
  real(pr), intent(out)    :: xen, yen
  real(pr)                 :: h, xei, xef, yei, yef

  h = (b-a)/real(n,pr)
  xei = x0
  yei = y0

  do i = 0, n-1

  	yef = yei+h*f(xei)
    xef = xei+h*g(yei)

	  xei = xef
	  yei = yef

  end do

  xen = xef
  yen = yef

end subroutine Eulerso

subroutine RK2so(a,b,n,x0,y0,x2n,y2n)
implicit none

  integer(4), intent(in)   :: n
  integer(4)               :: i
  real(pr), intent(in)     :: a, b, x0, y0
  real(pr), intent(out)    :: x2n, y2n
  real(pr)                 :: K1, K2, L1, L2, h, x2i, x2f, y2i, y2f

  h = (b-a)/real(n,pr)
  x2i = x0
  y2i = y0

  do i = 0, n-1

    K1 = h*g(y2i)
    L1 = h*f(x2i)
    K2 = h*g(y2i+L1)
    L2 = h*f(x2i+K1)

    x2f = x2i+(K1+K2)/2._pr
    y2f = y2i+(L1+L2)/2._pr

	  x2i = x2f
	  y2i = y2f

  end do

  x2n = x2f
  y2n = y2f

end subroutine RK2so

subroutine RK4so(a,b,n,x0,y0,x4n,y4n)
implicit none

  integer(4), intent(in)   :: n
  integer(4)               :: i
  real(pr), intent(in)     :: a, b, x0, y0
  real(pr), intent(out)    :: x4n, y4n
  real(pr)                 :: K1, K2, K3, K4, L1, L2, L3, L4, h, x4i, x4f, y4i, y4f


  h = (b-a)/real(n,pr)
  x4i = x0
  y4i = y0

  do i = 0, n-1

    K1 = h*g(y4i)
    L1 = h*f(x4i)
    K2 = h*g(y4i+L1/2._pr)
    L2 = h*f(x4i+k1/2._pr)
    K3 = h*g(y4i+L2/2._pr)
    L3 = h*f(x4i+k2/2._pr)
    K4 = h*g(y4i+L3)
    L4 = h*f(x4i+k3)

    x4f = x4i+(K1+K2*2._pr+K3*2._pr+K4)/6._pr
    y4f = y4i+(L1+L2*2._pr+L3*2._pr+L4)/6._pr

	  x4i = x4f
	  y4i = y4f

  end do

  x4n = x4f
  y4n = y4f

end subroutine RK4so

end module ODEssopr
