program p6pr
use precision, only: pr
use ODEssopr
implicit none

  1      format(A10,     2x, A11,            11x, A9,           13x, A9)
  2      format(I10,     2x, E22.16,          2x, E22.16,        2x, E22.16)

  integer(kind=4)       :: i, n
  real(pr)              :: a, b, h, x0, y0
  real(pr)              :: xe, ye, x2, y2, x4, y4
  real(pr)              :: pi, xf

  pi = acos(-1._pr)
  xf = sqrt(2._pr)*sin(10._pr+pi/4._pr)

  a = 0._pr
  b = 10._pr

  x0 = 1._pr
  y0 = 1._pr

  open(61, file = 'IntegralesError')
  write(61,1)  'Nº de ptos', 'Error Euler',      'Error RK2',        'Error RK4'

  do i = 1, 30, 1

    n = (2_4)**i
	  h = (b-a)/real(n,pr)

	  call Eulerso(a,b,n,x0,y0,xe,ye)
	  call   RK2so(a,b,n,x0,y0,x2,y2)
	  call   RK4so(a,b,n,x0,y0,x4,y4)

	  write(61,2)	n,     abs((xe-xf)/xf), abs((x2-xf)/xf), abs((x4-xf)/xf)

  end do

  close(61)

end program p6pr
