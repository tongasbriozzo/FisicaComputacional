program energia
use precision, only: pr
use ODEsso
implicit none

  1     format(A6,    18x, A5,    19x, A3,    21x, A3)
  2     format(E22.16, 2x, E22.16, 2x, E22.16, 2x, E22.16)

  integer               :: i, j, k, n
  real(pr)              :: a, b, h, x0, y0
  real(pr), allocatable :: t(:), xe(:), ye(:), x2(:), y2(:), x4(:), y4(:)
  real(pr)              :: pi, xf

  pi = acos(-1._pr)
  xf = sqrt(2._pr)*sin(10._pr+pi/4._pr)

  a = 0._pr
  b = 10._pr

  x0 = 1._pr
  y0 = 1._pr

  open(61, file='ec.dat')
  open(62, file='ep.dat')
  open(63, file='et.dat')
  write(61,1) 'Tiempo',   'Euler',    'RK2',      'RK4'
  write(62,1) 'Tiempo',   'Euler',    'RK2',      'RK4'
  write(63,1) 'Tiempo',   'Euler',    'RK2',      'RK4'

  i = 14

  n = 2**i
	h = (b-a)/real(n,pr)

	allocate(t(0:n),xe(0:n),ye(0:n),x2(0:n),y2(0:n),x4(0:n),y4(0:n))

	do j = 0, n
    t(j) = a+h*real(j,pr)
  end do

	call Eulerso(a,b,n,x0,y0,xe,ye)
	call   RK2so(a,b,n,x0,y0,x2,y2)
	call   RK4so(a,b,n,x0,y0,x4,y4)

	do k = 0, n, 1
	  write(61,2) t(k),      0.5_pr*xe(k)**2,      0.5_pr*x2(k)**2,      0.5_pr*x4(k)**2
    write(62,2) t(k),      0.5_pr*ye(k)**2,      0.5_pr*y2(k)**2,      0.5_pr*y4(k)**2
    write(63,2) t(k), abs(0.5_pr*(xe(k)**2+ye(k)**2)-1), abs(0.5_pr*(x2(k)**2+y2(k)**2)-1), abs(0.5_pr*(x4(k)**2+y4(k)**2)-1)
	end do

  close(61)
  close(62)
  close(63)

end program energia
