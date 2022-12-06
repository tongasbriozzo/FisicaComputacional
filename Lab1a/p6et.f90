program p6
use precision, only: pr
use ODEsso
implicit none

  1     format(A6,    18x, A8,    16x, A9)
  2     format(E22.16, 2x, E22.16, 2x, E22.16)

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

  open(61, file='eul.dat')
	open(62, file='rk2.dat')
  open(63, file='rk4.dat')
  write(61,1) 'Tiempo',   'Posición', 'Velocidad'
  write(62,1) 'Tiempo',   'Posición', 'Velocidad'
  write(63,1) 'Tiempo',   'Posición', 'Velocidad'

  i = 25
  n = 2**i
	h = (b-a)/real(n,pr)

	allocate(t(0:n),xe(0:n),ye(0:n))

	do j = 0, n
    t(j) = a+h*real(j,pr)
  end do

	call Eulerso(a,b,n,x0,y0,xe,ye)

  do k = 0, n, 32768
	  write(61,2) t(k), abs(xe(k)-x(t(k))), abs(ye(k)-v(t(k)))
	end do

  deallocate(t,xe,ye)

  i = 25
  n = 2**i
	h = (b-a)/real(n,pr)

	allocate(t(0:n),x2(0:n),y2(0:n))

	do j = 0, n
    t(j) = a+h*real(j,pr)
  end do

	call   RK2so(a,b,n,x0,y0,x2,y2)

  do k = 0, n, 32768
	  write(62,2) t(k), abs(x2(k)-x(t(k))), abs(y2(k)-v(t(k)))
	end do

  deallocate(t,x2,y2)

  i = 14
  n = 2**i
	h = (b-a)/real(n,pr)

	allocate(t(0:n),x4(0:n),y4(0:n))

	do j = 0, n
    t(j) = a+h*real(j,pr)
  end do

	call   RK4so(a,b,n,x0,y0,x4,y4)

  do k = 0, n, 16
	  write(63,2) t(k), abs(x4(k)-x(t(k))), abs(y4(k)-v(t(k)))
	end do

  deallocate(t,x4,y4)

  close(61)
	close(62)
  close(63)

contains

  function x(t)
  implicit none
    real(pr), intent(in) :: t
    real(pr)             :: x, pi

    pi = acos(-1._pr)
    x  = sqrt(2._pr)*sin(t+pi/4._pr)

  end function x

  function v(t)
  implicit none
    real(pr), intent(in) :: t
    real(pr)             :: v, pi

    pi = acos(-1._pr)
    v  = sqrt(2._pr)*cos(t+pi/4._pr)

  end function v

end program p6
