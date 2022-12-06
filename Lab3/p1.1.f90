program p1
use precision, only: pr
use random
use mzranmod
use mtmod
implicit none

  integer                    :: i, j, k, im, im2
  real(pr)                   :: x, x1, x2, x3, x4
  real(pr)                   :: u, u1, u2, u3, u4
  real(pr)                   :: m, d
  real(pr), dimension(1:100) :: y

  im  = 2147483647
  im2 = int(real(im)/2._pr)

  u1 = ran0(im2)
  u2 = ran2(im2)

  do k = 1, 10000

    u1 = ran0(int(u1)*im)
    u2 = ran2(int(u2)*im)
    u3 = rmzran()
    u4 = grnd()

    do i = 1, 100
      y(i) = 0._pr
    end do

    ! Asignamos los valores a y(i)
    do i = 1, 10000
      m = x(i)*100._pr
      do j = 1, 100
        d = m-real(j,pr)
        if (d<=0) then
          y(j) = y(j) + 1._pr
          exit
        end if
      end do
    end do

  end do

contains

  function f(u)
  implicit none
    real(pr), intent(in) :: u
    real(pr)             :: f
    f = -1._pr/u
  end function f

  function g1(u1,u2)
  implicit none
    real(pr), intent(in) :: u1, u2
    real(pr)             :: g1
    real(pr), parameter  :: pi=acos(-1._pr)
    g1 = sqrt(-2._pr*log(u1))*cos(2._pr*pi*u2)
  end function g1
  function g2(u1,u2)
  implicit none
    real(pr), intent(in) :: u1, u2
    real(pr)             :: g2
    real(pr), parameter  :: pi=acos(-1._pr)
    g2 = sqrt(-2._pr*log(u1))*sin(2._pr*pi*u2)
  end function g2

  function e(u,lambda)
  implicit none
    real(pr), intent(in) :: u, lambda
    real(pr)             :: e
    e = -log(u)/lambda
  end function e

end program p1
