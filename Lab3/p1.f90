program p1
use precision, only: pr
use random
use mzranmod
use mtmod
implicit none

  integer    :: i, j, k
  real(pr)   :: x, u, x1, x2, x3, x4

  do i = 1, 10000



  end do

contains

  function f(x)
  implicit none
    real(pr), intent(in) :: x
    real(pr)             :: f
    f = 1._pr/(x*x)
  end function f

  function g(x,sigma,mu)
  implicit none
    real(pr), intent(in) :: x, sigma, mu
    real(pr)             :: g
    real(pr), parameter  :: pi=acos(-1._pr)
    g = (exp(-0.5_pr*(x-mu)*(x-mu)/(sigma*sigma)))/sqrt(2._pr*pi*sigma*sigma)
  end function g

  function e(x,lambda)
  implicit none
    real(pr), intent(in) :: x, lambda
    real(pr)             :: e
    e = (exp(-x/lambda))/lambda
  end function e

end program p1
