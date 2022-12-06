module RUKU
use pcs
implicit none

contains

function f(x,y)
use pcs
implicit none

  real(dp), intent(in)     :: x, y(1:2)
  real(dp), dimension(1:2) :: f

  f(1) = -10._dp*sin(y(1))
  f(2) = y(1)
  
  
end function f

subroutine RK4(a,b,n,y0,y)
use pcs
implicit none

  integer, intent(in)                  :: n
  integer                              :: i
  real(dp), dimension(0:n)             :: x
  real(dp), intent(in)                 :: a, b
  real(dp), dimension(1:2), intent(in) :: y0
  real(dp), dimension(1:2,0:n)         :: y
  real(dp), dimension(1:2)             :: K1, K2, K3, K4
  real(dp)                             :: h
  
  
  y(:,0) = y0
  h = (b-a)/real(n,dp)
  
  do i = 0, n
  
    x(i) = a+h*real(i,dp)
  
  end do
  
  do i = 0, n-1
  
    K1 = f(x(i),y(i))
    K2 = f(x(i)+h/2._dp,y(i)+k1*h/2._dp)
    K3 = f(x(i)+h/2._dp,y(i)+k2*h/2._dp)
    K4 = f(x(i)+h,y(i)+k3*h)
    
    y(i+1) = y(i)+h*(K1+K2+K3+K4)/6._dp
  
  end do
  
end subroutine RK4

end module RUKU