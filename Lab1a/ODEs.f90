module ODEs
implicit none

integer, parameter :: sp=4, dp=8

contains

subroutine Euler(a,b,n,y0,y)
implicit none

  integer, intent(in)           :: n
  integer                       :: i
  real(kind=dp), dimension(0:n) :: x
  real(kind=dp), intent(in)     :: a, b, y0
  real(kind=dp), dimension(0:n) :: y
  real(kind=dp)                 :: h
  
  
  y(0) = y0
  h = (b-a)/real(n,kind=dp)
  
  do i = 0, n
  
    x(i) = a+h*real(i,kind=dp)
  
  end do
  
  do i = 0, n-1
  
    y(i+1) = y(i)+h*f(x(i),y(i))
  
  end do

end subroutine Euler

subroutine RK2(a,b,n,y0,y)
implicit none

  integer, intent(in)           :: n
  integer                       :: i
  real(kind=dp), dimension(0:n) :: x
  real(kind=dp), intent(in)     :: a, b, y0
  real(kind=dp), dimension(0:n) :: y
  real(kind=dp)                 :: K1, K2, h
  	  
  y(0) = y0
  h = (b-a)/real(n,kind=dp)
  
  do i = 0, n
  
    x(i) = a+h*real(i,kind=dp)
  
  end do
  
  do i = 0, n-1
  
    K1 = h*f(x(i),y(i))
    K2 = h*f(x(i)+h,y(i)+K1)
    
    y(i+1) = y(i)+(K1+K2)/2._dp
  
  end do

end subroutine RK2

subroutine RK4(a,b,n,y0,y)
implicit none

  integer, intent(in)           :: n
  integer                       :: i
  real(kind=dp), dimension(0:n) :: x
  real(kind=dp), intent(in)     :: a, b, y0
  real(kind=dp), dimension(0:n) :: y
  real(kind=dp)                 :: K1, K2, K3, K4, h
  
  
  y(0) = y0
  h = (b-a)/real(n,kind=dp)
  
  do i = 0, n
  
    x(i) = a+h*real(i,kind=dp)
  
  end do
  
  do i = 0, n-1
  
    K1 = f(x(i),y(i))
    K2 = f(x(i)+h/2._dp,y(i)+k1*h/2._dp)
    K3 = f(x(i)+h/2._dp,y(i)+k2*h/2._dp)
    K4 = f(x(i)+h,y(i)+k3*h)
    
    y(i+1) = y(i)+h*(K1+K2+K3+K4)/6._dp
  
  end do
  
end subroutine RK4

end module ODEs
