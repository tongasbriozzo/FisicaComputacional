program test
use precision, only: pr
implicit none

  real(pr), parameter :: a=12.55592670560186_pr, b=-21.29665330186477, c=12.04072659926291
  integer :: i
  real(pr) :: T, x

  open(100,file='test.d')
  do i = 1, 500
    x=real(i,pr)/500._pr
    T=a*x**3 + b*x**2 + c*x
    write(100,*) i, T
  end do
  close(100)

end program
