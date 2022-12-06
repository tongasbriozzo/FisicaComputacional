program test
use precision, only: pr
use gasdev
implicit none

  integer                         :: i, j, k
  integer, parameter              :: N=256
  real(pr), dimension(1:N)        :: v1, v2, v3

  call gasdev_v(v1)
  call gasdev_v(v2)
  call gasdev_v(v3)

  open(600,file='tst.d')
  do i = 1, N
    write(600,*) v1(i), v2(i), v3(i)
  end do
  close(600)

end program test
