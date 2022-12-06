program tc
use precision, only: pr
implicit none

  integer                  :: i, j, k
  real(pr), dimension(1:9) :: tce, tcm
  real(pr), dimension(1:9) :: te,  tm

  tce = 0._pr
  tcm = 0._pr

  open(451,file='ac1.d',status='old')
  open(452,file='ac2.d',status='old')
  open(453,file='ac3.d',status='old')
  open(454,file='ac4.d',status='old')
  open(455,file='ac5.d',status='old')
  open(456,file='ac6.d',status='old')
  open(457,file='ac7.d',status='old')
  open(458,file='ac8.d',status='old')
  open(459,file='ac9.d',status='old')

  open(461,file='te.d')
  open(462,file='tm.d')

  do i=1, 1000
    read(451,*) k, te(1), tm(1)
    read(452,*) k, te(2), tm(2)
    read(453,*) k, te(3), tm(3)
    read(454,*) k, te(4), tm(4)
    read(455,*) k, te(5), tm(5)
    read(456,*) k, te(6), tm(6)
    read(457,*) k, te(7), tm(7)
    read(458,*) k, te(8), tm(8)
    read(459,*) k, te(9), tm(9)
    tce = tce + te
    tcm = tcm + tm
    write(461,*) i, tce(1), tce(2), tce(3), tce(4), tce(5), tce(6), tce(7), tce(8), tce(9)
    write(462,*) i, tcm(1), tcm(2), tcm(3), tcm(4), tcm(5), tcm(6), tcm(7), tcm(8), tcm(9)
  end do

end program tc
