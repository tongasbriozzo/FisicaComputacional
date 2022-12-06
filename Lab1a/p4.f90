program p4
implicit none
!use pcs, only:pr
!use	gaussmod

  integer            :: i, n
  integer, parameter ::	dp=8
  real(dp)           :: In, a, b, G, h


!10
  1     format(A10,     2x,  A54)
  2     format(I10,     2x,  E22.16)

  In = exp(1._dp)-(1._dp)
  G  = 0._dp
  a  = 0._dp
  b  = 1._dp

!20
  open(41,file='cgl16')
  write(41,1) 'Divisiones', 'Gauss'

  do i = 1, 16, 1

    n = i
	  h = (b-a)/real(n,dp)

	  call gauss(n,a,b,G)

	  write(41,2)	n,           abs((G-In)/In)

  end do

CONTAINS
!40
       subroutine gauss(npts, a, b, IGL)
       implicit none
       integer, intent(in)     :: npts
       real (dp), intent(in)   :: a,b
       real (dp)               :: x(0:npts),w(0:npts)
       real (dp)               :: t,t1,pp,p1,p2,p3,eps
       real (dp)               :: pi
!50
	   real (dp), intent(out)  :: IGL
       integer                 :: m,i,j

	   IGL = 0._dp
       pi = acos(-1._dp)
       m = 0; i = 0; j = 0
       t = 0._dp;  t1 = 0._dp; pp = 0._dp;
       p1 = 0._dp; p2 = 0._dp; p3 = 0._dp
!60
       eps = 1.5e-16_dp
       m = (npts + 1)/2

       do i=1,m + 1

          t = cos(pi*(dble(i) - 0.25_dp)/(dble(npts+1) + 0.5_dp) )
          t1 = 1._dp

          do while( (abs(t - t1) ) >= eps)

             p1 = 1._dp ; p2 = 0._dp
!70
             do j=1,npts + 1
                p3 = p2; p2 = p1
                p1 = ( (2._dp*dble(j) - 1._dp)*t*p2 - dble(j-1)*p3)/(dble(j) )
             enddo

             pp = dble(npts+1)*(t*p1 - p2)/(t*t - 1._dp)
             t1 = t
!80
             t = t1 - p1/pp

          enddo

 !         x(i - 1) = - t;
 !         x(npts + 1 - i) = t
 !         w(i - 1) = 2._dp/( (1._dp - t*t)*pp*pp)
 !         w(npts + 1 - i) = w(i - 1)
!90
          x(i - 1) = - t*(b - a)/2._dp + (b + a)/2._dp;
          x(npts + 1 - i) = t*(b - a)/2._dp + (b + a)/2._dp
          w(i - 1) = 2._dp/( (1._dp - t*t)*pp*pp)*(b - a)/2._dp
          w(npts + 1 - i) = w(i - 1)

       enddo

          do i=0, npts
!            x(i) = x(i)*(b - a)/2._dp + (b + a)/2._dp
!            w(i) = w(i)*(b - a)/2._dp
			 IGL = IGL + exp(x(i))*w(i)
          end do

       return
       end subroutine gauss

end program p4
