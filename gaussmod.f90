       module gaussmod
       use precision, only: pr
       CONTAINS
!
       subroutine gauss(npts, a, b, x, w)
       implicit none
       integer, intent(in)     :: npts
       real (pr), intent(in)   :: a,b
       real (pr), intent(out)  :: x(0:npts),w(0:npts)
       real (pr)               :: t,t1,pp,p1,p2,p3,eps
       real (pr)               :: pi
       real (pr)               :: xi
       integer                 :: m,i,j

       pi = acos(-1._pr)
       m = 0; i = 0; j = 0
       t = 0._pr;  t1 = 0._pr; pp = 0._pr;
       p1 = 0._pr; p2 = 0._pr; p3 = 0._pr
       eps = 1.5e-16_pr
       m = (npts + 1)/2

       do i=1,m + 1

          t = cos(pi*(dble(i) - 0.25_pr)/(dble(npts+1) + 0.5_pr) )
          t1 = 1._pr

          do while( (abs(t - t1) ) >= eps)

             p1 = 1._pr ; p2 = 0._pr

             do j=1,npts + 1
                p3 = p2; p2 = p1
                p1 = ( (2._pr*dble(j) - 1._pr)*t*p2 - dble(j-1)*p3)/(dble(j) )
             end do

             pp = dble(npts+1)*(t*p1 - p2)/(t*t - 1._pr)
             t1 = t
             t = t1 - p1/pp

          end do

          x(i - 1) = - t;
          x(npts + 1 - i) = t
          w(i - 1) = 2._pr/( (1._pr - t*t)*pp*pp)
          w(npts + 1 - i) = w(i - 1)

       enddo

          do i=0, npts
             x(i) = x(i)*(b - a)/2._pr + (b + a)/2._pr
             w(i) = w(i)*(b - a)/2._pr
          enddo

       return
       end subroutine gauss

       end module gaussmod
