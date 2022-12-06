program p3
implicit none

  integer			 :: i, j, k, n
  integer, parameter ::	sp=4, dp=8
  real(kind=dp)		 ::	In, T, S, h, x, a, b, w1, w2, w3, w4, w5

  1     format(A10,     4x,  A8,    16x,  A7)
  2     format(I12,     2x,  E22.16, 2x,  E22.16)

  In = exp(1._dp)-(1._dp)
  a  = 0._dp
  b  = 1._dp
  
  open(31,file='its')
  write(31,1) 'Divisiones', 'Trapecio', 'Simpson'

  do i = 2, 29, 1

    n = 2**i
	h = (b-a)/real(n)

	w1 = h
	w2 = h/(2._dp)
	w3 = h*(4._dp)/(3._dp)
	w4 = h*(2._dp)/(3._dp) 
	w5 = h/(3._dp)

	T = 0._dp
	S = 0._dp

	do j = 1, n-1
	  x = real(j)*h
	  T = T + w1*exp(x)
	end do
	T = T +	w2*(exp(a)+exp(b))

	do k = 1, n-1, 2
	  x = real(k)*h
	  S = S + w3*exp(x)
	end do 
	do k = 2, n-2, 2
	  x = real(k)*h
	  S = S + w4*exp(x)
	end do
	S = S + w5*(exp(a)+exp(b))

	write(31,2)	n,      abs((T-In)/In),   abs((S-In)/In)
	  
  end do
	  
  close(31)

end program p3