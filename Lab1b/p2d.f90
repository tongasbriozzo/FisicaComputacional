program p2d
use precision, only: pr
implicit none

! Este programa calcula el diagrama de órbitas de la ec. logística para distintos intervalos de r

  integer                     :: i, j
  real(pr)                    :: r
  real(pr), dimension(1:600)  :: x
  real(pr), dimension(1:1500) :: y

  open(241,file='do1.dat')  ! Aquí, obtenemos los datos para la primer parte del problema, con 3.4<r<4

  do i = 0, 6000 ! Dividimos el intervalo en 6000 secciones idénticas

    ! Inicializamos x y r
    r    = real(34000+i,pr)/10000._pr
    x(1) = 0.6_pr

    ! Hacemos evolucionar en el tiempo a la ec. logística
    do j = 1, 299
      x(j+1) = r*x(j)*(1._pr-x(j))
    end do
    ! Descartamos el transitorio de 300 pasos y consideramos los siguientes 600
    x(1) = x(300)
    do j = 1, 599
      x(j+1) = r*x(j)*(1._pr-x(j))
    end do

    ! Anotamos los x(i) obtenidos para cada r
    do j = 1, 600
      write(241,*) r, x(j)
    end do

  end do

  close(241)

  ! Ahora fino

  open(242,file='do2.dat')  ! Aquì, obtenemos los datos para la segunda parte del problema, con 3.847<r<3.8568

  do i = 0, 9800  ! Dividimos el intervalo en 9800 secciones idénticas

    ! Inicializamos r e y
    r    = real(3847000+i,pr)/1000000._pr
    y(1) = 0.6_pr

    ! Hacemos evolucionar en el tiempo la ec. logística
    do j = 1, 299
      y(j+1) = r*y(j)*(1._pr-y(j))
    end do
    ! Descartamos el transitorio de 300 pasos y consideramos los siguientes 1500
    y(1) = y(300)
    do j = 1, 1499
      y(j+1) = r*y(j)*(1._pr-y(j))
    end do

    ! Anotamos los y(i) obtenidos para cada r
    do j = 1, 1500
      write(242,*) r, y(j)
    end do

  end do

  close(242)

  ! Terminó el problema, lo que sigue es por diversión

  open(243,file='do3.dat')  ! Aquí, anotamos los datos para r en el intervalo de 1 a 4
  open(244,file='do4.dat')  ! Aquí, anotaremos y(i+1) en función de y(i)

  do i = 0, 10000  ! dividimos el intervalo en 10000 secciones idénticas

    ! Inicializamos r e y
    r    = 1._pr + real(3*i,pr)/10000._pr
    y(1) = 0.6_pr

    ! Hacemos evolucionar la ec. logística en el tiempo
    do j = 1, 299
      y(j+1) = r*y(j)*(1._pr-y(j))
    end do
    ! Descartamos el transitorio de 300 pasos y consideramos los siguientes 1500
    y(1) = y(300)
    do j = 1, 1499
      y(j+1) = r*y(j)*(1._pr-y(j))
      if (i==10000) then  ! Aquì, para r=4, anotamos y(i+1) en función de y(i)
        write(244,*) y(j), y(j+1)
      end if
    end do

    ! Anotamos los y(i) obtenidos para cada r
    do j = 1, 1500
      write(243,*) r, y(j)
    end do

  end do

  close(243)
  close(244)

end program p2d
