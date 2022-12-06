program p2e
use precision, only: pr
implicit none

! Este programa calcula el exponente de Lyapunov en función de r

  integer                     :: i, j
  real(pr)                    :: r, l, x

  ! Abrimos el documento de texto
  open(25,file='exl.dat')

  ! dividimos el intervalo de r entre 0 y 4 en 40000 secciones idénticas
  do i = 0, 40000
    ! Inicializamos el parámetro variable r, el exponente de Lyapunov l y la posición inicial x(0)
    r = real(i,pr)/10000._pr
    l = 0._pr
    x = 0.6_pr
    ! Hacemos evolucionar x en el tiempo. No usamos vectores para ahorrar memoria
    do j = 1, 300
      x = r*x*(1._pr-x)
    end do
    ! descartamos el transitorio de 300 pasos y consideramos los siguientes 5000
    do j = 1, 5000
      x = r*x*(1._pr-x)
      ! Sumamos los términos de l
      l = l + log(abs(r-r*x*2._pr))/5000._pr
    end do
    ! Anotamos el valor obtenido para el exponente de Lyapunov en función de r
    write(25,*) r, l
  end do

  close(25)

end program p2e
