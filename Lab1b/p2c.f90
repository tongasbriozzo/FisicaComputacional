program p2c
use precision, only: pr
implicit none

! Este programa da una función y(x) para graficar el histograma de x_i con r = 4

  integer                      :: i, j, n
  real(pr)                     :: r, m, d
  real(pr), dimension(1:10000) :: x
  real(pr), dimension(1:100)   :: y

  ! inicializamos las variables
  r    = 4._pr
  x(1) = 0.6_pr

  ! Evolicionamos el sistema en el tiempo
  do i = 1, 299
    x(i+1) = r*x(i)*(1._pr-x(i))
  end do
  ! Descartamos el transitorio de 300 pasos y calculamos los siguientes 10000
  x(1) = x(300)
  do i = 1, 9999
    x(i+1) = r*x(i)*(1._pr-x(i))
  end do

  ! Inicializamos los valores de la función y
  do i = 1, 100
    y(i) = 0._pr
  end do

  ! Asignamos los valores a y(i)
  do i = 1, 10000
    m = x(i)*100._pr
    do j = 1, 100
      d = m-real(j,pr)
      if (d<=0) then
        y(j) = y(j) + 1._pr
        exit
      end if
    end do
  end do

  ! abrimos el documento de texto y anotamos los valores de y(i)
  open(23,file='his.dat')
  do i = 1, 100
    write(23,*) real(i-1,pr)/100._pr, 0 ! El cero al inicio de cada bin es para unir el histograma con líneas
    ! Para darle continuidad al gráfico, dividimos cada bin en 100 segmentos identicos, cada uno de los cuales corresponde a un mismo y(i)
    do j = 0, 100
      write(23,*) real(100*(i-1)+j,pr)/10000._pr, y(i)/10000._pr  ! Notese que y(i) está normalizada, es decir, indica probabilidad
    end do
    write(23,*) real(i,pr)/100._pr, 0 ! El cero al final de cada bin es para unir el histograma con lineas
  end do
  close(23)

end program p2c
