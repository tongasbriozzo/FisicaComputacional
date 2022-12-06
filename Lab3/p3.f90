! ______________________________________________________________________________
!
! MONTE CALRLO
!
! Este programa emplea el método de Monte Carlo para realizar la integral de
! x^expo entre los límites de integración a y b. Para esto generamos numeros
! pseudoaleatorios mediante los métodos ran0, ran2, Marsaglia y Zaman y Mersenne
! Twister. Los resultados numericos son comparados con la sulución analítica
! para obtener los errores relativos. El procedimiento se repite para distintas
! cantidades de puntos a evaluar.
!
! Gaston Briozzo, FaMAFyC, UNC, 06/05/2021
! ______________________________________________________________________________

program p3
use precision, only: pr
use randommod, only: ran0, ran2
use mzranmod,  only: rmzran
use mtmod,     only: grnd
implicit none

  integer                  :: ni, i, j, k, nn
  integer                  :: expo, sem1, sem2
  real(pr)                 :: a, b, Ie
  real(pr), dimension(1:4) :: In0, In1, In2, In23, In3
  real(pr), dimension(1:4) :: Er0, Er1, Er2, Er23, Er3
  real(pr), dimension(1:4) :: Vz0, Vz1, Vz2, Vz23, Vz3  ! Aclaración: las variables Vz y sg se emplean tando para calcular
  real(pr), dimension(1:4) :: sg0, sg1, sg2, sg23, sg3  ! las componentes de la varianza como para representar su valor final
  real(pr), dimension(1:4) :: x

! Primero, definimos los formatos que emplearemos en los documentos de texto
  1 format (A1, 9x, A5,    20x, A5,    20x, A5,    20x, A5)
  2 format (I7, 3x, E23.16, 2x, E23.16, 2x, E23.16, 2x, E23.16)
  3 format (A1, 9x, A5,    20x, A5,    20x, A5,    20x, A5,    20x, A5,    20x, A5,    20x, A5,    20x, A5)
  4 format (I7, 3x, E23.16, 2x, E23.16, 2x, E23.16, 2x, E23.16, 2x, E23.16, 2x, E23.16, 2x, E23.16, 2x, E23.16)

! Inicializamos los parámetros
! Límites de integración
  a = 0._pr
  b = 1._pr
! Exponente
  expo = 3
! Integral analítica
  Ie = 1._pr/(real(expo,pr)+1._pr)
! Semillas
  sem1 = 1073741827  ! Es primo
  sem2 = 2073741823  ! Es primo

! Abrimos los documentos de texto y describimos los datos
  open(330, file='is0.d')
  open(331, file='is1.d')
  open(332, file='is2.d')
  open(3323,file='is23.d')
  open(333, file='is3.d')
  open(340, file='va0.d')
  open(341, file='va1.d')
  open(342, file='va2.d')
  open(3423,file='va23.d')
  open(343, file='va3.d')
  write(330,1)  'N', 'Er_r0', 'Er_r2', 'Er_mz', 'Er_mt'
  write(331,1)  'N', 'Er_r0', 'Er_r2', 'Er_mz', 'Er_mt'
  write(332,1)  'N', 'Er_r0', 'Er_r2', 'Er_mz', 'Er_mt'
  write(3323,1) 'N', 'Er_r0', 'Er_r2', 'Er_mz', 'Er_mt'
  write(333,1)  'N', 'Er_r0', 'Er_r2', 'Er_mz', 'Er_mt'
  write(340,3)  'N', 'Va_r0', 'sg_r0', 'Va_r2', 'sg_r2', 'Va_mz', 'sg_mz', 'Va_mt', 'sg_mt'
  write(341,3)  'N', 'Va_r0', 'sg_r0', 'Va_r2', 'sg_r2', 'Va_mz', 'sg_mz', 'Va_mt', 'sg_mt'
  write(342,3)  'N', 'Va_r0', 'sg_r0', 'Va_r2', 'sg_r2', 'Va_mz', 'sg_mz', 'Va_mt', 'sg_mt'
  write(3423,3) 'N', 'Va_r0', 'sg_r0', 'Va_r2', 'sg_r2', 'Va_mz', 'sg_mz', 'Va_mt', 'sg_mt'
  write(343,3)  'N', 'Va_r0', 'sg_r0', 'Va_r2', 'sg_r2', 'Va_mz', 'sg_mz', 'Va_mt', 'sg_mt'

! ______________________________________________________________________________
! Primero probamos con números aleatorios de distribución homogenea entre 0 y 1
! Luego, tomamos la distribución p(x)=(k+1)x^k, con k=1,2
! ______________________________________________________________________________

! Realizamos 10**6 integrales, variando nn (numero de puntos a evaluar)
  do ni = 0, 5
  do i = 1, 10
    nn = i*(10**ni)                 ! Definimos cuantos puntos evaluaremos
!   Integrales numéricas
    In0  = 0._pr
    In1  = 0._pr
    In2  = 0._pr
    In23 = 0._pr
    In3  = 0._pr
!   Varianzas
    Vz0  = 0._pr
    Vz1  = 0._pr
    Vz2  = 0._pr
    Vz23 = 0._pr
    Vz3  = 0._pr
!   Sigma
    sg0  = 0._pr
    sg1  = 0._pr
    sg2  = 0._pr
    sg23 = 0._pr
    sg3  = 0._pr
    do j = 1, nn              ! Creamos nn números aleatórios por método
      x(1) = ran0(sem1)
      x(2) = ran2(sem2)
      x(3) = rmzran()
      x(4) = grnd()
      do k = 1, 4
!       Evaluamos las primitivas y las integramos
        In0(k)  = In0(k)  + f(x(k),expo)
        In1(k)  = In1(k)  + f(y(x(k),1._pr), expo)/p(y(x(k),1._pr), 1._pr)
        In2(k)  = In2(k)  + f(y(x(k),2._pr), expo)/p(y(x(k),2._pr), 2._pr)
        In23(k) = In23(k) + f(y(x(k),2.3_pr),expo)/p(y(x(k),2.3_pr),2.3_pr)
        In3(k)  = In3(k)  + f(y(x(k),3._pr), expo)/p(y(x(k),3._pr), 3._pr)
!       Término de la varianza, cuadrado de la suma
        Vz0(k)  = Vz0(k)  + f(x(k),expo)
        Vz1(k)  = Vz1(k)  + f(y(x(k),1._pr), expo)/p(y(x(k),1._pr), 1._pr)
        Vz2(k)  = Vz2(k)  + f(y(x(k),2._pr), expo)/p(y(x(k),2._pr), 2._pr)
        Vz23(k) = Vz23(k) + f(y(x(k),2.3_pr),expo)/p(y(x(k),2.3_pr),2.3_pr)
        Vz3(k)  = Vz3(k)  + f(y(x(k),3._pr), expo)/p(y(x(k),3._pr), 3._pr)
!       Término de la varianza, suma de cuadrados
        sg0(k)  = sg0(k)  + (f(x(k),expo))**2
        sg1(k)  = sg1(k)  + (f(y(x(k),1._pr), expo)/p(y(x(k),1._pr), 1._pr))**2
        sg2(k)  = sg2(k)  + (f(y(x(k),2._pr), expo)/p(y(x(k),2._pr), 2._pr))**2
        sg23(k) = sg23(k) + (f(y(x(k),2.3_pr),expo)/p(y(x(k),2.3_pr),2.3_pr))**2
        sg3(k)  = sg3(k)  + (f(y(x(k),3._pr), expo)/p(y(x(k),3._pr), 3._pr))**2
      end do
    end do
!   Obtenemos las integrales numéricas
    In0(:)  = (b-a)*In0(:)/real(nn,pr)
    In1(:)  = (b-a)*In1(:)/real(nn,pr)
    In2(:)  = (b-a)*In2(:)/real(nn,pr)
    In23(:) = (b-a)*In23(:)/real(nn,pr)
    In3(:)  = (b-a)*In3(:)/real(nn,pr)
!   Calculamos el error relativo
    Er0(:)  = abs(In0(:) -Ie)/Ie
    Er1(:)  = abs(In1(:) -Ie)/Ie
    Er2(:)  = abs(In2(:) -Ie)/Ie
    Er23(:) = abs(In23(:)-Ie)/Ie
    Er3(:)  = abs(In3(:) -Ie)/Ie
!   Calculamos la varianza
    Vz0(:)  = -Vz0(:)*Vz0(:)/real(nn,pr)/real(nn,pr)   + sg0(:)/real(nn,pr)
    Vz1(:)  = -Vz1(:)*Vz1(:)/real(nn,pr)/real(nn,pr)   + sg1(:)/real(nn,pr)
    Vz2(:)  = -Vz2(:)*Vz2(:)/real(nn,pr)/real(nn,pr)   + sg2(:)/real(nn,pr)
    Vz23(:) = -Vz23(:)*Vz23(:)/real(nn,pr)/real(nn,pr) + sg23(:)/real(nn,pr)
    Vz3(:)  = -Vz3(:)*Vz3(:)/real(nn,pr)/real(nn,pr)   + sg3(:)/real(nn,pr)
!   Calsulamos los sigmas
    sg0(:)  = sqrt(Vz0(:))
    sg1(:)  = sqrt(Vz1(:))
    sg2(:)  = sqrt(Vz2(:))
    sg23(:) = sqrt(Vz23(:))
    sg3(:)  = sqrt(Vz3(:))
!   Anotamos los resultados
!   Anotamos los errores relativos
    write(330,2)  nn, Er0(1),  Er0(2),  Er0(3),  Er0(4)
    write(331,2)  nn, Er1(1),  Er1(2),  Er1(3),  Er1(4)
    write(332,2)  nn, Er2(1),  Er2(2),  Er2(3),  Er2(4)
    write(3323,2) nn, Er23(1), Er23(2), Er23(3), Er23(4)
    write(333,2)  nn, Er3(1),  Er3(2),  Er3(3),  Er3(4)
!   Anotamos las varianzas y sigmas
    write(340,4)  nn, Vz0(1),  sg0(1),  Vz0(2),  sg0(2),  Vz0(3),  sg0(3),  Vz0(4),  sg0(4)
    write(341,4)  nn, Vz1(1),  sg1(1),  Vz1(2),  sg1(2),  Vz1(3),  sg1(3),  Vz1(4),  sg1(4)
    write(342,4)  nn, Vz2(1),  sg2(1),  Vz2(2),  sg2(2),  Vz2(3),  sg2(3),  Vz2(4),  sg2(4)
    write(3423,4) nn, Vz23(1), sg23(1), Vz23(2), sg23(2), Vz23(3), sg23(3), Vz23(4), sg23(4)
    write(343,4)  nn, Vz3(1),  sg3(1),  Vz3(2),  sg3(2),  Vz3(3),  sg3(3),  Vz3(4),  sg3(4)
!   Anotamos en terminar para supervizar el funcionamiento del programa
    write(*,2) nn, Er0(1),  Er0(2),  Er0(3),  Er0(4)
    write(*,2) nn, Er1(1),  Er1(2),  Er1(3),  Er1(4)
    write(*,2) nn, Er2(1),  Er2(2),  Er2(3),  Er2(4)
    write(*,2) nn, Er23(1), Er23(2), Er23(3), Er23(4)
    write(*,2) nn, Er3(1),  Er3(2),  Er3(3),  Er3(4)
  end do
  end do

  close(330)
  close(331)
  close(332)
  close(3323)
  close(333)
  close(340)
  close(341)
  close(342)
  close(3423)
  close(343)

CONTAINS

  function f(xx,expo)
  implicit none
    real(pr), intent(in) :: xx
    integer, intent(in)  :: expo
    real(pr)             :: f
    f = xx**expo
  end function f

  function p(xx,kk)
  implicit none
    real(pr), intent(in) :: xx, kk
    real(pr)             :: p
    p = real(kk+1,pr)*xx**kk
  end function p

  function y(x,kk)
  implicit none
    real(pr), intent(in) :: x, kk
    real(pr)             :: y
    y = x**(1._pr/(real(kk,pr)+1._pr))
  end function y

end program p3
